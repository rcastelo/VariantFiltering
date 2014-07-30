## constructor
MafDb <- function(conn) {
  .MafDb$new(conn=conn)
}

## extractor of the columns containing MAF values
setMethod("knownVariantsMAFcols", signature(mafdb="MafDb"),
          function(mafdb) {
            conn <- dbConn(mafdb)
            sql <- "PRAGMA table_info(knownVariants)"
            res <- .dbEasyQuery(conn, sql)$name
            res <- res[grep("AF", res)]
            res
          })

## method for fetching variants in a MafDb database object
setMethod("fetchKnownVariantsByID", signature(mafdb="MafDb"),
          function(mafdb, varID) {

            conn <- dbConn(mafdb)
            if (any(is.na(varID))) {
              warning("'NA' variant IDs have been removed")
              varID <- varID[!is.na(varID)]
            }

            sql <- sprintf("SELECT * FROM knownVariants WHERE varID in (%s)",
                           paste0('"', varID, '"', collapse=","))
            res <- .dbEasyQuery(conn, sql)
            idxAFcols <- grep("AF", colnames(res))
            for (i in idxAFcols)
              res[[i]] <- decodeRAW2AF(sapply(res[[i]], charToRaw))

            ## return a data frame with varID values in the same order as they were queried
            df <- data.frame(varID=varID)
            res <- merge(df, res, all.x=TRUE)
            res <- res[match(df$varID, res$varID), ]
            res
          })

setMethod("dbConn", "MafDb",
          function(x) x$conn)

setMethod("keytypes", "MafDb",
          function(x) return("varID"))

setMethod("keys", "MafDb",
          function(x, keytype) {
            if (missing(keytype))
              keytype <- "varID"

            if (keytype != "varID")
              stop("Only 'varID' is a valid keytype.")

            res <- .dbEasyQuery(dbConn(x), "SELECT varID FROM knownVariants")
            unique(as.character(res[!is.na(res)]))
          })

setMethod("columns", "MafDb",
          function(x)
            dbListFields(dbConn(x), "knownVariants"))

setMethod("select", "MafDb",
          function(x, keys, columns, keytype) {
            if (missing(keytype)) keytype <- "varID"
            if (missing(columns)) columns <- columns(x)

            if (any(!(columns %in% columns(x))))
              stop(sprintf("Columns %s do not form part of this annotation package.",
                           paste(columns[!(columns %in% columns(x))], collapse=", ")))

            res <- fetchKnownVariantsByID(x, keys)
            res[, columns]
          })

##
## private functions
##

load_taginfo <- function(mafdb) {
  conn <- dbConn(mafdb)
  sql <- "SELECT value FROM metadata WHERE name='Data source tag'"
  res <- .dbEasyQuery(conn, sql)
  res[[1]]
}

## adapted from dbEasyQuery() in GenomicFeatures/R/utils.R
.dbEasyQuery <- function(conn, SQL, j0=NA) {
  data0 <- dbGetQuery(conn, SQL)

  if (is.na(j0))
    return(data0)

  ## Needed to deal properly with data frame with 0 column ("NULL data
  ## frames with 0 rows") returned by RSQLite when the result of a SELECT
  ## query has 0 row
 if (nrow(data0) == 0L)
   character(0)
  else
    data0[[j0]]
}

## adapted from dbEasyPreparedQuery() in GenomicFeatures/R/utils.R
dbEasyPreparedQuery <- function(conn, SQL, bind.data) {
  ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
  ## when the nb of rows to insert is 0, hence the early bail out.
  if (nrow(bind.data) == 0L)
    return()
  dbBeginTransaction(conn)
  dbGetPreparedQuery(conn, SQL, bind.data)
  dbCommit(conn)
}

## code and decode allele frequencies into a single-byte raw type, where AF values >= 0.01 are rounded to 2 digits.
## AF values < 0.01 are rounded to 3 digits and AF values < 0.001 are rounded to 4 digits. The latter
## (0, 0.0001, 0.0002, ..., 0.0009) are stored as raw byte values 1 to 10, the next precision level
## (0.001, 0.002, ... 0.009) are stored as raw byte values 11 to 19 and the lower precision level
## (0.01, 0.02, ..., 0.98, 0.99, 1.0) are stored as the raw byte values 20 to 119
## because NAs are encoded by the raw byte value 0, and this corresponds by default to the null string
## when raw byte values are coerced into char, which will be necessary when they are stored in the database.
## NAs will be recoded to the highest possible raw byte value of 255
codeAF2RAW <- function(x) {
  maskNAs <- is.na(x)
  z <- x[!maskNAs]
  maskLevel1 <- z < 0.001
  maskLevel2 <- z >= 0.001 & z < 0.01
  maskLevel3 <- z >= 0.01

  z[maskLevel1] <- round(z[maskLevel1], digits=4)
  z[maskLevel2] <- round(z[maskLevel2], digits=3)
  z[maskLevel3] <- round(z[maskLevel3], digits=2)
  
  maskLevel1 <- z < 0.001
  maskLevel2 <- z >= 0.001 & z < 0.01
  maskLevel3 <- z >= 0.01

  z[maskLevel1] <- z[maskLevel1] * 10000 + 1     ## zero is coded as 1 to avoid the 'nul' byte meaning throughout
  z[maskLevel2] <- z[maskLevel2] * 1000 + 9 + 1  ## as for instance when coercing to char to store these bytes
  z[maskLevel3] <- z[maskLevel3] * 100 + 18 + 1  ## in a mysql database, this implies all values are shifted by 1

  x[!maskNAs] <- z
  x[maskNAs] <- 255 ## code NAs as raw byte value 255
  as.raw(x)
}

decodeRAW2AF <- function(x) {
  x <- as.numeric(x)
  maskNAs <- x == 255 ## decode raw byte value 255 as NAs
  z <- x[!maskNAs]

  maskLevel1 <- z < 10
  maskLevel2 <- z > 9 & z < 20
  maskLevel3 <- z > 19

  z[maskLevel1] <- (z[maskLevel1]-1) / 10000
  z[maskLevel2] <- (z[maskLevel2]-9-1) / 1000
  z[maskLevel3] <- (z[maskLevel3]-18-1) / 100

  x[!maskNAs] <- z
  x[maskNAs] <- NA
  x
}

