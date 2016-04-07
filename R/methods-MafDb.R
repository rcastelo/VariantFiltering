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

setMethod("fetchKnownVariantsByID", signature(mafdb="MafDb"),
          function(mafdb, varID) {
            .Defunct("snpid2maf")
          })

## method for fetching variants in a MafDb database object
setMethod("snpid2maf", signature(mafdb="MafDb"),
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

            res <- snpid2maf(x, keys)
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
.dbEasyPreparedQuery <- function(conn, SQL, bind.data) {
  ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
  ## when the nb of rows to insert is 0, hence the early bail out.
  if (nrow(bind.data) == 0L)
    return()
  dbBegin(conn)
  dbGetPreparedQuery(conn, SQL, bind.data)
  dbCommit(conn)
}

## based on the precision requirements of MAF data, particularly from ExAC MAF data,
## we the lowest observed allele frequency value is 8.236e-06, we
## code and decode allele frequencies into a single-byte raw type, where:
##
## AF values < 0.00001 are rounded to 6 digits;
## AF values >= 0.00001 & values < 0.0001 are rounded to 5 digits;
## AF values >= 0.0001 & values < 0.001 are rounded to 4 digits;
## AF values >= 0.001 & values < 0.01 are rounded to 3 digits;
## AF values >= 0.01 & values <= 1 are rounded to 2 digits;
##
## and
##
## AF values >= 0.01 & values <= 1 rounded to 2 digits are stored as raw byte values 1 to 100
## AF values >= 0.001 & values < 0.01 rounded to 3 digits are stored as raw byte values 101 to 109
## AF values >= 0.0001 & values < 0.001 rounded to 4 digits are stored as raw byte values 111 to 119
## AF values >= 0.00001 & values < 0.0001 rounded to 5 digits are stored as raw byte values 121 to 129
## AF values < 0.00001 rounded to 6 digits are stored as raw byte values 130 to 139 (include 0 in the lowest range)
##
## because NAs are encoded by the raw byte value 0, and this corresponds by default to the null string
## when raw byte values are coerced into char, which will be necessary when they are stored in the database,
## NAs will be recoded to the highest possible raw byte value of 255

codeAF2RAW <- function(x) {
  maskNAs <- is.na(x)
  z <- x[!maskNAs]

  maskLevel1 <- z >= 0.01
  maskLevel2 <- z >= 0.001 & z < 0.01
  maskLevel3 <- z >= 0.0001 & z < 0.001
  maskLevel4 <- z >= 0.00001 & z < 0.0001
  maskLevel5 <- z < 0.00001

  z[maskLevel1] <- round(z[maskLevel1], digits=2)
  z[maskLevel2] <- round(z[maskLevel2], digits=3)
  z[maskLevel3] <- round(z[maskLevel3], digits=4)
  z[maskLevel4] <- round(z[maskLevel4], digits=5)
  z[maskLevel5] <- round(z[maskLevel5], digits=6)
  
  ## recompute again masks to deal with rounding 0.0095 to 0.01, etc.
  maskLevel1 <- z >= 0.01
  maskLevel2 <- z >= 0.001 & z < 0.01
  maskLevel3 <- z >= 0.0001 & z < 0.001
  maskLevel4 <- z >= 0.00001 & z < 0.0001
  maskLevel5 <- z < 0.00001

  ## code AFs into integer numbers
  z[maskLevel1] <- z[maskLevel1] * 100
  z[maskLevel2] <- z[maskLevel2] * 1000    + 100
  z[maskLevel3] <- z[maskLevel3] * 10000   + 110
  z[maskLevel4] <- z[maskLevel4] * 100000  + 120
  z[maskLevel5] <- z[maskLevel5] * 1000000 + 130 ## zero is coded here as 130

  x[!maskNAs] <- z
  x[maskNAs] <- 255 ## code NAs as raw byte value 255
  as.raw(x)
}

decodeRAW2AF <- function(x) {
  x <- as.numeric(x)
  maskNAs <- x == 255 ## decode raw byte value 255 as NAs
  z <- x[!maskNAs]

  maskLevel1 <- z <= 100
  maskLevel2 <- z > 100 & z <= 110 ## use <= as temporary fix for the rounding 0.0095 to 0.01 case
  maskLevel3 <- z > 110 & z <= 120 ## use <= as temporary fix for the rounding 0.0095 to 0.01 case
  maskLevel4 <- z > 120 & z < 130
  maskLevel5 <- z >= 130

  z[maskLevel1] <- z[maskLevel1] / 100
  z[maskLevel2] <- (z[maskLevel2]-100) / 1000
  z[maskLevel3] <- (z[maskLevel3]-110) / 10000
  z[maskLevel4] <- (z[maskLevel4]-120) / 100000
  z[maskLevel5] <- (z[maskLevel5]-130) / 1000000

  x[!maskNAs] <- z
  x[maskNAs] <- NA
  x
}
