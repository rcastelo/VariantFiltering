## constructor
GenePhylostrataDb <- function(organism, download_url, download_date,
                              source_pub, data_pkgname, data_dirpath) {
  data_serialized_objnames <- "genePhylostrata"
  data_cache <- new.env(hash=TRUE, parent=emptyenv())
  load(file.path(data_dirpath, "humanEnsGenePhylostrata.rda"), envir=data_cache)

  gps <- get("humanEnsGenePhylostrata", envir=data_cache)
  gps <- unique(gps[, c("TaxID", "OldestPhylostratum", "Description")])
  strataNames <- gps[gps$OldestPhylostratum, "Description"]
  strataTaxID <- gps[gps$OldestPhylostratum, "TaxID"]

  new("GenePhylostrataDb", organism=organism,
                           download_url=download_url,
                           download_date=download_date,
                           source_pub=source_pub,
                           data_pkgname=data_pkgname,
                           data_dirpath=data_dirpath,
                           data_serialized_objnames=data_serialized_objnames,
                           strataNames=strataNames,
                           strataTaxID=strataTaxID,
                           .data_cache=data_cache)
}

## accessors

setMethod("organism", "GenePhylostrataDb", function(x) x@organism)

setMethod("genePhylostratum", c("GenePhylostrataDb", "missing"),
          function(object, ids) {
            get("humanEnsGenePhylostrata", envir=object@.data_cache)
          })

setMethod("genePhylostratum", c("GenePhylostrataDb", "character"),
          function(object, ids) {
            n <- length(ids)
            res <- data.frame(EntrezID=rep(NA_character_, n),
                              TaxID=rep(NA_character_, n),
                              OldestPhylostratum=rep(NA_integer_, n),
                              Description=rep(NA_character_, n),
                              stringsAsFactors=FALSE)

            idsNotNA <- ids[!is.na(ids)]
            if (length(idsNotNA) > 0) {
              ## if the input identifers are not Ensembl, Entrez ID is assumed
              if (substr(idsNotNA[1], 1, 4) != "ENSG") {
                eg2ens <- get("entrezToPhylostrataEnsGene", envir=object@.data_cache)
                idsNotNA <- eg2ens[idsNotNA]
              }

              gps <- genePhylostratum(object)
              res[!is.na(ids), ] <- gps[idsNotNA, ]
            }

            res
          })

setMethod("genePhylostrata", "GenePhylostrataDb",
          function(object) {
            df <- data.frame(TaxID=object@strataTaxID, Description=object@strataNames,
                             stringsAsFactors=FALSE)
            rownames(df) <- 1:nrow(df)
            df
          })


## show method
setMethod("show", "GenePhylostrataDb",
          function(object) {
            cat(class(object), "object for", organism(object), "\n")
            cat("  Source URL: ", object@download_url, "\n")
            cat("  Download date: ", object@download_date, "\n")
            cat("  Source publication: ", object@source_pub, "\n")
          })
