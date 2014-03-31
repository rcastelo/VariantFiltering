.isPkgLoaded <- function(name) {
   (paste("package:", name, sep="") %in% search()) ## || (name %in% loadedNamespaces()) <- problematic with cleanEx()
}


setMethod("VariantFilteringParam", signature(vcfFilenames="character"),
          function(vcfFilenames, pedFilename=character(),
                   orgdb="org.Hs.eg.db",
                   txdb="TxDb.Hsapiens.UCSC.hg19.knownGene",
                   snpdb="SNPlocs.Hsapiens.dbSNP.20120608",
                   radicalAAchangeFilename=file.path(system.file("extdata", package="VariantFiltering"),
                                                     "AA_chemical_properties_HanadaGojoboriLi2006.tsv"),
                   allTranscripts=FALSE,
                   otherAnnotations=c("MafDb.ESP6500SI.V2.SSA137.dbSNP138",
                                      "MafDb.ALL.wgs.phase1.release.v3.20101123",
                                      "PolyPhen.Hsapiens.dbSNP131",
                                      "SIFT.Hsapiens.dbSNP137",
                                      "phastCons100way.UCSC.hg19",
                                      "humanGenesPhylostrata"),
                   filterTag=NA_character_) {

            ## store call to reproducing it later
            callobj <- match.call()
            callstr <- gsub(".local", "VariantFilteringParam", deparse(callobj))
            callstr <- gsub("= vcfFilenames", sprintf("= c(%s)", paste(sprintf("\"%s\"", vcfFilenames), collapse=", ")), callstr)

            ## check if input VCF, PED and radical AA change files exist
            tryCatch({
              .io_check_exists(c(vcfFilenames, pedFilename, radicalAAchangeFilename))
            }, error=function(err) {
                 stop(conditionMessage(err), call.=FALSE)
            })

            maskGz <- grepl("vcf.gz$", vcfFilenames)

            ## files that are no bgzipped should get bgzipped for the indexTabis() call below
            if (any(!maskGz)) {
              maskVcf <- grepl("vcf$", vcfFilenames[!maskGz])
              if (!all(maskVcf))
                stop(sprintf("file(s) %s have no .vcf extension.", paste(vcfFilenames[!maskGz][!maskVcf], collapse=", ")))
              sapply(vcfFilenames[!maskGz], function(f) {
                                              message(sprintf("Tabix compressing with bgzip %s", f))
                                              bgzip(f, overwrite=TRUE)
                                            })
              vcfFilenames[!maskGz] <- paste0(vcfFilenames[!maskGz], ".gz")
            }

            maskTbi <- .io_exists_mask(paste0(vcfFilenames, ".tbi"))
            ## files with no tabix index should get built their tabix index
            if (any(!maskTbi)) {
              for (f in vcfFilenames) {
                tryCatch({
                  message(sprintf("Tabix indexing %s", f))
                  indexTabix(f, format="vcf")
                }, error=function(err) {
                  stop(conditionMessage(err), call.=TRUE)
                })
              }
            }

            tfl <- list()
            for (i in seq(along=vcfFilenames))
              tfl[[i]] <- TabixFile(vcfFilenames[i])
            tfl <- do.call(TabixFileList, tfl)

            ## read radical amino acid change matrix
            radicalAAchangeMatrix <- readAAradicalChangeMatrix(radicalAAchangeFilename)

            ## check that the given annotation packages are installed and can be loaded,
            ## load them and save the annotation object into the argument
            if (!is.character(orgdb) && length(orgdb) != 1)
              stop("argument 'orgdb' should contain the name of an 'OrgDb' gene-centric annotation package.")
            else {
              if (!orgdb %in% installed.packages()[, "Package"])
                stop(sprintf("please install the Bioconductor package %s.", orgdb))
              if (!.isPkgLoaded(orgdb)) {
                message("Loading gene-centric annotation package ", orgdb)
                if (!suppressPackageStartupMessages(require(orgdb, character.only=TRUE)))
                  stop(sprintf("package %s could not be loaded.", orgdb))
              }
              tryCatch({
                orgdb <- get(orgdb)
              }, error=function(err) {
                stop(sprintf("The gene annotation package %s should automatically load an 'OrgDb' object with the same name as the package."))
              })

              if (!is(orgdb, "OrgDb"))
                stop(sprintf("The object loaded with name %s is not an 'OrgDb' object.", orgdb))
            }

            if (!is.character(txdb) && length(txdb) != 1)
              stop("argument 'txdb' should contain the name of a 'TxDb' transcript-centric annotation package.")
            else {
              if (!txdb %in% installed.packages()[, "Package"])
                stop(sprintf("please install the Bioconductor package %s.", txdb))
              if (!.isPkgLoaded(txdb)) {
                message("Loading transcript-centric annotation package ", txdb)
                if (!suppressPackageStartupMessages(require(txdb, character.only=TRUE)))
                  stop(sprintf("package %s could not be loaded.", txdb))
              }
              tryCatch({
                txdb <- get(txdb)
              }, error=function(err) {
                stop(sprintf("The gene annotation package %s should automatically load a 'TranscriptDb' object with the same name as the package."))
              })
              if (!is(txdb, "TranscriptDb"))
                stop(sprintf("The object loaded with name %s is not a 'TranscriptDb' object.", txdb))
            }

            if (!is.character(snpdb) && length(snpdb) != 1)
              stop("argument 'snpdb' should contain the name of a 'SNPlocs' SNP-centric annotation package.")
            else {
              if (!snpdb %in% installed.packages()[, "Package"])
                stop(sprintf("please install the Bioconductor package %s.", snpdb))
              if (!.isPkgLoaded(snpdb)) {
                message("Loading SNP-centric annotation package ", snpdb)
                if (!suppressPackageStartupMessages(require(snpdb, character.only=TRUE)))
                  stop(sprintf("package %s could not be loaded.", snpdb))
              }
              tryCatch({
                snpdb <- get(snpdb)
              }, error=function(err) {
                stop(sprintf("The gene annotation package %s should automatically load a 'SNPlocs' object with the same name as the package."))
              })
              if (!is(snpdb, "SNPlocs"))
                stop(sprintf("The object loaded with name %s is not a 'SNPlocs' object.", snpdb))
            }

            if (!is.logical(allTranscripts))
              stop("argument 'allTranscripts' should be logical.")

            if (!is.character(otherAnnotations))
              stop("argument 'otherAnnotations' should contain names of objects or installed packages that can be used by VariantFiltering to annotate genetic variants.")

            otherannotations <- list()
            for (name in otherAnnotations)  {
              if (!exists(name)) {
                if (!name %in% installed.packages()[, "Package"])
                  stop(sprintf("please install the Bioconductor package %s.", name))
                if (!.isPkgLoaded(name)) {
                  message("Loading annotation package ", name)
                  if (!suppressPackageStartupMessages(require(name, character.only=TRUE)))
                    stop(sprintf("package %s could not be loaded.", name))
                }
              } else
                message("Fetching annotation object ", name)
              otherannotations[[name]] <- get(name)
            }

            new("VariantFilteringParam", callObj=callobj, callStr=callstr, vcfFiles=tfl,
                pedFilename=pedFilename, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
                radicalAAchangeFilename=radicalAAchangeFilename, radicalAAchangeMatrix=radicalAAchangeMatrix,
                otherAnnotations=otherannotations, allTranscripts=allTranscripts, filterTag=filterTag)
          })

setMethod("show", signature(object="VariantFilteringParam"),
          function(object) {
            cat("\nVariantFiltering parameter object\n\n")
            cat("  VCF file(s):")
            for (f in path(object$vcfFiles))
              cat(sprintf(" %s", basename(f)))
            cat(sprintf("\n  PED file: %s\n", basename(object$pedFilename)))
            cat(sprintf("  Gene-centric annotation package: %s\n", object$orgdb$packageName))
            cat(sprintf("  Transcript-centric annotation package: %s\n", object$txdb$packageName))
            cat(sprintf("  SNP-centric annotation package: %s (%s %s)\n",
                        object$snpdb@data_pkgname, provider(object$snpdb), releaseName(object$snpdb)))
            cat(sprintf("  Radical/Conservative AA changes file: %s\n", basename(object$radicalAAchangeFilename)))
            cat(sprintf("  Other annotation pkg/obj: %s\n",
                        paste(names(object$otherAnnotations),
                              collapse=",\n                            ")))
            cat(sprintf("  All transcripts: %s\n", object$allTranscripts))
            cat(sprintf("  Filter tag: %s\n\n", object$filterTag))
          })

setMethod("names", signature(x="VariantFilteringParam"),
          function(x) {
            n <- slotNames(x)
            n[-grep(c("callObj", "callStr"), n)]
          })

setMethod("$", signature(x="VariantFilteringParam"),
          function(x, name) {
            if (is.na(match(name, slotNames(x))))
              stop("unknown VariantFilteringParam slot. Use names() to find out which are the valid ones.")

            slot(x, name)
          })


## some utility functions copied from Rsamtools/R/utilities.R

.ppath <- function(tag, filepath) {
  wd <- options('width')[[1]] - nchar(tag) - 6
  if (0L == length(filepath) || nchar(filepath) < wd)
    return(sprintf("%s: %s\n", tag, filepath))
  bname <- basename(filepath)
  wd1 <- wd - nchar(bname)
  dname <- substr(dirname(filepath), 1, wd1)
  sprintf("%s: %s...%s%s\n",
          tag, dname, .Platform$file.sep, bname)
}

.io_check_exists <- function(file) {
  if (!length(file))
    stop("'file' is length(0)")
  idx <- !grepl("^(ftp)|(http)://", file)
  if (!all(sapply(file[idx], file.exists))) {
    msg <- paste(sprintf("'%s'", file[idx]), collapse="\n  ")
    stop("file(s) do not exist:\n  ", msg)
  }
}

.io_exists_mask <- function(file) {
  if (!length(file))
    stop("'file' is length(0)")
  idx <- !grepl("^(ftp)|(http)://", file)
  mask <- sapply(file[idx], file.exists)
  mask
}

.normalizePath <- function(path) {
    idx <- !grepl("^(ftp)|(http)://", path)
    ## expand ~/, but don't chase links (i.e., don't normalizePath())
    path[idx] <- path.expand(path[idx])
    path
}
