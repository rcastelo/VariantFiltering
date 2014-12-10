setMethod("VariantFilteringParam", signature(vcfFilenames="character"),
          function(vcfFilenames, pedFilename=character(),
                   bsgenome="BSgenome.Hsapiens.UCSC.hg19",
                   orgdb="org.Hs.eg.db",
                   txdb="TxDb.Hsapiens.UCSC.hg19.knownGene",
                   snpdb="SNPlocs.Hsapiens.dbSNP.20120608",
                   spliceSiteMatricesFilenames=c(file.path(system.file("extdata", package="VariantFiltering"),
                                                           "hsap.donors.hcmc10_15_1.ibn"),
                                                 file.path(system.file("extdata", package="VariantFiltering"),
                                                           "hsap.acceptors.hcmc10_15_1.ibn")),
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

            ## check if input VCF, PED, splice site and radical AA change files exist
            tryCatch({
              .io_check_exists(c(vcfFilenames, pedFilename, spliceSiteMatricesFilenames, radicalAAchangeFilename))
            }, error=function(err) {
                 stop(conditionMessage(err), call.=FALSE)
            })

            maskGz <- grepl("vcf.bgz$", vcfFilenames)

            ## files that are no bgzipped should get bgzipped for the indexTabix() call below
            if (any(!maskGz)) {
              maskVcf <- grepl("vcf$", vcfFilenames[!maskGz])
              if (!all(maskVcf))
                stop(sprintf("file(s) %s have no .vcf extension.", paste(vcfFilenames[!maskGz][!maskVcf], collapse=", ")))
              sapply(vcfFilenames[!maskGz], function(f) {
                                              message(sprintf("Tabix compressing with bgzip %s", f))
                                              bgzip(f, overwrite=TRUE)
                                            })
              vcfFilenames[!maskGz] <- paste0(vcfFilenames[!maskGz], ".bgz")
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

            sampleNames <- character()
            seqinfos <- vector(mode="list", length=length(vcfFilenames))
            tfl <- list()
            for (i in seq(along=vcfFilenames)) {
              hd <- scanVcfHeader(vcfFilenames[i])
              seqinfos[[i]] <- seqinfo(hd)
              sampleNames <- c(sampleNames, samples(hd))
              tfl[[i]] <- TabixFile(vcfFilenames[i])
            }
            tfl <- do.call(TabixFileList, tfl)

            ## check that the genome information is identical for all VCFs with genome information
            wh <- which(sapply(seqinfos, length) > 0)
            for (i in wh)
              if (!identical(seqinfos[[i]], seqinfos[[wh[1]]]))
                stop("Genome version for VCF file %s is different from VCF file %s\n",
                     basename(vcfFilenames[i]), basename(vcfFilenames[wh[1]]))

            ## those VCFs without genome information take the genome information from
            ## VCFs with genome information (if any)
            wh0 <- which(sapply(seqinfos, length) == 0)
            if (length(wh0) > 0 && length(wh) > 0) {
              seqinfos[wh0] <- seqinfos[wh[1]]
              warning(sprintf("Genome information missing in VCF file(s) %s is taken from VCF file %s.",
                              paste(sapply(vcfFilenames[wh0], basename), collapse=", "),
                              basename(vcfFilenames[wh[1]])))
            }

            ## read splice site matrices
            spliceSiteMatrices <- list()
            if (any(!is.na(spliceSiteMatricesFilenames))) {
              if (class(spliceSiteMatricesFilenames) != "character" || length(spliceSiteMatricesFilenames) != 2)
                stop("'spliceSiteMatricesFilenames' should be a character vector with two filenames of donor and acceptor splice site matrices, respectively.")
              spliceSiteMatrices <- list(wmDonorSites=readWm(spliceSiteMatricesFilenames[1]),
                                         wmAcceptorSites=readWm(spliceSiteMatricesFilenames[2]))
            }

            ## read radical amino acid change matrix
            radicalAAchangeMatrix <- readAAradicalChangeMatrix(radicalAAchangeFilename)

            ## check that the given annotation packages are installed and can be loaded,
            ## load them and save the annotation object into the argument

            bsgenome <- .loadAnnotationPackageObject(bsgenome, "begenome", "BSgenome")

            orgdb <- .loadAnnotationPackageObject(orgdb, "orgdb", "OrgDb")

            txdb <- .loadAnnotationPackageObject(txdb, "txdb", "TxDb")

            ## when no VCF has genome information, this information is taken from the TxDb package
            wh0 <- which(sapply(seqinfos, length) == 0)
            if (length(wh0) > 0) {
              seqinfos[wh0] <- seqinfo(txdb)
              warning(sprintf("No genome information available from any VCF file. This information will be taken from the transcript-centric package %s, thus assuming a genome version %s with %s chromosome nomenclature\n",
                      txdb$packageName, unique(genome(txdb)), seqlevelsStyle(txdb)))
            }

            snpdb <- .loadAnnotationPackageObject(snpdb, "snpdb", "SNPlocs")

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

            new("VariantFilteringParam", callObj=callobj, callStr=callstr, vcfFiles=tfl, seqInfos=seqinfos,
                sampleNames=sampleNames, pedFilename=pedFilename, bsgenome=bsgenome, orgdb=orgdb, txdb=txdb,
                snpdb=snpdb, spliceSiteMatricesFilenames=spliceSiteMatricesFilenames,
                spliceSiteMatrices=spliceSiteMatrices, radicalAAchangeFilename=radicalAAchangeFilename,
                radicalAAchangeMatrix=radicalAAchangeMatrix, otherAnnotations=otherannotations,
                allTranscripts=allTranscripts, filterTag=filterTag)
          })

setMethod("show", signature(object="VariantFilteringParam"),
          function(object) {
            cat("\nVariantFiltering parameter object\n\n")
            cat("  VCF file(s):")
            for (f in path(object$vcfFiles))
              cat(sprintf(" %s", basename(f)))
            sampleNames <- ifelse(length(object$sampleNames) <= 4, paste(object$sampleNames, collapse=", "),
                                  paste(paste(head(object$sampleNames, n=4), collapse=", "), "...", sep=", "))
            cat(sprintf("\n  Genome version(s):"))
            for (i in seq(along=object$seqInfos))
              if (length(object$seqInfos[[i]]) > 0)
                cat(sprintf(" %s(%s)", paste(unique(genome(object$seqInfos[[i]])), collapse=","),
                            seqlevelsStyle(object$seqInfos[[i]])))
              else
                cat(" NA")

            cat(sprintf("\n  Number of individuals: %d (%s)\n", length(object$sampleNames), sampleNames))
            if (length(object$pedFilename) > 0)
              cat(sprintf("  PED file: %s\n", basename(object$pedFilename)))
            cat(sprintf("  Genome-centric annotation package: %s (%s %s %s)\n",
                        object$bsgenome@pkgname, provider(object$bsgenome),
                        providerVersion(object$bsgenome), releaseName(object$bsgenome)))
            cat(sprintf("  SNP-centric annotation package: %s (%s %s)\n",
                        object$snpdb@data_pkgname, provider(object$snpdb), releaseName(object$snpdb)))
            if (length(object$txdb$packageName) > 0)
              cat(sprintf("  Transcript-centric annotation package: %s\n", object$txdb$packageName))
            else
              cat(sprintf("  Transcript-centric annotation table: %s\n", metadata(object$txdb)[grep("Table", metadata(object$txdb)$name), "value"]))
            cat(sprintf("  Gene-centric annotation package: %s\n", object$orgdb$packageName))
            cat(sprintf("  Splice site matrices: %s\n", paste(basename(object$spliceSiteMatricesFilenames), collapse=", ")))
            cat(sprintf("  Radical/Conservative AA changes: %s\n", basename(object$radicalAAchangeFilename)))
            cat(sprintf("  Other annotation pkg/obj: %s\n",
                        paste(names(object$otherAnnotations),
                              collapse=",\n                            ")))
            cat(sprintf("  All transcripts: %s\n", object$allTranscripts))
            cat(sprintf("  Filter tag: %s\n\n", object$filterTag))
          })

setMethod("names", signature(x="VariantFilteringParam"),
          function(x) {
            n <- slotNames(x)
            n[!n %in% c("callObj", "callStr")]
          })

setMethod("$", signature(x="VariantFilteringParam"),
          function(x, name) {
            if (is.na(match(name, slotNames(x))))
              stop("unknown VariantFilteringParam slot. Use names() to find out which are the valid ones.")

            slot(x, name)
          })


## private functions

.isPkgLoaded <- function(name) {
   (paste("package:", name, sep="") %in% search()) ## || (name %in% loadedNamespaces()) <- problematic with cleanEx()
}


.loadAnnotationPackageObject <- function(pkgName, argName, pkgType) {

  callobj <- match.call()
  annotObj <- NULL

  if (is.character(pkgName)) {
    if (!pkgName %in% installed.packages()[, "Package"])
      stop(sprintf("please install the Bioconductor package %s.", pkgName))
    if (!.isPkgLoaded(pkgName)) {
      message("Loading ", pkgType, " annotation package ", pkgName)
      if (!suppressPackageStartupMessages(require(pkgName, character.only=TRUE)))
        stop(sprintf("package %s could not be loaded.", pkgName))
    }
    tryCatch({
      annotObj <- get(pkgName)
    }, error=function(err) {
      stop(sprintf("The gene annotation package %s should automatically load an %s object with the same name as the package.", pkgName, pkgType))
    })
  } else if (class(pkgName) != pkgType)
    stop(sprintf("argument '%s' should either contain the name of an '%s' annotation package or be an '%s' annotation object itself.",
         argName, pkgType, pkgType))
  else
    annotObj <- pkgName

  if (!is(annotObj, pkgType))
    stop(sprintf("The object loaded with name %s is not an '%s' object.",
                 ifelse(is.character(pkgName), pkgName, gettext(callobj)[2])), pkgType)

  annotObj
}

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
