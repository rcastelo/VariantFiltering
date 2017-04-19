VariantFilteringParam <- function(vcfFilenames, pedFilename=NA_character_,
                                  bsgenome="BSgenome.Hsapiens.1000genomes.hs37d5",
                                  orgdb="org.Hs.eg.db",
                                  txdb="TxDb.Hsapiens.UCSC.hg19.knownGene",
                                  snpdb="SNPlocs.Hsapiens.dbSNP144.GRCh37",
                                  weightMatricesFilenames=NA,
                                  weightMatricesLocations=rep(list(variantLocations()), length(weightMatricesFilenames)),
                                  weightMatricesStrictLocations=rep(list(FALSE), length(weightMatricesFilenames)),
                                  radicalAAchangeFilename=file.path(system.file("extdata", package="VariantFiltering"),
                                                                    "AA_chemical_properties_HanadaGojoboriLi2006.tsv"),
                                  codonusageFilename=file.path(system.file("extdata", package="VariantFiltering"),
                                                               "humanCodonUsage.txt"),
                                  geneticCode=getGeneticCode("SGC0"),
                                  allTranscripts=FALSE,
                                  regionAnnotations=list(CodingVariants(), IntronVariants(),
                                                         FiveSpliceSiteVariants(), ThreeSpliceSiteVariants(),
                                                         PromoterVariants(), FiveUTRVariants(), ThreeUTRVariants()),
                                  otherAnnotations=c("MafDb.1Kgenomes.phase1.hs37d5",
                                                     "PolyPhen.Hsapiens.dbSNP131",
                                                     "SIFT.Hsapiens.dbSNP137",
                                                     "phastCons100way.UCSC.hg19",
                                                     "humanGenesPhylostrata"),
                                  geneKeytype=NA_character_,
                                  yieldSize=NA_integer_) {

  ## store call for reproducing it later
  callobj <- match.call()
  callstr <- gsub(".local", "VariantFilteringParam", deparse(callobj))
  callstr <- gsub("= vcfFilenames", sprintf("= c(%s)", paste(sprintf("\"%s\"", vcfFilenames), collapse=", ")), callstr)

  if (length(vcfFilenames) > 1) {
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with vcftools, do the variant calling simultaneously on all samples, or proceed analyzing each file separately.")
  } else if (length(vcfFilenames) < 1)
    stop("Missing VCF filename.")

  ## check if input VCF, PED, splice site and radical AA change files exist
  tryCatch({
    .io_check_exists(c(vcfFilenames, radicalAAchangeFilename, codonusageFilename))
  }, error=function(err) {
       stop(conditionMessage(err), call.=FALSE)
  })

  if (!is.na(pedFilename) && length(pedFilename) == 1)
    tryCatch({
      .io_check_exists(pedFilename)
    }, error=function(err) {
         stop(conditionMessage(err), call.=FALSE)
    })
  else
    pedFilename <- NA_character_

  if (any(!is.na(weightMatricesFilenames)))
    tryCatch({
      .io_check_exists(weightMatricesFilenames)
    }, error=function(err) {
         stop(conditionMessage(err), call.=FALSE)
    })
  else
    weightMatricesFilenames <- NA_character_

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
  for (i in seq_along(vcfFilenames)) {
    hd <- scanVcfHeader(vcfFilenames[i])
    seqinfos[[i]] <- seqinfo(hd)
    sampleNames <- c(sampleNames, samples(hd))
    tfl[[i]] <- TabixFile(vcfFilenames[i], yieldSize=yieldSize)
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

  ## read weight matrices
  weightMatrices <- list()
  if (any(!is.na(weightMatricesFilenames))) {
    if (class(weightMatricesFilenames) != "character")
      stop("'weightMatricesFilenames' should be a character vector.")
    if (any(!is.na(weightMatricesLocations)) && class(unlist(weightMatricesLocations, use.names=FALSE)) != "character")
      stop("elements of 'weightMatricesLocations' should be a character vector.")
    if (any(!is.na(weightMatricesLocations)) && length(weightMatricesFilenames) > 1 && class(weightMatricesLocations) != "list")
      stop("'weightMatricesLocations' should be a list of character vectors when more than one weight matrix filename is given to 'weightMatricesFilenames'")
    if (any(!is.na(weightMatricesLocations)) && length(weightMatricesFilenames) == 1 && class(weightMatricesLocations) != "list")
      weightMatricesLocations <- list(weightMatricesLocations)
    if (any(!is.na(weightMatricesStrictLocations)) && length(weightMatricesFilenames) > 1 && class(weightMatricesStrictLocations) != "list")
      stop("'weightMatricesStrictLocations' should be a list of logical vectors when more than one weight matrix filename is given to 'weightMatricesFilenames'")
    if (any(!is.na(weightMatricesStrictLocations)) && length(weightMatricesFilenames) == 1 && class(weightMatricesStrictLocations) != "list")
      weightMatricesStrictLocations <- list(weightMatricesStrictLocations)

    weightMatrices <- mapply(function(fname, loc, sloc) readWm(fname, loc, sloc),
                             as.list(weightMatricesFilenames), weightMatricesLocations, weightMatricesStrictLocations)
  }

  ## read radical amino acid change matrix
  radicalAAchangeMatrix <- readAAradicalChangeMatrix(radicalAAchangeFilename)

  ## read codon usage table
  codonusageTable <- read.table(file=codonusageFilename, sep=";")
  codonusageTable <- do.call("names<-", list(codonusageTable[[3]], codonusageTable[[1]]))
  
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

  snpdblst <- as.list(snpdb)
  names(snpdblst) <- snpdb
  ## assume the first bit of the name of the package contains the class name
  snpdb <- lapply(as.list(snpdblst),
                  function(pkg) .loadAnnotationPackageObject(pkg, "snpdb",
                                                             strsplit(pkg, ".", fixed=TRUE)[[1]][1]))

  if (!is.logical(allTranscripts))
    stop("argument 'allTranscripts' should be logical.")

  if (!is.character(otherAnnotations))
    stop("argument 'otherAnnotations' should contain names of objects or installed packages that can be used by VariantFiltering to annotate genetic variants.")

  ## check that object determining region annotations belong to the 'VariantType' class
  if (!is.list(regionAnnotations))
    regionAnnotations <- list(regionAnnotations)

  mask <- sapply(regionAnnotations, function(x) !is(x, "VariantType"))
  if (any(mask))
    stop("The argument 'regionAnnotations' should contain a list of one or more 'VariantType' objects. See the help page of CodingVariants() from the 'VariantAnnotation' package.")
  names(regionAnnotations) <- sapply(regionAnnotations, class)
  regionAnnotations <- SimpleList(regionAnnotations)

  ## fetch other annotation packages
  otherannotations <- list()
  for (name in otherAnnotations)  {
    if (!exists(name)) {
      if (!name %in% installed.packages(noCache=TRUE)[, "Package"])
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

  ## set the empty sequence ontology graph (sequence variant) and its matrix of descendants
  gSO <- sequence_variant.gSOXP
  gSOdmat <- .buildDescendantsMatrix(gSO)
  gSOamat <- .buildAncestorsMatrix(gSO)

  ## build quality filters contained in the input VCF
  ## somehow sometimes filters described in the VCF header
  ## do not show up in the softFilterMatrix of the VRanges object
  ## to deal with this we check whether the column exists and if not
  ## just return a mask of truth logical values
  qualityFilterDescriptions <- fixed(scanVcfHeader(tfl[[1]]))$FILTER
  qualityFilterNames <- make.names(rownames(qualityFilterDescriptions), unique=TRUE)
  rownames(qualityFilterDescriptions) <- qualityFilterNames
  qfilters <- sapply(qualityFilterNames,
                     function(qfname) {
                       f <- sprintf("function(x) { sfm <- softFilterMatrix(allVariants(x, groupBy=\"nothing\")) ; mask <- rep(TRUE, nrow(sfm)) ; if (!is.na(match(qfname, colnames(sfm)))) mask <- sfm[, \"%s\"] ; mask }", qfname)
                       eval(parse(text=f))
                     })
  qualityFR <- FilterRules(qfilters)

  ## add default functional annotation filters
  defaultFR <- FilterRules(.defaultFilters)

  ## initially default functional annotation filters are not active
  active(defaultFR) <- FALSE

  defaultFR <- c(qualityFR, defaultFR)
  defaultFilterDesc <- rbind(qualityFilterDescriptions, .defaultFilterDescriptions)

  new("VariantFilteringParam", callObj=callobj, callStr=callstr, vcfFiles=tfl, seqInfos=seqinfos,
      sampleNames=sampleNames, pedFilename=pedFilename, bsgenome=bsgenome, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
      gSO=gSO, gSOdmat=gSOdmat, gSOamat=gSOamat, weightMatrices=weightMatrices,
      radicalAAchangeFilename=radicalAAchangeFilename, radicalAAchangeMatrix=radicalAAchangeMatrix,
      codonusageFilename=codonusageFilename, codonusageTable=codonusageTable, geneticCode=geneticCode,
      regionAnnotations=regionAnnotations, otherAnnotations=otherannotations, allTranscripts=allTranscripts,
      filters=defaultFR, filterDescriptions=defaultFilterDesc, qualityFilterNames=qualityFilterNames,
      cutoffs=.defaultCutoffs, geneKeytype=geneKeytype, yieldSize=as.integer(yieldSize))
}

## functions with default parameters

spliceSiteMatricesHuman <- function() {
  c(file.path(system.file("extdata", package="VariantFiltering"), "hsap.donors.hcmc10_15_1.ibn"),
    file.path(system.file("extdata", package="VariantFiltering"), "hsap.acceptors.hcmc10_15_1.ibn"))
}

## show method

setMethod("show", signature(object="VariantFilteringParam"),
          function(object) {
            cat("\nVariantFiltering parameter object\n\n")
            cat("  VCF file(s):")
            for (f in path(object$vcfFiles))
              cat(sprintf(" %s", basename(f)))
            sampleNames <- ifelse(length(object$sampleNames) <= 4, paste(object$sampleNames, collapse=", "),
                                  paste(paste(head(object$sampleNames, n=4), collapse=", "), "...", sep=", "))
            cat(sprintf("\n  Genome version(s):"))
            for (i in seq_along(object$seqInfos))
              if (length(object$seqInfos[[i]]) > 0) {
                ## sometimes the 'genome' tag comes flanked by \" :-? let's strip that out for pure aesthetics
                genomestr <- paste(gsub("\\\"", "", unique(genome(object$seqInfos[[i]]))), collapse=",")
                ## we subset to the first element of the value returned by seqlevelsStyle()
                ## to deal with cases in which only a subset of chrosomes is contained in the
                ## input VCF (typically for teaching/example/illustration purposes) which
                ## matches more than one chromosome style
                cat(sprintf(" %s(%s)", genomestr, seqlevelsStyle(object$seqInfos[[i]])[1]))
              } else
                cat(" NA")

            cat(sprintf("\n  Number of individuals: %d (%s)\n", length(object$sampleNames), sampleNames))
            if (!is.na(object$pedFilename) && length(object$pedFilename) > 0)
              cat(sprintf("  PED file: %s\n", basename(object$pedFilename)))
            cat(sprintf("  Genome-centric annotation package: %s (%s %s %s)\n",
                        object$bsgenome@pkgname, provider(object$bsgenome),
                        providerVersion(object$bsgenome), releaseName(object$bsgenome)))
            cat(sprintf("  Variant-centric annotation package: %s (%s %s)\n",
                        names(object$snpdb), sapply(object$snpdb, provider), sapply(object$snpdb, releaseName)),
                sep="")
            if (length(object$txdb$packageName) > 0)
              cat(sprintf("  Transcript-centric annotation package: %s\n", object$txdb$packageName))
            else
              cat(sprintf("  Transcript-centric annotation table: %s\n", metadata(object$txdb)[grep("Table", metadata(object$txdb)$name), "value"]))
            cat(sprintf("  Gene-centric annotation package: %s\n", object$orgdb$packageName))
            if (length(object$weightMatrices) > 0)
              cat(sprintf("  Weight matrices: %s\n", paste(sapply(object$weightMatrices, wmName), collapse=", ")))
            cat(sprintf("  Radical/Conservative AA changes: %s\n", basename(object$radicalAAchangeFilename)))
            cat(sprintf("  Codon usage table: %s\n", basename(object$codonusageFilename)))
            cat(sprintf("  Regions to annotate: %s\n", paste(names(object$regionAnnotations), collapse=", ")))
            if (length(object$otherAnnotations) > 0)
              cat(sprintf("  Other annotation pkg/obj: %s\n",
                          paste(names(object$otherAnnotations),
                                collapse=",\n                            ")))
            cat(sprintf("  All transcripts: %s\n", object$allTranscripts))
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

## getters

setMethod("filters", signature(x="VariantFilteringParam"),
          function (x) {
            x@filters
          })

setMethod("cutoffs", signature(x="VariantFilteringParam"),
          function (x) {
            x@cutoffs
          })

setMethod("sog", signature(x="VariantFilteringParam"),
          function(x) {
            x@gSO
          })

setMethod("sodmat", signature(x="VariantFilteringParam"),
          function(x) {
            x@gSOdmat
          })

setMethod("soamat", signature(x="VariantFilteringParam"),
          function(x) {
            x@gSOamat
          })

## private functions

.isPkgLoaded <- function(name) {
   (paste("package:", name, sep="") %in% search()) ## || (name %in% loadedNamespaces()) <- problematic with cleanEx()
}


.loadAnnotationPackageObject <- function(pkgName, argName, pkgType) {

  callobj <- match.call()
  annotObj <- NULL

  if (is.character(pkgName)) {
    if (!pkgName %in% installed.packages(noCache=TRUE)[, "Package"])
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
