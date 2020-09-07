##
## helper functions not exported
##

## harmonize sequence style between variants, annotations and genome sequence
## because variants, annotations and genome sequence may originate from the
## same genome under different nomenclatures (e.g., GRCh37.p13, hg19, hs37d5, etc.)
## we'll assume that genome sequence is the same when there is at least one
## sequence that matches in name (up to the known "styles") and length.
.matchSeqinfo <- function(variantsGR, txdb, bsgenome) {

  stopifnot(class(variantsGR) == "VRanges") ## QC
  stopifnot(class(txdb) == "TxDb") ## QC
  stopifnot(class(bsgenome) == "BSgenome") ## QC

  ## set temporarily the sequence style of variants to that of
  ## the genome and discard chromosomes that do not have the same
  ## length. this is basically to handle situations such as with
  ## NCBI GRCh37 / h19 vs b37 where the MT chromosome has different lengths
  seqlevelsStyle(variantsGR) <- seqlevelsStyle(bsgenome)[1]
  slenVcf <- seqlengths(variantsGR)
  slenBSgenome <- seqlengths(bsgenome)
  commonChr <- intersect(names(slenVcf), names(slenBSgenome))
  slenVcf <- slenVcf[commonChr]
  slenBSgenome <- slenBSgenome[commonChr]
  if (any(slenVcf != slenBSgenome)) {
    if (sum(slenVcf != slenBSgenome) == 1) {
      message(sprintf("Chromosome %s has different lengths between the input VCF and the BSgenome pakage. This chromosome will be discarded from further analysis", paste(commonChr[which(slenVcf != slenBSgenome)], collapse=", ")))
    } else {
      message(sprintf("Chromosomes %s have different lengths between the input VCF and the BSgenome package. These chromosomes will be discarded from further analysis", paste(commonChr[which(slenVcf != slenBSgenome)], collapse=", ")))
    }
    if (sum(slenVcf == slenBSgenome) == 0)
      stop("None of the chromosomes in the input VCF file has the same length as the chromosomes in the input BSgenome package. The genome reference sequence employed to generate the VCF file was probably different from the one in the input BSgenome package.")
    variantsGR <- keepSeqlevels(variantsGR, commonChr[slenVcf == slenBSgenome], pruning.mode="coarse")
    commonChr <- commonChr[slenVcf == slenBSgenome]
  }

  ## set the genome information to the one of the genome package, which should be complete
  message(sprintf("Assuming the genome build of the input variants is %s.",
                  unique(genome(bsgenome)[commonChr])))
  seqinfo(variantsGR, new2old=match(commonChr, seqlevels(variantsGR))) <- seqinfo(bsgenome)[commonChr]

  ## set the sequence style of variants to the one of annotations
  message(sprintf("Switching to the %s chromosome-name style from the transcript-centric annotation package.",
                  seqlevelsStyle(txdb)))
  seqlevelsStyle(variantsGR) <- seqlevelsStyle(txdb)[1]
  commonChr <- intersect(seqlevels(variantsGR), seqlevels(txdb))

  ## subset further the sequences to analyze to those whose length also matches the ones of annotations
  slenVcf <- seqlengths(variantsGR)[commonChr]
  slenTxDb <- seqlengths(txdb)[commonChr]
  if (any(slenVcf != slenTxDb)) {
    if (sum(slenVcf != slenTxDb) == 1) {
      message(sprintf("Chromosome %s has different lengths between the input VCF and the input TxDb pakage. This chromosome will be discarded from further analysis", paste(commonChr[which(slenVcf != slenTxDb)], collapse=", ")))
    } else {
      message(sprintf("Chromosomes %s have different lengths between the input VCF and the input TxDb package. These chromosomes will be discarded from further analysis", paste(commonChr[which(slenVcf != slenTxDb)], collapse=", ")))
    }
    if (sum(slenVcf == slenTxDb) == 0)
      stop("None of the chromosomes in the input VCF file has the same length as the chromosomes in the input TxDb package. The genome reference sequence employed to generate the VCF file was probably different from the one in the input TxDb package.")
    variantsGR <- keepSeqlevels(variantsGR, commonChr[slenVcf == slenTxDb], pruning.mode="coarse")
    commonChr <- commonChr[slenVcf == slenTxDb]
  }

  ## inform the user the genome build is going to be the one of the
  ## variants and genome package, in case it does not match the one of the annotations
  ## the genome build is going to be the one of the annotations assuming it's the same
  if (any(is.na(genome(txdb)[commonChr])))
    message(sprintf("Assuming the genome build of transcript-centric annotations is %s.",
                    unique(genome(variantsGR)[commonChr])))
  else if (any(genome(variantsGR)[commonChr] != genome(txdb)[commonChr])) {
    message(sprintf("Assumming %s represent the same genome build.",
                    paste(c(unique(genome(variantsGR)[commonChr]), unique(genome(txdb)[commonChr])),
                          collapse=" and ")))
    seqinfo(variantsGR, new2old=match(commonChr, seqlevels(variantsGR))) <- seqinfo(txdb)[commonChr]
  }

  ## discard scaffold sequences
  message("Discarding scaffold sequences.")
  variantsGR <- keepStandardChromosomes(variantsGR)

  variantsGR
}

## return the 5'->3' strand of 'alleles' according to
## the 'LOCSTRAND' column in the input 'VRanges' object 'variantsVR'
## the returned object is of the same class as the input 'alleles' object
.adjustForStrandSense <- function(variantsVR, alleles){
  stopifnot(is(variantsVR, "VRanges")) ## QC
  stopifnot(!is.null(variantsVR$LOCSTRAND)) ## QC

  adjustedAlleles <- alleles
  nstrand <- variantsVR$LOCSTRAND == "-"
  if (any(nstrand)) {
    if (class(alleles) == "DNAStringSet")
      adjustedAlleles[nstrand] <- reverseComplement(alleles[nstrand])
    else if (class(alleles) == "DNAStringSetList") {
      adjustedAlleles[nstrand] <- relist(reverseComplement(unlist(alleles[nstrand])), alleles[nstrand])
    } else if (class(alleles) == "character" || class(alleles) == "characterRle") {
      nstrand <- as.vector(nstrand)
      adjustedAlleles[nstrand] <- as.character(reverseComplement(DNAStringSet(alleles[nstrand])))
    } else {
        stop(".adjustForStrandSense: argument 'alleles' should be either a 'DNAStringSet' object, a 'DNAStringSetList' object, a 'characterRle' object or a 'character' vector.")
    }
  }

  adjustedAlleles
}

## get cDNA position
## it expects that gr has a transcript identifier column called TXID
## resulting from a previous call to locateVariants()
.cDNAloc <- function(gr, txdb) {
  if (!"TXID" %in% colnames(mcols(gr)))
    stop("cDNApos: metadata column 'TXID' not found in input GRanges object.")

  exonsbytx <- exonsBy(txdb, by="tx") ## CONSIDER USING THE CACHE HERE !!!
  map <- mapToTranscripts(gr, exonsbytx, ignore.strand=FALSE)
  qolap <- map$xHits
  res <- gr[qolap]
  mcols(res) <- append(values(res),
                       DataFrame(cDNALOC=ranges(map), TXID2=map$transcriptsHits))
  ## restrict results to cDNA positions of query TXID
  ## this is necessary with overlaping transcripts
  res <- res[res$TXID == res$TXID2]

  start(res$cDNALOC)

  ## this could probably be done faster by doing something like
  ## cDNA <- mapToTranscripts(gr[!duplicated(ranges(gr))], exonsbytx)
  ## cDNA <- ranges(cDNA)[togroup(gr$QUERYID)]
}

.readPEDfile <- function(pedFilename) {

  if (!file.exists(pedFilename))
    stop(sprintf("could not open the PED file %s.", pedFilename))

  pedDf <- read.table(pedFilename, header=FALSE, stringsAsFactors=FALSE)
  
  if (dim(pedDf)[2]==6) {
     
      pedDf <- pedDf[, 1:6]
      colnames(pedDf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Sex", "Phenotype") 
      
  } else if (dim(pedDf)[2]==7) {
      
      pedDf <- pedDf[, 1:7]
      colnames(pedDf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Sex", "Phenotype", "Age") 
      message("Age is present as the 7th column in the PED file")
      
  } else { stop("Not standard PED file format. Required columns: ",
                "'FamilyID', 'IndividualID', 'FatherID', 'MotherID', 'Sex', 'Phenotype'",
                "Optionally another 'Age' column can be attached.") }
  
  if (is.character(pedDf$Sex)) {
    if (all(sort(unique(pedDf$Sex)) == c("f", "m")))
      pedDf$Sex <- match(pedDf$Sex, c("m", "f"))
    else
      stop("Character values in the 'Sex' column of the PED file are not valid. ",
           "Please use only 'm' for males and 'f' for females.")
  }

  if (is.character(pedDf$Phenotype)) {
    if (all(unique(pedDf$Phenotype) %in% c("a", "h")))
      pedDf$Phenotype <- match(pedDf$Phenotype, c("h", "a"))
    else if (all(unique(pedDf$Phenotype) %in% c("a", "h", "u")))
      pedDf$Phenotype <- match(pedDf$Phenotype, c("u", "h", "a")) - 1L
    else
      stop("Character values in the 'Phenotype' column of the PED file are not valid. ",
           "Please use only 'a' for affected, 'h' for healthy and 'u' for unknown.")
  }

  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedDf$Phenotype == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")

  maskFatherIDs <- pedDf$FatherID != 0
  maskMotherIDs <- pedDf$MotherID != 0
  if (any(maskFatherIDs != maskMotherIDs))
    stop("Father and mother identifiers must be both either zero or non-zero.")

  if (any(duplicated(pedDf$IndividualID))) {
    pedDf$IndividualID <- paste(pedDf$FamilyID, pedDf$IndividualID, sep="_")
    pedDf$FatherID[maskFatherIDs] <- paste(pedDf$FamilyID[maskFatherIDs],
                                           pedDf$FatherID[maskFatherIDs], sep="_")
    pedDf$MotherID[maskMotherIDs] <- paste(pedDf$FamilyID[maskMotherIDs],
                                           pedDf$MotherID[maskMotherIDs], sep="_")
    if (!duplicated(pedDf$IndividualID))
      message("Individual identifiers have been made unique.")
    else
      stop("Individual identifiers are not unique, even when combined with ",
           "family identifiers.")
  }

  pedDf
}

## get the operating system
## adapted from http://www.r-bloggers.com/identifying-the-os-from-r/
.getOS <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else {
    ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

## load precomputed results for the vignette
## this is useful to meet building times in
## the windows build machines

.loadPrecomputedVariantFilteringResults <- function() {
  uind <- readRDS(system.file("extdata", "uind.rds", package="VariantFiltering"))
  uind
}
