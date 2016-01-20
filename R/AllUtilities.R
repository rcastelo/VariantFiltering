##
## helper functions not exported
##

## harmonize sequence style between variants, annotations and genome sequence
.matchSeqinfo <- function(variantsGR, txdb, bsgenome) {

  ## set temporarily the sequence style of variants to that of
  ## the genome and discard chromosomes that do not have the same
  ## length. this is basically to handle situations such as with
  ## NCBI GRCh37 / h19 vs b37 where the MT chromosome has different lengths
  seqlevelsStyle(variantsGR) <- seqlevelsStyle(bsgenome)
  slenVcf <- seqlengths(variantsGR)
  slenBSgenome <- seqlengths(bsgenome)
  commonChr <- intersect(names(slenVcf), names(slenBSgenome))
  slenVcf <- slenVcf[commonChr]
  slenBSgenome <- slenBSgenome[commonChr]
  if (any(slenVcf != slenBSgenome)) {
    if (sum(slenVcf != slenBSgenome) == 1) {
      warning(sprintf("Chromosome %s has different lengths between the input VCF and the input BSgenome pakage. This chromosome will be discarded from further analysis", paste(commonChr[which(slenVcf != slenBSgenome)], collapse=", ")))
    } else {
      warning(sprintf("Chromosomes %s have different lengths between the input VCF and the input BSgenome package. These chromosomes will be discarded from further analysis", paste(commonChr[which(slenVcf != slenBSgenome)], collapse=", ")))
    }
    if (sum(slenVcf == slenBSgenome) == 0)
      stop("None of the chromosomes in the input VCF file has the same length as the chromosomes in the input BSgenome package. The genome reference sequence employed to generate the VCF file was probably different from the one in the input BSgenome package.")
    variantsGR <- keepSeqlevels(variantsGR, commonChr[slenVcf == slenBSgenome])
    commonChr <- commonChr[slenVcf == slenBSgenome]
  }

  ## set the genome build to the one of the genome package
  message(sprintf("Assuming the genome build of the input variants is %s.",
                  unique(genome(bsgenome)[commonChr])))
  seqinfo(variantsGR, new2old=match(commonChr, seqlevels(variantsGR))) <- seqinfo(bsgenome)[commonChr]

  ## set the sequence style of variants to the one of annotations
  message(sprintf("Switching to the %s chromosome-name style from the transcript-centric annotation package.",
                  seqlevelsStyle(txdb)))
  seqlevelsStyle(variantsGR) <- seqlevelsStyle(txdb)
  commonChr <- intersect(seqlevels(variantsGR), seqlevels(txdb))

  ## inform the user the genome build is going to be the one of the
  ## variants and genome package, in case it does not match the one of the annotations
  if (any(is.na(genome(txdb)[commonChr])))
    warning(sprintf("Assuming the genome build of transcript-centric annotations is %s.",
                    unique(genome(variantsGR)[commonChr])))
  else if (any(genome(variantsGR)[commonChr] != genome(txdb)[commonChr])) {
    warning(sprintf("Assumming %s represent the same genome build.",
                    paste(c(unique(genome(variantsGR)[commonChr]), unique(genome(txdb)[commonChr])),
                          collapse=" and ")))
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
  pedDf <- pedDf[, 1:6]
  colnames(pedDf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Gender", "Phenotype")

  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedDf$Phenotype  == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")

  pedDf
}
