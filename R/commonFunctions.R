
####
# Multi-Sample
###################

one_ind_ms <- function(vcf, genot, ind, filterTag=NA_character_) {
    candidates <- vcf2GR(vcf, genot, ind[1, 2], filterTag)
    
    return(candidates)
} 

two_ind_ms <- function(vcf, genot, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR(vcf, genot, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR(vcf, genot, ind[2, 2], filterTag)
  
  realcommonind <- sharedVariants(gr_ind1, gr_ind2)
  candidates <- gr_ind1[realcommonind] 
  
  return(candidates)
}

three_ind_ms <- function(vcf, genot, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR(vcf, genot, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR(vcf, genot, ind[2, 2], filterTag)
  gr_ind3 <- vcf2GR(vcf, genot, ind[3, 2], filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind <- sharedVariants(filter1, gr_ind3)
  candidates <- filter1[realcommonind]
  
  return(candidates)
}

four_ind_ms <- function(vcf, genot, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR(vcf, genot, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR(vcf, genot, ind[2, 2], filterTag)
  gr_ind3 <- vcf2GR(vcf, genot, ind[3, 2], filterTag)
  gr_ind4 <- vcf2GR(vcf, genot, ind[4, 2], filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind <- sharedVariants(filter2, gr_ind4)
  candidates <- filter2[realcommonind]
  
  return(candidates)
}

five_ind_ms <- function(vcf, genot, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR(vcf, genot, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR(vcf, genot, ind[2, 2], filterTag)
  gr_ind3 <- vcf2GR(vcf, genot, ind[3, 2], filterTag)
  gr_ind4 <- vcf2GR(vcf, genot, ind[4, 2], filterTag)
  gr_ind5 <- vcf2GR(vcf, genot, ind[5, 2], filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind3 <- sharedVariants(filter2, gr_ind4)
  filter3 <- filter2[realcommonind3]
  
  realcommonind <- sharedVariants(filter3, gr_ind5)
  candidates <- filter3[realcommonind]
  
  return(candidates)
}

###

# these ones accept more than one genotipic status

one_ind_ms_2opt <- function(vcf, genot1, genot2, ind, filterTag=NA_character_) {
  candidates <- vcf2GR_2options(vcf, genot1, genot2, ind[1, 2], filterTag)
  
  return(candidates)
} 

two_ind_ms_2opt <- function(vcf, genot1, genot2, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR_2options(vcf, genot1, genot2, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR_2options(vcf, genot1, genot2, ind[2, 2], filterTag)
  
  realcommonind <- sharedVariants(gr_ind1, gr_ind2)
  candidates <- gr_ind1[realcommonind] 
  
  return(candidates)
}

three_ind_ms_2opt <- function(vcf, genot1, genot2, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR_2options(vcf, genot1, genot2, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR_2options(vcf, genot1, genot2, ind[2, 2], filterTag)
  gr_ind3 <- vcf2GR_2options(vcf, genot1, genot2, ind[3, 2], filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind <- sharedVariants(filter1, gr_ind3)
  candidates <- filter1[realcommonind]
  
  return(candidates)
}

four_ind_ms_2opt <- function(vcf, genot1, genot2, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR_2options(vcf, genot1, genot2, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR_2options(vcf, genot1, genot2, ind[2, 2], filterTag)
  gr_ind3 <- vcf2GR_2options(vcf, genot1, genot2, ind[3, 2], filterTag)
  gr_ind4 <- vcf2GR_2options(vcf, genot1, genot2, ind[4, 2], filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind <- sharedVariants(filter2, gr_ind4)
  candidates <- filter2[realcommonind]
  
  return(candidates)
}

five_ind_ms_2opt <- function(vcf, genot1, genot2, ind, filterTag=NA_character_) {
  gr_ind1 <- vcf2GR_2options(vcf, genot1, genot2, ind[1, 2], filterTag)
  gr_ind2 <- vcf2GR_2options(vcf, genot1, genot2, ind[2, 2], filterTag)
  gr_ind3 <- vcf2GR_2options(vcf, genot1, genot2, ind[3, 2], filterTag)
  gr_ind4 <- vcf2GR_2options(vcf, genot1, genot2, ind[4, 2], filterTag)
  gr_ind5 <- vcf2GR_2options(vcf, genot1, genot2, ind[5, 2], filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind3 <- sharedVariants(filter2, gr_ind4)
  filter3 <- filter2[realcommonind3]
  
  realcommonind <- sharedVariants(filter3, gr_ind5)
  candidates <- filter3[realcommonind]
  
  return(candidates)
}

#####################

####
# Unique sample
#######################

one_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd), genomeInfo)
  
  candidates <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  
  return(candidates)
} 

two_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  
  gr_ind1 <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  gr_ind2 <- vcf2GR(vcf_ind2, genot, 1, filterTag)
  
  realcommonind <- sharedVariants(gr_ind1, gr_ind2)
  candidates <- gr_ind1[realcommonind] 
  
  return(candidates)
} 

three_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeInfo)
  
  gr_ind1 <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  gr_ind2 <- vcf2GR(vcf_ind2, genot, 1, filterTag)
  gr_ind3 <- vcf2GR(vcf_ind3, genot, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind <- sharedVariants(filter1, gr_ind3)
  candidates <- filter1[realcommonind]
  
  return(candidates)
} 

four_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeInfo)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeInfo)
  
  gr_ind1 <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  gr_ind2 <- vcf2GR(vcf_ind2, genot, 1, filterTag)
  gr_ind3 <- vcf2GR(vcf_ind3, genot, 1, filterTag)
  gr_ind4 <- vcf2GR(vcf_ind4, genot, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind <- sharedVariants(filter2, gr_ind4)
  candidates <- filter2[realcommonind]
  
  return(candidates)
} 

five_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeInfo)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeInfo)
  vcf_ind5 <- readVcf(unlist(pwd[5]), genomeInfo)
  
  gr_ind1 <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  gr_ind2 <- vcf2GR(vcf_ind2, genot, 1, filterTag)
  gr_ind3 <- vcf2GR(vcf_ind3, genot, 1, filterTag)
  gr_ind4 <- vcf2GR(vcf_ind4, genot, 1, filterTag)
  gr_ind5 <- vcf2GR(vcf_ind5, genot, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind3 <- sharedVariants(filter2, gr_ind4)
  filter3 <- filter2[realcommonind3]
  
  realcommonind <- sharedVariants(filter3, gr_ind5)
  candidates <- filter3[realcommonind]
  
  return(candidates)
} 

# these ones accept more than one genotipic status

one_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd), genomeInfo)
  
  candidates <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  
  return(candidates)
} 

two_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  
  gr_ind1 <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  gr_ind2 <- vcf2GR_2options(vcf_ind2, genot1, genot2, 1, filterTag)
  
  realcommonind <- sharedVariants(gr_ind1, gr_ind2)
  candidates <- gr_ind1[realcommonind] 
  
  return(candidates)
}

three_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeInfo)
  
  gr_ind1 <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  gr_ind2 <- vcf2GR_2options(vcf_ind2, genot1, genot2, 1, filterTag)
  gr_ind3 <- vcf2GR_2options(vcf_ind3, genot1, genot2, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind <- sharedVariants(filter1, gr_ind3)
  candidates <- filter1[realcommonind]
  
  return(candidates)
}

four_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeInfo)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeInfo)
  
  gr_ind1 <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  gr_ind2 <- vcf2GR_2options(vcf_ind2, genot1, genot2, 1, filterTag)
  gr_ind3 <- vcf2GR_2options(vcf_ind3, genot1, genot2, 1, filterTag)
  gr_ind4 <- vcf2GR_2options(vcf_ind4, genot1, genot2, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind <- sharedVariants(filter2, gr_ind4)
  candidates <- filter2[realcommonind]
  
  return(candidates)
}

five_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeInfo) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeInfo)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeInfo)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeInfo)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeInfo)
  vcf_ind5 <- readVcf(unlist(pwd[5]), genomeInfo)
  
  gr_ind1 <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  gr_ind2 <- vcf2GR_2options(vcf_ind2, genot1, genot2, 1, filterTag)
  gr_ind3 <- vcf2GR_2options(vcf_ind3, genot1, genot2, 1, filterTag)
  gr_ind4 <- vcf2GR_2options(vcf_ind4, genot1, genot2, 1, filterTag)
  gr_ind5 <- vcf2GR_2options(vcf_ind5, genot1, genot2, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind2 <- sharedVariants(filter1, gr_ind3)
  filter2 <- filter1[realcommonind2]
  
  realcommonind3 <- sharedVariants(filter2, gr_ind4)
  filter3 <- filter2[realcommonind3]
  
  realcommonind <- sharedVariants(filter3, gr_ind5)
  candidates <- filter3[realcommonind]
  
  return(candidates)
}

##############

# this one creates a loop that discards shared variants with controls one by one
pullout_shared_us <- function(gr_affected, control_list, filterTag=NA_character_, genomeInfo) {
  for (i in unlist(control_list)) {
    ind <- readVcf(i, genomeInfo)
    gr_ind <- vcf2GR_2options(ind, "0/1", "1/1", 1, filterTag)
    
    realcommonind <- sharedVariants(gr_affected, gr_ind)
    gr_affected <- gr_affected[-realcommonind]
  }
  return(gr_affected)
}


############

## match chromosomes between an input GRanges object and a TxDb object
## using functionality from the GenomeInfoDb package. Further, only
## chromosomes in the input GRanges object that also occur in the input
## TxDb object will be retained. Among those, their lengths will be compared
## and those with different lengths will be discarded. This may happen between
## identical versions of the genome obtained from different sites, as with
## the mitochondrial chromosome from UCSC (older version) and b37 (newer version)
## (see https://wiki.dnanexus.com/Scientific-Notes/human-genome), or between
## different versions of the genome. When no chromosome match its length between
## the input GRanges object and the input TxDb object, an error is prompt.

matchChromosomes <- function(variantsGR, txdb) {
  message("Discarding scaffold sequences.")
  variantsGR <- keepStandardChromosomes(variantsGR)
  message(sprintf("Switching to %s chromosome-names style.", seqlevelsStyle(txdb)))
  seqlevelsStyle(variantsGR) <- seqlevelsStyle(txdb)
  genome(variantsGR) <- genome(txdb)
  ## genome(variantsGR) <- genome(txdb)[intersect(names(genome(variantsGR)), names(genome(txdb)))]
  slenVcf <- seqlengths(variantsGR)
  slenTxDb <- seqlengths(txdb)
  commonChr <- intersect(names(slenVcf), names(slenTxDb))
  slenVcf <- slenVcf[commonChr]
  slenTxDb <- slenTxDb[commonChr]
  if (any(slenVcf != slenTxDb)) {
    if (sum(slenVcf != slenTxDb) == 1) {
      warning(sprintf("Chromosome %s has different lengths between input VCF and transcript-centric annotations. This chromosome will be discarded from further analysis", paste(commonChr[which(slenVcf != slenTxDb)], collapse=", ")))
    } else {
      warning(sprintf("Chromosomes %s have different lengths between input VCF and transcript-centric annotations. These chromosomes will be discarded from further analysis", paste(commonChr[which(slenVcf != slenTxDb)], collapse=", ")))
    }
    if (sum(slenVcf == slenTxDb) == 0)
      stop("None of the chromosomes in the input VCF file has the same length as the chromosomes in the transcript-centric annotations. Versions of the genome reference sequence may be different between the input VCF file and the transcript-centric annotation package.")
    variantsGR <- keepSeqlevels(variantsGR, commonChr[slenVcf == slenTxDb])
  }

  variantsGR
}


sharedVariants <- function(query1, query2) {
  overgen12 <- findOverlaps(query1, query2, type="equal", ignore.strand=FALSE)
  commongen12 <- queryHits(overgen12)
  maskcarr12 <- unlist(query1[queryHits(overgen12)]$ALT) == unlist(query2[subjectHits(overgen12)]$ALT)
  realcommongen12 <- commongen12[maskcarr12]
  realcommongen12
}

vcf2GR <- function(vcf, genot, ind, filterTag){
  maskgenvcf <- geno(vcf)$GT[, ind] == genot
  genvcf <- vcf[maskgenvcf, ]
  gr <- rowData(genvcf)
  if (sum(elementLengths(gr$ALT) > 1) >= 1) {
    biallelic <- which(elementLengths(gr$ALT) > 1)
    gr <- gr[-biallelic]
  }
  if (!is.na(filterTag))
    gr <- gr[gr$FILTER %in% filterTag, ]  
  gr
}


vcf2GR_2options <- function(vcf, genot1, genot2, ind, filterTag){
  maskgenvcf <- geno(vcf)$GT[, ind] == genot1 | geno(vcf)$GT[, ind] == genot2
  genvcf <- vcf[maskgenvcf, ]
  gr <- rowData(genvcf)
  if (sum(elementLengths(gr$ALT) > 1) >= 1) {
    biallelic <- which(elementLengths(gr$ALT) > 1)
    gr <- gr[-biallelic]
  }
  if (!is.na(filterTag))
    gr <- gr[gr$FILTER %in% filterTag, ]
  gr
}

## this should be removed
retrieveLostInfo <- function(GRbeforeLocateVariants, GRafterLocateVariants) {
  overgen12 <- findOverlaps(GRbeforeLocateVariants, GRafterLocateVariants, type="equal", ignore.strand=FALSE)
  commongen12 <- queryHits(overgen12)
  GRbeforeLocateVariants[commongen12]
}

## adjustForStrandSense returns the 5'->3' strand of sequences
## stored in DNAStringSet or DNAStringSetList objects 
## the returned object is of the same class as the input object
adjustForStrandSense <- function(variantsGR, alleles){
  adjustedAlleles <- alleles
  nstrand <- strand(variantsGR) == "-"
  if (any(nstrand)) {
    if (class(alleles) == "DNAStringSet")
      adjustedAlleles[nstrand] <- reverseComplement(alleles[nstrand])
    else {
      if (class(alleles) == "DNAStringSetList")
        adjustedAlleles[nstrand] <- relist(reverseComplement(unlist(alleles[nstrand])),
                                           alleles[nstrand])
      else
        stop("adjustForStrandSense: argument 'alleles' should be either a 'DNAStringSet' or a 'DNAStringSetList' object.")
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

  exonsbytx <- exonsBy(txdb, by="tx")
  map <- mapCoords(gr, exonsbytx, ignore.strand=FALSE, elt.hits=TRUE)
  qolap <- map$queryHits
  res <- gr[qolap]
  mcols(res) <- append(values(res),
                       DataFrame(cDNALOC=ranges(map), TXID2=map$subjectHits))
  ## restrict results to cDNA positions of query TXID
  ## this is necessary with overlaping transcripts
  res <- res[res$TXID == res$TXID2]

  start(res$cDNALOC)
}
