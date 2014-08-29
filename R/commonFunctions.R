
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

one_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd), genomeVersion)
  
  candidates <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  
  return(candidates)
} 

two_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  
  gr_ind1 <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  gr_ind2 <- vcf2GR(vcf_ind2, genot, 1, filterTag)
  
  realcommonind <- sharedVariants(gr_ind1, gr_ind2)
  candidates <- gr_ind1[realcommonind] 
  
  return(candidates)
} 

three_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeVersion)
  
  gr_ind1 <- vcf2GR(vcf_ind1, genot, 1, filterTag)
  gr_ind2 <- vcf2GR(vcf_ind2, genot, 1, filterTag)
  gr_ind3 <- vcf2GR(vcf_ind3, genot, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind <- sharedVariants(filter1, gr_ind3)
  candidates <- filter1[realcommonind]
  
  return(candidates)
} 

four_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeVersion)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeVersion)
  
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

five_ind_us <- function(pwd, genot, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeVersion)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeVersion)
  vcf_ind5 <- readVcf(unlist(pwd[5]), genomeVersion)
  
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

one_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd), genomeVersion)
  
  candidates <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  
  return(candidates)
} 

two_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  
  gr_ind1 <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  gr_ind2 <- vcf2GR_2options(vcf_ind2, genot1, genot2, 1, filterTag)
  
  realcommonind <- sharedVariants(gr_ind1, gr_ind2)
  candidates <- gr_ind1[realcommonind] 
  
  return(candidates)
}

three_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeVersion)
  
  gr_ind1 <- vcf2GR_2options(vcf_ind1, genot1, genot2, 1, filterTag)
  gr_ind2 <- vcf2GR_2options(vcf_ind2, genot1, genot2, 1, filterTag)
  gr_ind3 <- vcf2GR_2options(vcf_ind3, genot1, genot2, 1, filterTag)
  
  realcommonind1 <- sharedVariants(gr_ind1, gr_ind2)
  filter1 <- gr_ind1[realcommonind1]
  
  realcommonind <- sharedVariants(filter1, gr_ind3)
  candidates <- filter1[realcommonind]
  
  return(candidates)
}

four_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeVersion)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeVersion)
  
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

five_ind_us_2opt <- function(pwd, genot1, genot2, filterTag=NA_character_, genomeVersion) {
  vcf_ind1 <- readVcf(unlist(pwd[1]), genomeVersion)
  vcf_ind2 <- readVcf(unlist(pwd[2]), genomeVersion)
  vcf_ind3 <- readVcf(unlist(pwd[3]), genomeVersion)
  vcf_ind4 <- readVcf(unlist(pwd[4]), genomeVersion)
  vcf_ind5 <- readVcf(unlist(pwd[5]), genomeVersion)
  
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
pullout_shared_us <- function(gr_affected, control_list, filterTag=NA_character_, genomeVersion) {
  for (i in unlist(control_list)) {
    ind <- readVcf(i, genomeVersion)
    gr_ind <- vcf2GR_2options(ind, "0/1", "1/1", 1, filterTag)
    
    realcommonind <- sharedVariants(gr_affected, gr_ind)
    gr_affected <- gr_affected[-realcommonind]
  }
  return(gr_affected)
}


############

## match chromosome names between the queryGRanges and a TxDb object
## by renaming to the UCSC nomenclature any of the two input objects. This is
## done by simply pasting the 'chr' prefix. However, this does not handle other
## mismatches between unlocalized and unplaced sequences, alternate loci sequences
## and mitochondrial sequences (M vs MT). The latter, though, seem to be different
## between UCSC and b37 (see https://wiki.dnanexus.com/Scientific-Notes/human-genome),
## and therefore, it may be sensible to leave it unmatched.

## note that although the levels of the GRanges object are by returning the
## updated object the levels of the input 'txdb' are renamed through the argument
## itself (i.e., by reference)
matchChromosomeNames <- function(variantsGRanges, txdb) {

  if (class(variantsGRanges) != "GRanges")
    stop("matchChromosomeNames: argument 'variantsGRanges' is not a GRanges object.")

  if (class(txdb) != "TxDb")
    stop("matchChromosomeNames: argument 'txdb' is not a TxDb object.")

  if (identical(seqlevels(variantsGRanges), seqlevels(txdb)))
    return(variantsGRanges)

  vcfUCSC <- TRUE
  if (substr(seqlevels(variantsGRanges), 1, 3)[1] != "chr")
    vcfUCSC <- FALSE
  
  txdbUCSC <- TRUE
  if (unique(substr(seqlevels(txdb), 1, 3)) != "chr")
    txdbUCSC <- FALSE
  
  if (vcfUCSC == TRUE && txdbUCSC == FALSE) {
    message("Renaming UCSC chromosome names from variants into b37 (1000 Genomes Project) nomenclature (mitochondrial, unlocalized, unplaced and alternate loci sequences may not be properly renamed)")
    txdb <- renameSeqlevels(txdb, paste0("chr", seqlevels(txdb)))
  } else if (vcfUCSC == FALSE && txdbUCSC == TRUE) {
    message("Renaming b37 (1000 Genomes Project) chromosome names from variants into UCSC nomenclature (mitochondrial, unlocalized, unplaced and alternate loci sequences may not be properly renamed)")
    variantsGRanges <- renameSeqlevels(variantsGRanges, paste0("chr", seqlevels(variantsGRanges)))
  }

  if (!identical(seqlevels(variantsGRanges), seqlevels(txdb))) {
    message("Discarding variants whose chromosome names do not match the UCSC nomenclature")
    seqlevels(variantsGRanges) <- intersect(seqlevels(variantsGRanges), seqlevels(txdb))
  }

  ## complete the seqinfo data
  isCircular(variantsGRanges) <- isCircular(txdb)[names(isCircular(variantsGRanges))]

  variantsGRanges
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


vcf2GR_2options<- function(vcf, genot1, genot2, ind, filterTag){
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
  map <- mapCoords(gr, exonsbytx, elt.Hits=TRUE)
  eolap <- map$eltHits
  qolap <- map$queryHits
  txids <- rep(names(exonsbytx), elementLengths(exonsbytx))[eolap]
  res <- gr[qolap]
  mcols(res) <- append(values(res),
                       DataFrame(cDNALOC=ranges(map), TXID2=as.integer(txids)))
  ## restrict results to cDNA positions of query TXID
  res <- res[res$TXID == res$TXID2]

  start(res$cDNALOC)
}
