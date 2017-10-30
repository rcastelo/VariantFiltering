
setMethod("autosomalRecessiveHeterozygous", signature(param="VariantFilteringParam"),
          function(param, svparam=ScanVcfParam(), BPPARAM=bpparam("SerialParam")) {

  ## store call for reproducing it later
  callobj <- match.call()
  callstr <- gsub(".local", "autosomalRecessiveHeterozygous", deparse(callobj))

  ## fetch necessary parameters
  vcfFiles <- param$vcfFiles
  ped <- param$pedFilename
  seqInfos <- param$seqInfos
  txdb <- get(param$txdb)
  bsgenome <- get(param$bsgenome)
  sampleNames <- param$sampleNames

  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))

  if (length(vcfFiles) > 1)
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with software such as vcftools or GATK, or do the variant calling simultaneously on all samples, or proceed analyzing each file separately.")
  else if (length(vcfFiles) < 1)
    stop("A minimum of 1 VCF file has to be provided.")
  
  if (is.na(ped))
    stop("Please specify a PED file name when building the parameter object.")

  if (!file.exists(ped))
    stop(sprintf("could not open the PED file %s.", ped))

  ## REVISE THIS !!
  ## if (allTranscripts == TRUE) {
  ##    message("Compound heterozygous analysis cannot handle more than one transcript per gene.\n")
  ##    message("Automatically setting allTranscripts to FALSE")
  ##    allTranscripts <- FALSE
  ## }
   
  pedf <- read.table(ped, header=FALSE, stringsAsFactors=FALSE)
  pedf <- pedf[, 1:6]
  colnames(pedf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Gender", "Phenotype")
  
  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedf$Phenotype == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")

  unaff <- pedf[pedf$Phenotype == 1, ]
  aff <- pedf[pedf$Phenotype == 2, ]
  
  motherID <- unique(aff$MotherID)
  if (length(motherID) > 1) {
    stop("The mother of all the affected individuals has to be the same. Please check out the PED file.")
  } else if (length(motherID) < 1)
    stop("One individual has to be set as mother of the affected individual(s).")
  carrierMother <- pedf[pedf$IndividualID == as.character(motherID), ]
  
  fatherID <- unique(aff$FatherID)
  if (length(fatherID) > 1) {
    stop("The father of all the affected individuals has to be the same. Please check out the PED file.")
  } else if (length(fatherID) < 1)
    stop("One individual has to be set as father of the affected individual(s)")
  carrierFather <- pedf[pedf$IndividualID == as.character(fatherID), ]
  
  mom_comphet <- dad_comphet <- NULL

  annotationCache <- new.env() ## cache annotations when using VariantAnnotation::locateVariants()
  annotated_variants <- VRanges()
  open(vcfFiles[[1]])
  n.var <- 0
  while(nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]], param=svparam))) {

    ## insert an index for each variant in the VCF file, and indicator
    ## values for heterozygous variants in affected and one of the parents only
    info(header(vcf)) <- rbind(info(header(vcf)),
                               DataFrame(Number=rep(1, 3), Type=rep("Integer", 3),
                                         Description=c("Variant index in the VCF file.",
                                                       "Heterozygous in affected and mother only.",
                                                       "Heterozygous in affected and father only."),
                                         row.names=c("VCFIDX", "COMPHETMOTHER", "COMPHETFATHER")))
    info(vcf)$VCFIDX <- (n.var+1):(n.var+nrow(vcf))
    varIDs <- rownames(vcf)

    n.var <- n.var + nrow(vcf)

    ## restrict upfront variants to those in autosomal chromosomes
    ## we subset to the first element of the value returned by seqlevelsStyle()
    ## to deal with cases in which only a subset of chromosomes is contained in
    ## the input VCF (typically for teaching/example/illustration purposes) which
    ## matches more than one chromosome style, or because Ensembl is identical to NCBI for human :\
    autosomalMask <- seqnames(vcf) %in% extractSeqlevelsByGroup(organism(bsgenome),
                                                                seqlevelsStyle(vcf)[1],
                                                                group="auto")
    vcf <- vcf[autosomalMask, ]

    ## heterozygous mask for *all* affected individuals
    affectedMask <- geno(vcf)$GT[, aff$IndividualID, drop=FALSE] == "0/1"
    affectedMask <- apply(affectedMask, 1, all)

    ## heterozygous mask for mother
    motherHetMask <- geno(vcf)$GT[, motherID, drop=FALSE] == "0/1"

    ## heterozygous mask for father
    fatherHetMask <- geno(vcf)$GT[, fatherID, drop=FALSE] == "0/1"

    ## homozygous alternative mask for *any* of the unaffected individuals
    unaffHomMask <- geno(vcf)$GT[, unaff$IndividualID, drop=FALSE] == "1/1"
    unaffHomMask <- apply(unaffHomMask, 1, any)

    ## candidate variants are heterozygous in affected individuals,
    ## heterozygous in one of the parents but not in the other and
    ## are not homozygous in the unaffected individuals
    compHetMotherMask <- affectedMask & motherHetMask & !fatherHetMask & !unaffHomMask
    compHetFatherMask <- affectedMask & !motherHetMask & fatherHetMask & !unaffHomMask
    compHetMask <- compHetMotherMask | compHetFatherMask

    ## annotate separately father and mother contributions to the
    ## compound heterozygous segregation
    info(vcf)$COMPHETMOTHER <- compHetMotherMask + 0L
    info(vcf)$COMPHETFATHER <- compHetFatherMask + 0L

    ## filter out variants that do not segregate as an autosomal recessive heterozygous trait
    ## still a final filter needs to be applied after annotating the variants to genes,
    ## to select those that occur within a common gene where at least one comes from each parent
    vcf <- vcf[compHetMask, ]

    if (any(compHetMask)) {
      ## coerce the VCF object to a VRanges object
      variants <- as(vcf, "VRanges")

      ## since the conversion of VCF to VRanges strips the VCF ID field, let's put it back
      variants$VARID <- varIDs[variants$VCFIDX]

      ## harmonize Seqinfo data between variants, annotations and reference genome
      variants <- .matchSeqinfo(variants, txdb, bsgenome)

      ## annotate variants
      annotated_variants <- c(annotated_variants, annotationEngine(variants, param, annotationCache,
                                                                   BPPARAM=BPPARAM))
    }

    message(sprintf("%d variants processed", n.var))
  }
  close(vcfFiles[[1]])
  
  if (length(annotated_variants) > 0) {
    ## candidate variants should occur within a common gene where at least one comes from each parent
    geneMask <- !is.na(annotated_variants$GENEID)
    annotated_variants <- annotated_variants[geneMask] ## discard variants wo/ gene annotation

    vcfidxbygene <- split(annotated_variants$VCFIDX, annotated_variants$GENEID) ## group variants by gene
    vcfidxbygene <- lapply(vcfidxbygene, ## within each gene, select variants contributed from different parents
                           function(vcfidx, mothervcfidx, fathervcfidx) {
                             motherMask <- vcfidx %in% mothervcfidx
                             fatherMask <- vcfidx %in% fathervcfidx       
                             xorMask <- (motherMask | fatherMask) & !(motherMask & fatherMask)
                             if (sum(xorMask & motherMask) > 1 & sum(xorMask & fatherMask) > 1)
                               vcfidx <- vcfidx[xorMask]
                             else                    ## if there is no variant contributed from
                               vcfidx <- integer(0)  ## any of the parents, discard everything

                             vcfidx
                           }, annotated_variants$VCFIDX[annotated_variants$COMPHETMOTHER == 1L],
                           annotated_variants$VCFIDX[annotated_variants$COMPHETFATHER == 1L])
    elen <- elementNROWS(vcfidxbygene)
    vcfidxbygene <- vcfidxbygene[elen > 1] ## discard variants found alone in a gene
    compHetMask <- annotated_variants$VCFIDX %in% unique(unlist(vcfidxbygene, use.names=FALSE))

    ## select the final set of variant segregating as a compound heterozygous trait
    annotated_variants <- annotated_variants[compHetMask]

  }

  gSO <- annotateSO(annotated_variants, sog(param))
  annotated_variants <- addSOmetadata(annotated_variants)
  locMask <- conMask <- varTypMask <- logical(0)

  if (length(annotated_variants) > 0) {

    locMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$LOCATION)),
                                       levels(annotated_variants$LOCATION)))
    conMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$CONSEQUENCE)),
                                       levels(annotated_variants$CONSEQUENCE)))
    varTypMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$TYPE)),
                                          levels(annotated_variants$TYPE)))
  } else
    warning("No variants segregate following an autosomal recessive heterozygous inheritance model.")

  MAFpopMask <- NA
  if ("MafDb" %in% param$otherAnnotationsClass) {
    ## assume AF columns are those containing AF[A-Z]+ and being of class 'numeric'
    cnAF <- colnames(mcols(annotated_variants))
    colsclasses <- sapply(mcols(annotated_variants), class)
    cnAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

  annoGroups <- list()
  if (!is.null(mcols(mcols(annotated_variants))$TAB)) {
    mask <- !is.na(mcols(mcols(annotated_variants))$TAB)
    annoGroups <- split(colnames(mcols(annotated_variants))[mask],
                      mcols(mcols(annotated_variants))$TAB[mask])
  }
  bsgenome <- get(param$bsgenome)
  genomeDescription <- GenomeDescription(organism=organism(bsgenome),
                                         common_name=commonName(bsgenome),
                                         provider=provider(bsgenome),
                                         provider_version=providerVersion(bsgenome),
                                         release_date=releaseDate(bsgenome),
                                         release_name=releaseName(bsgenome),
                                         seqinfo=seqinfo(bsgenome))

  new("VariantFilteringResults", callObj=callobj, callStr=callstr,
      genomeDescription=genomeDescription, inputParameters=param,
      activeSamples=sampleNames, inheritanceModel="autosomal recessive heterozygous",
      variants=annotated_variants, bamViews=BamViews(), gSO=gSO, filters=filters(param),
      cutoffs=cutoffs(param), annoGroups=annoGroups, dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantTypeMask=varTypMask, locationMask=locMask, consequenceMask=conMask, aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minScore5ss=NA_real_, minScore3ss=NA_real_, minCUFC=0)
})
