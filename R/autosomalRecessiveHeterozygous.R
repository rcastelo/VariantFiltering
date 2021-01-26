
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

  if (any(!exists(as.character(substitute(BPPARAM)))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))

  if (length(vcfFiles) > 1)
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with software such as vcftools or GATK, or do the variant calling simultaneously on all samples, or proceed analyzing each file separately.")
  else if (length(vcfFiles) < 1)
    stop("A minimum of 1 VCF file has to be provided.")
  
  if (is.na(ped))
    stop("Please specify a PED file name when building the parameter object.")

  pedDf <- .readPEDfile(ped)

  ## REVISE THIS !!
  ## if (allTranscripts == TRUE) {
  ##    message("Compound heterozygous analysis cannot handle more than one transcript per gene.\n")
  ##    message("Automatically setting allTranscripts to FALSE")
  ##    allTranscripts <- FALSE
  ## }
   
  ## assuming Phenotype == 2 means affected and Phenotype == 1 means unaffected
  if (sum(pedDf$Phenotype == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")

  unaff <- pedDf[pedDf$Phenotype == 1, ]
  aff <- pedDf[pedDf$Phenotype == 2, ]
  
  motherID <- unique(aff$MotherID)
  if (length(motherID) > 1) {
    stop("The mother of all the affected individuals has to be the same. Please check out the PED file.")
  } else if (length(motherID) < 1)
    stop("One individual has to be set as mother of the affected individual(s).")
  carrierMother <- pedDf[pedDf$IndividualID == as.character(motherID), ]
  
  fatherID <- unique(aff$FatherID)
  if (length(fatherID) > 1) {
    stop("The father of all the affected individuals has to be the same. Please check out the PED file.")
  } else if (length(fatherID) < 1)
    stop("One individual has to be set as father of the affected individual(s)")
  carrierFather <- pedDf[pedDf$IndividualID == as.character(fatherID), ]
  
  mom_comphet <- dad_comphet <- NULL

  annotationCache <- new.env() ## cache annotations when using VariantAnnotation::locateVariants()
  annotated_variants <- VRanges()
  metadata(mcols(annotated_variants)) <- list(cutoffs=CutoffsList(), sortings=CutoffsList())

  open(vcfFiles[[1]])
  n.var <- 0
  flag <- TRUE
  while(flag && nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]], param=svparam))) {

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
    vcf <- vcf[autosomalMask, , drop=FALSE]

    gt <- geno(vcf)$GT

    ## further restrict affected and unaffected individuals to
    ## those who have been genotyped
    gtind <- colnames(gt)
    unaff <- unaff[unaff$IndividualID %in% gtind, , drop=FALSE]
    aff <- aff[aff$IndividualID %in% gtind, , drop=FALSE]
    if (nrow(aff) == 0)
      stop("No affected individuals have genotypes.")
    if (!fatherID %in% gtind)
      stop("Father is not genotyped.")
    if (!motherID %in% gtind)
      stop("Mother is not genotyped.")

    ## heterozygous mask for *all* affected individuals
    affectedMask <- gt[, aff$IndividualID, drop=FALSE] == "0/1" ||
                    gt[, aff$IndividualID, drop=FALSE] == "0|1" ||
                    gt[, aff$IndividualID, drop=FALSE] == "1|0"
    affectedMask <- apply(affectedMask, 1, all)

    ## heterozygous mask for mother
    motherHetMask <- gt[, motherID, drop=FALSE] == "0/1" ||
                     gt[, motherID, drop=FALSE] == "0|1" ||
                     gt[, motherID, drop=FALSE] == "1|0"

    ## heterozygous mask for father
    fatherHetMask <- gt[, fatherID, drop=FALSE] == "0/1" ||
                     gt[, fatherID, drop=FALSE] == "0|1" ||
                     gt[, fatherID, drop=FALSE] == "1|0"

    ## homozygous alternative mask for *any* of the unaffected individuals
    unaffHomMask <- gt[, unaff$IndividualID, drop=FALSE] == "1/1" ||
                    gt[, unaff$IndividualID, drop=FALSE] == "1|1"
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
      if (length(annotated_variants) > 0)
        annotated_variants <- c(annotated_variants, annotationEngine(variants, param, annotationCache,
                                                                     BPPARAM=BPPARAM))
      else
        annotated_variants <- annotationEngine(variants, param, annotationCache,
                                               BPPARAM=BPPARAM)
    }

    if (length(vcfWhich(svparam)) > 0) ## temporary fix to keep this looping
      flag <- FALSE                    ## structure with access through genomic ranges

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

  if (length(annotated_variants) == 0)
    warning("No variants segregate following an autosomal recessive heterozygous inheritance model.")

  annoGroups <- list()
  if (!is.null(mcols(mcols(annotated_variants))$TAB)) {
    mask <- !is.na(mcols(mcols(annotated_variants))$TAB)
    annoGroups <- split(colnames(mcols(annotated_variants))[mask],
                      mcols(mcols(annotated_variants))$TAB[mask])
  }

  ## add functional annotation filters generated by the annotation engine
  funFilters <- FilterRules(lapply(metadata(mcols(annotated_variants))$filters,
                                   function(f) { environment(f) <- baseenv() ; f}))
  active(funFilters) <- FALSE ## by default, functional annotation filters are inactive
  flt <- c(filters(param), funFilters)
  fltMd <- rbind(filtersMetadata(param),
                 DataFrame(Description=sapply(metadata(mcols(annotated_variants))$filters, 
                                              attr, "description"),
                           AnnoGroup=sapply(metadata(mcols(annotated_variants))$filters,
                                            attr, "TAB")))
  cutoffs <- metadata(mcols(annotated_variants))$cutoffs
  sortings <- metadata(mcols(annotated_variants))$sortings
  bsgenome <- get(param$bsgenome)

  new("VariantFilteringResults", callObj=callobj, callStr=callstr,
      genomeDescription=as(bsgenome, "GenomeDescription"), inputParameters=param,
      activeSamples=sampleNames, inheritanceModel="autosomal recessive heterozygous",
      variants=annotated_variants, bamViews=BamViews(), gSO=gSO, filters=flt,
      filtersMetadata=fltMd, cutoffs=cutoffs, sortings=sortings, annoGroups=annoGroups,
      minScore5ss=NA_real_, minScore3ss=NA_real_)
})
