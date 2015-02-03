setMethod("unrelatedIndividuals", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  callobj <- match.call()
  callstr <- deparse(callobj)
  vcfFiles <- param$vcfFiles
  seqInfos <- param$seqInfos
  txdb <- param$txdb
  bsgenome <- param$bsgenome
  
  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))
  
  if (length(vcfFiles) > 1) {
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with vcftools, do the variant calling simultaneously on all samples, or proceed analysing each file separately.")
  } else if (length(vcfFiles) < 1)
    stop("A minimum of 1 vcf file has to be provided")

  
  annotated_variants <- GRanges()
  open(vcfFiles[[1]])
  n.var <- 0
  while(nrow(vcf <- readVcf(vcfFiles[[1]], genome=seqInfos[[1]]))) {

    variants <- rowData(vcf)

    variants <- .matchSeqinfo(variants, txdb, bsgenome)

    variants <- annotationEngine(variants, param, BPPARAM=BPPARAM)

    homalt <- geno(vcf)$GT == "1/1" | geno(vcf)$GT == "1|1"
    colnames(homalt) <- paste0("homALT_", colnames(homalt))
    homalt <- homalt[variants$IDX, , drop=FALSE]
    rownames(homalt) <- NULL
    mcols(variants) <- cbind(mcols(variants), DataFrame(homalt))
    homalt <- DataFrame(HomozygousAlt=CharacterList(apply(homalt, 1,
                                                          function(mask, snames)
                                                            snames[mask],
                                                          rownames(colData(vcf)))))
    mcols(variants) <- cbind(mcols(variants), DataFrame(homalt))

    het <- geno(vcf)$GT == "0/1" | geno(vcf)$GT == "0|1" |
           geno(vcf)$GT == "1/0" | geno(vcf)$GT == "1|0"
    colnames(het) <- paste0("het_", colnames(het))
    het <- het[variants$IDX, , drop=FALSE]
    rownames(het) <- NULL
    mcols(variants) <- cbind(mcols(variants), DataFrame(het))
    het <- DataFrame(Heterozygous=CharacterList(apply(het, 1,
                                                      function(mask, snames)
                                                         snames[mask],
                                                       rownames(colData(vcf)))))
    mcols(variants) <- cbind(mcols(variants), DataFrame(het))

    homref <- geno(vcf)$GT == "0/0" | geno(vcf)$GT == "0|0"
    colnames(homref) <- paste0("homREF_", colnames(homref))
    homref <- homref[variants$IDX, , drop=FALSE]
    rownames(homref) <- NULL
    mcols(variants) <- cbind(mcols(variants), DataFrame(homref))
    homref <- DataFrame(HomozygousRef=CharacterList(apply(homref, 1,
                                                          function(mask, snames)
                                                            snames[mask],
                                                          rownames(colData(vcf)))))
    mcols(variants) <- cbind(mcols(variants), DataFrame(homref))

    annotated_variants <- c(annotated_variants, variants)

    n.var <- n.var + nrow(vcf)
    message(sprintf("%d variants processed", n.var))
  }
  close(vcfFiles[[1]])

  MAFpopMask <- NA
  if ("MafDb" %in% sapply(param$otherAnnotations, class)) {
    ## assume AF columns are those containing AF[A-Z]+ and being of class 'numeric'
    cnAF <- colnames(mcols(annotated_variants))
    colsclasses <- sapply(mcols(annotated_variants), class)
    nAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

  new("VariantFilteringResultsUI", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="unrelated individuals", variants=annotated_variants,
      indselected=NA_character_, selectgene=NA_character_, dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantType="Any", aaChangeType="Any", MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
