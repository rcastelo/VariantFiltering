setMethod("unrelatedIndividuals", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  callobj <- match.call()
  callstr <- deparse(callobj)
  input_list <- as.list(path(param$vcfFiles))
  genomeInfo <- param$seqInfos[[1]]
  txdb <- param$txdb
  
  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))
  
  if (class(txdb) != "TxDb")
    stop("argument 'txdb' should be a 'TxDb' object (see GenomicFeatures package)\n")
  
  message("Reading input VCF file into main memory.")
  vcf1 <- readVcf(unlist(input_list), genomeInfo)
  unrelated <- rowData(vcf1)

  unrelated <- matchChromosomes(unrelated, txdb)

  ##########################
  ##                      ##
  ##      ANNOTATION      ##
  ##                      ##
  ##########################
  
  unrelated_annotated <- annotationEngine(unrelated, param, BPPARAM=BPPARAM)

  homalt <- geno(vcf1)$GT == "1/1" | geno(vcf1)$GT == "1|1"
  colnames(homalt) <- paste0("homALT_", colnames(homalt))
  homalt <- homalt[unrelated_annotated$IDX, , drop=FALSE]
  rownames(homalt) <- NULL
  mcols(unrelated_annotated) <- cbind(mcols(unrelated_annotated), DataFrame(homalt))
  homalt <- DataFrame(HomozygousAlt=CharacterList(apply(homalt, 1,
                                                        function(mask, snames)
                                                          snames[mask],
                                                        rownames(colData(vcf1)))))
  mcols(unrelated_annotated) <- cbind(mcols(unrelated_annotated), DataFrame(homalt))

  het <- geno(vcf1)$GT == "0/1" | geno(vcf1)$GT == "0|1" |
         geno(vcf1)$GT == "1/0" | geno(vcf1)$GT == "1|0"
  colnames(het) <- paste0("het_", colnames(het))
  het <- het[unrelated_annotated$IDX, , drop=FALSE]
  rownames(het) <- NULL
  mcols(unrelated_annotated) <- cbind(mcols(unrelated_annotated), DataFrame(het))
  het <- DataFrame(Heterozygous=CharacterList(apply(het, 1,
                                                    function(mask, snames)
                                                       snames[mask],
                                                     rownames(colData(vcf1)))))
  mcols(unrelated_annotated) <- cbind(mcols(unrelated_annotated), DataFrame(het))

  homref <- geno(vcf1)$GT == "0/0" | geno(vcf1)$GT == "0|0"
  colnames(homref) <- paste0("homREF_", colnames(homref))
  homref <- homref[unrelated_annotated$IDX, , drop=FALSE]
  rownames(homref) <- NULL
  mcols(unrelated_annotated) <- cbind(mcols(unrelated_annotated), DataFrame(homref))
  homref <- DataFrame(HomozygousRef=CharacterList(apply(homref, 1,
                                                        function(mask, snames)
                                                          snames[mask],
                                                        rownames(colData(vcf1)))))
  mcols(unrelated_annotated) <- cbind(mcols(unrelated_annotated), DataFrame(homref))

  ##########################
  ##                      ##
  ## BUILD RESULTS OBJECT ##
  ##                      ##
  ##########################

  cnAF <- colnames(mcols(unrelated_annotated))
  cnAF <- cnAF[grep("AF", cnAF)]
  MAFpopMask <- rep(TRUE, length(cnAF))
  names(MAFpopMask) <- cnAF

  new("VariantFilteringResultsUI", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="unrelated individuals", variants=unrelated_annotated,
      indselected=NA_character_, selectgene=NA_character_, dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantType="Any", aaChangeType="Any", MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
