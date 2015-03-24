
setMethod("autosomalRecessiveHeterozygous", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  stop("this is being updated and cannot be used at the moment.")

  ## store call for reproducing it later
  callobj <- match.call()
  callstr <- gsub(".local", "autosomalRecessiveHeterozygous", deparse(callobj))

  ## fetch necessary parameters
  vcfFiles <- param$vcfFiles
  ped <- param$pedFilename
  seqInfos <- param$seqInfos
  txdb <- param$txdb
  bsgenome <- param$bsgenome
  sampleNames <- param$sampleNames

  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))

  if (length(vcfFiles) > 1)
    stop("More than one input VCF file is currently not supported. Please either merge the VCF files into a single one with software such as vcftools or GATK, or do the variant calling simultaneously on all samples, or proceed analyzing each file separately.")
  else if (length(vcfFiles) < 1)
    stop("A minimum of 1 VCF file has to be provided.")
  
  ## if (allTranscripts == TRUE) {
  ##    message("Compound heterozygous analysis cannot handle more than one transcript per gene.\n")
  ##    message("Automatically setting allTranscripts to FALSE")
  ##    allTranscripts <- FALSE
  ## }
   
  pedf <- read.table(ped, header=FALSE, stringsAsFactors=FALSE)
  pedf <- pedf[, 1:6]
  colnames(pedf) <- c("FamilyID", "IndividualID", "FatherID", "MotherID", "Gender", "Phenotype")
  
  if (sum(pedf$Phenotype == 2) < 1)
    stop("No affected individuals detected. Something is wrong with the PED file.")

  aff <- pedf[pedf$Phenotype == 2, ]
  
  carr1 <- unique(aff$MotherID)
  if (length(carr1) > 1) {
    stop("The mother of all the affected individuals has to be the same. Please check out the PED file.")
  } else if (length(carr1) < 1)
    stop("One individual has to be set as mother of the affected individual(s).")
  row_carr1 <- pedf[pedf$IndividualID == as.character(carr1), ]
  
  carr2 <- unique(aff$FatherID)
  if (length(carr2) > 1) {
    stop("The father of all the affected individuals has to be the same. Please check out the PED file.")
  } else if (length(carr2) < 1)
    stop("One individual has to be set as father of the affected individual(s)")
  row_carr2 <- pedf[pedf$IndividualID == as.character(carr2), ]
  
  ## CONTINUE HERE !!!!!
  # analysis
    
  mom_comphet <- dad_comphet <- NULL

  if (multiSample) {
    vcf_vcf1 <- readVcf(unlist(input_list), genomeInfo)
    
    gr_carr1 <- one_ind_ms(vcf_vcf1, "0/1", row_carr1, filterTag)
    if (length(gr_carr1) < 1) {
      stop("Unable to find the data about the mother inside de VCF file. Please check that the id is the same for both files")
    }
    gr_carr2 <- one_ind_ms(vcf_vcf1, "0/1", row_carr2, filterTag)
    if (length(gr_carr2) < 1) {
      stop("Unable to find the data about the father inside de VCF file. Please check that the id is the same for both files")
    }
    
    
    affected <- switch(nrow(aff)+1,
                       stop("No affected individuals detected. Something might be wrong with the .ped file..."),
                       one_ind_ms(vcf_vcf1, "0/1", aff, filterTag),
                       two_ind_ms(vcf_vcf1, "0/1", aff, filterTag),
                       three_ind_ms(vcf_vcf1, "0/1", aff, filterTag),
                       four_ind_ms(vcf_vcf1, "0/1", aff, filterTag),
                       five_ind_ms(vcf_vcf1, "0/1", aff, filterTag))
    
    realcommongen12 <- sharedVariants(gr_carr1, gr_carr2)
    realcommongen21 <- sharedVariants(gr_carr2, gr_carr1)
    
    # i have the shared ones. I want all the others.
    carr1 <- gr_carr1[-realcommongen12]
    carr2 <- gr_carr2[-realcommongen21]
    
    # now, a filter to discard same changes of one parent wich appears in the other AS HOMOZYGOUS, therefore, are not responsible
    # for the phenotype (this parent would present an het change in one candidate and a hom change in the other one, passing the first filters and also being compound heterozygous himself)
    homodiscard1 <- rowRanges(vcf_vcf1[geno(vcf_vcf1)$GT[, row_carr1[, 2]] == "1/1", ])
    homodiscard2 <- rowRanges(vcf_vcf1[geno(vcf_vcf1)$GT[, row_carr2[, 2]] == "1/1", ])
    
    mom_comphet <- carr1[!names(carr1) %in% names(homodiscard2), ]
    dad_comphet <- carr2[!names(carr2) %in% names(homodiscard1), ]
    
    mom_comphet <- matchChromosomes(mom_comphet, txdb)
    dad_comphet <- matchChromosomes(dad_comphet, txdb)
    
  } else {
    
    input_list_ind <- sapply(input_list, function(x) gsub("\\.vcf$", "", gsub("\\.vcf\\.gz$", "", x, ignore.case=TRUE), ignore.case=TRUE))
    
    input_list_father_vector <- c()
    input_list_mother_vector <- c()
    input_list_aff_vector <- c()
    
    for (i in 1:length(input_list_ind)) {
      if (input_list_ind[i] %in% aff[, 3]) {
        input_list_father_vector <-  i
      } else if (input_list_ind[i]%in% aff[, 4]) {
        input_list_mother_vector <-  i
      } else if (input_list_ind[i]%in% aff[, 2]) {
        input_list_aff_vector <- c(input_list_aff_vector, i)
      }
    }
    
    input_list_father <- input_list[input_list_father_vector]
    input_list_mother <- input_list[input_list_mother_vector]
    input_list_aff <- input_list[input_list_aff_vector] 
  
    vcf_carrier1 <- readVcf(unlist(input_list_mother), genomeInfo)
    vcf_carrier2 <- readVcf(unlist(input_list_father), genomeInfo)
    
    affected <- switch(length(input_list_aff)+1,
                       stop("No affected individuals detected. Something might be wrong with the .ped file..."),
                       one_ind_us(input_list_aff, "0/1", filterTag, genomeInfo),
                       two_ind_us(input_list_aff, "0/1", filterTag, genomeInfo),
                       three_ind_us(input_list_aff, "0/1", filterTag, genomeInfo),
                       four_ind_us(input_list_aff, "0/1", filterTag, genomeInfo),
                       five_ind_us(input_list_aff, "0/1", filterTag, genomeInfo))
    
    gr_carr1 <- vcf2GR(vcf_carrier1, "0/1", filterTag)
    gr_carr2 <- vcf2GR(vcf_carrier2, "0/1", filterTag)
    
  
    # now we want to discard shared changes on the same position, because it will lead to a compound het
    # in one of the healthy parents
    
    # in few steps we verificate that exactly same changes are present in both individuals, checking same alleles and so on
    # and then, we discard them
    
    
    realcommongen12 <- sharedVariants(gr_carr1, gr_carr2)
    realcommongen21 <- sharedVariants(gr_carr2, gr_carr1)
    
    # i have the shared ones. Now i have to select all the others.
    carr1 <- gr_carr1[-realcommongen12]
    carr2 <- gr_carr2[-realcommongen21]
    
    # now, a filter to discard same changes of one parent wich appears in the other AS HOMOZYGOUS, therefore, are not responsible
    # for the phenotype (this parent would present an het change in one candidate and a hom change in the other one, passing the first filters and also being compound heterozygous himself)
    homodiscard1 <- rowRanges(vcf_carrier1[geno(vcf_carrier1)$GT[, 1] == "1/1", ])
    homodiscard2 <- rowRanges(vcf_carrier2[geno(vcf_carrier2)$GT[, 1] == "1/1", ])
    
    mom_comphet <- carr1[!names(carr1) %in% names(homodiscard2), ]
    dad_comphet <- carr2[!names(carr2) %in% names(homodiscard1), ]
    
    mom_comphet <- matchChromosomes(mom_comphet, txdb)
    dad_comphet <- matchChromosomes(dad_comphet, txdb)
  }
  
  affected <- .matchSeqinfo(affected, txdb, bsgenome)

  ##########################
  ##                      ##
  ##      ANNOTATION      ##
  ##                      ##
  ##########################
  
  affected_annotated <- annotationEngine(affected, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
                                         radicalAAchangeMatrix=radicalAAchangeMatrix,
                                         otherAnnotations=otherAnnotations,
                                         allTranscripts=allTranscripts, BPPARAM=BPPARAM)

  afflist <- affected_annotated[which(duplicated(affected_annotated$GENE))]
  afflist <- afflist[!is.na(afflist$GENE)]
  
  realcommonAffMom <- sharedVariants(afflist, mom_comphet)
  affMom_comphet <- afflist[realcommonAffMom]

  realcommonAffDad <- sharedVariants(afflist, dad_comphet)
  affDad_comphet <- afflist[realcommonAffDad]

  # now only lefts checking if there is at least one change per gene per parent 
  
  mom_contrib <- affMom_comphet[affMom_comphet$GENEID %in% affDad_comphet$GENEID]
  dad_contrib <- affDad_comphet[affDad_comphet$GENEID %in% affMom_comphet$GENEID]

  mcols(mom_contrib) <- cbind(mcols(mom_contrib), DataFrame(SOURCE=rep(as.character(aff[, 4]), length(mom_contrib))))
  mcols(dad_contrib) <- cbind(mcols(dad_contrib), DataFrame(SOURCE=rep(as.character(aff[, 3]), length(dad_contrib))))

 parents_contrib <- c(mom_contrib, dad_contrib)

 parents_contrib_sorted <- parents_contrib[order(parents_contrib$GENE)]
  
  ##########################
  ##                      ##
  ## BUILD RESULTS OBJECT ##
  ##                      ##
  ##########################

  locMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$LOCATION)),
                                     levels(annotated_variants$LOCATION)))
  conMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$CONSEQUENCE)),
                                     levels(annotated_variants$CONSEQUENCE)))
  varTypMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$TYPE)),
                                        levels(annotated_variants$TYPE)))
  MAFpopMask <- NA
  if ("MafDb" %in% sapply(param$otherAnnotations, class)) {
    ## assume AF columns are those containing AF[A-Z]+ and being of class 'numeric'
    cnAF <- colnames(mcols(annotated_variants))
    colsclasses <- sapply(mcols(annotated_variants), class)
    cnAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

  gSO <- sequence_variant.gSOXP
  nodeDataDefaults(gSO, "varIdx") <- integer(0)

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      activeSamples=sampleNames, inheritanceModel="autosomal recessive heterozygous",
      variants=parents_contrib_sorted, bamViews=BamViews(), gSO=gSO, dbSNPflag=NA_character_, OMIMflag=NA_character_,
      variantTypeMask=varTypMask, locationMask=locMask, consequenceMask=conMask, aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_, minCUFC=0)
})
