
setMethod("autosomalRecessiveHeterozygous", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  callobj <- match.call()
  callstr <- deparse(callobj)
  input_list <- as.list(path(param$vcfFiles))
  ped <- param$pedFilename
  orgdb <- param$orgdb
  txdb <- param$txdb
  snpdb <- param$snpdb
  radicalAAchangeMatrix <- param$radicalAAchangeMatrix
  allTranscripts <- param$allTranscripts
  otherAnnotations <- param$otherAnnotations
  filterTag <- param$filterTag

  genomeVersion <- unique(genome(txdb))
  # now the version of the human genome that will be called will depend on the TxDb version used for the annotation
  
  if (class(txdb) != "TranscriptDb")
    stop("argument 'txdb' should be a 'TranscriptDb' object (see GenomicFeatures package)\n")
  if (allTranscripts == TRUE) {
     message("Compound heterozygous analysis cannot handle more than one transcript per gene.\n")
     message("Automatically setting allTranscripts to FALSE")
     allTranscripts <- FALSE
  }
   
  if (length(input_list) > 1) {
    multiSample <- FALSE
  } else if (length(input_list) == 1) {
    multiSample <- TRUE
  } else {
    stop("A minimum of 1 vcf file have to be provided")
  }
  
  pedf <- read.table(ped, header=F)
  
  aff <- pedf[pedf[, 6] == 2, ]
  aff_ind <- as.character(aff[, 2])
  
  carr1 <- unique(aff[, 4])
  if (length(carr1) > 1) {
    stop("The mother of all the affected individuals has to be the same. Please check out the ped file")
  } else if (length(carr1) < 1) {
    stop("One individual has to be set as mother of the affected individuals")
  }
  row_carr1 <- pedf[pedf[, 2] == as.character(carr1), ]
  
  carr2 <- unique(aff[, 3])
  if (length(carr2) > 1) {
    stop("The father of all the affected individuals has to be the same. Please check out the ped file")
  } else if (length(carr2) < 1) {
    stop("One individual has to be set as father of the affected individuals")
  }
  
  row_carr2 <- pedf[pedf[, 2] == as.character(carr2), ]
  
  
  # analysis
    
  if (multiSample) {
    vcf_vcf1 <- readVcf(unlist(input_list), genomeVersion)
    
    gr_carr1 <- one_ind_ms(vcf_vcf1, "0/1", row_carr1, filterTag)
    if (length(gr_carr1) < 1) {
      stop("Unable to find the data about the mother inside de VCF file. Please check that the id is the same for both files")
    }
    gr_carr2 <- one_ind_ms(vcf_vcf1, "0/1", row_carr2, filterTag)
    if (length(gr_carr2) < 1) {
      stop("Unable to find the data about the father inside de VCF file. Please check that the id is the same for both files")
    }
    
    
    affected <- switch (nrow(aff)+1,
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
    homodiscard1 <- rowData(vcf_vcf1[geno(vcf_vcf1)$GT[, row_carr1[, 2]] == "1/1", ])
    homodiscard2 <- rowData(vcf_vcf1[geno(vcf_vcf1)$GT[, row_carr2[, 2]] == "1/1", ])
    
    mom_comphet_n <- carr1[!names(carr1) %in% names(homodiscard2), ]
    dad_comphet_n <- carr2[!names(carr2) %in% names(homodiscard1), ]
    
    mom_comphet <- matchChromosomeNames(mom_comphet_n, txdb)
    dad_comphet <- matchChromosomeNames(dad_comphet_n, txdb)
    
    
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
  
    vcf_carrier1 <- readVcf(unlist(input_list_mother), genomeVersion)
    vcf_carrier2 <- readVcf(unlist(input_list_father), genomeVersion)
    
    affected <- switch (length(input_list_aff)+1,
                        stop("No affected individuals detected. Something might be wrong with the .ped file..."),
                        one_ind_us(input_list_aff, "0/1", filterTag, genomeVersion),
                        two_ind_us(input_list_aff, "0/1", filterTag, genomeVersion),
                        three_ind_us(input_list_aff, "0/1", filterTag, genomeVersion),
                        four_ind_us(input_list_aff, "0/1", filterTag, genomeVersion),
                        five_ind_us(input_list_aff, "0/1", filterTag, genomeVersion))

  
    
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
    homodiscard1 <- rowData(vcf_carrier1[geno(vcf_carrier1)$GT[, 1] == "1/1", ])
    homodiscard2 <- rowData(vcf_carrier2[geno(vcf_carrier2)$GT[, 1] == "1/1", ])
    
    mom_comphet_n <- carr1[!names(carr1) %in% names(homodiscard2), ]
    dad_comphet_n <- carr2[!names(carr2) %in% names(homodiscard1), ]
    
    mom_comphet <- matchChromosomeNames(mom_comphet_n, txdb)
    dad_comphet <- matchChromosomeNames(dad_comphet_n, txdb)
    
  }
  

  affected <- matchChromosomeNames(affected, txdb)

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

  cnAF <- colnames(mcols(parents_contrib_sorted))
  cnAF <- cnAF[grep("AF", cnAF)]
  MAFpopMask <- rep(TRUE, length(cnAF))
  names(MAFpopMask) <- cnAF

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="autosomal recessive heterozygous", variants=parents_contrib_sorted,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any", aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
