setMethod("xLinked", signature(param="VariantFilteringParam"),
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
 
  if (length(input_list) > 1) {
    multiSample <- FALSE
  } else if (length(input_list) == 1) {
    multiSample <- TRUE
  } else {
    stop("A minimum of 1 vcf file have to be provided")
  }
  
  pedf <- read.table(ped, header=F)
  
  unaff_males <- pedf[pedf[, 6] == 1 & pedf[, 5] == 1, ]
  carrier_female <-  pedf[pedf[, 6] == 1 & pedf[, 5] == 2, ]
  
  aff_males <- pedf[pedf[, 6] == 2 & pedf[, 5] == 1, ]
  
  
  if (multiSample) {
    vcf_vcf1 <- readVcf(unlist(input_list), genomeVersion) 
    
    carrierFemales <- switch (nrow(carrier_female),
                              one_ind_ms(vcf_vcf1, "0/1", carrier_female, filterTag),
                              two_ind_ms(vcf_vcf1, "0/1", carrier_female, filterTag),
                              three_ind_ms(vcf_vcf1, "0/1", carrier_female, filterTag),
                              four_ind_ms(vcf_vcf1, "0/1", carrier_female, filterTag),
                              five_ind_ms(vcf_vcf1, "0/1", carrier_female, filterTag))
    
    unaff_males <- switch (nrow(unaff_males),
                          one_ind_ms(vcf_vcf1, "0/0", unaff_males, filterTag),
                          two_ind_ms(vcf_vcf1, "0/0", unaff_males, filterTag),
                          three_ind_ms(vcf_vcf1, "0/0", unaff_males, filterTag),
                          four_ind_ms(vcf_vcf1, "0/0", unaff_males, filterTag),
                          five_ind_ms(vcf_vcf1, "0/0", unaff_males, filterTag))
    
    affected <- switch (nrow(aff_males)+1,
                        stop("No affected males detected. Something might be wrong with the .ped file..."),
                        one_ind_ms(vcf_vcf1, "1/1", aff_males, filterTag),
                        two_ind_ms(vcf_vcf1, "1/1", aff_males, filterTag),
                        three_ind_ms(vcf_vcf1, "1/1", aff_males, filterTag),
                        four_ind_ms(vcf_vcf1, "1/1", aff_males, filterTag),
                        five_ind_ms(vcf_vcf1, "1/1", aff_males, filterTag))

    
    if (length(unaff_males) < 1) {
      affected_filtered <- affected
    } else {
      
      realaffected <- sharedVariants(affected, unaff_males)
      affected_filtered <- affected[realaffected]
    }
    
  } else {
    
    input_list_ind <- sapply(input_list, function(x) gsub("\\.vcf$", "", gsub("\\.vcf\\.gz$", "", x, ignore.case=TRUE), ignore.case=TRUE))
    
    input_list_unaff_males_vector <- c()
    input_list_carrier_female_vector <- c()
    input_list_aff_males_vector <- c()
    
    for (i in 1:length(input_list_ind)) {
      if (input_list_ind[i] %in% unaff_males[, 2]) {
        input_list_unaff_males_vector <-  c(input_list_unaff_males_vector, i)
      } else if (input_list_ind[i]%in% carrier_female[, 2]) {
        input_list_carrier_female_vector <-  c(input_list_carrier_female_vector, i)
      } else if (input_list_ind[i]%in% aff_males[, 2]) {
        input_list_aff_males_vector <- c(input_list_aff_males_vector, i)
      }
    }
    
    input_list_unaff_males <- input_list[input_list_unaff_males_vector]
    input_list_carrier_female <- input_list[input_list_carrier_female_vector]
    input_list_aff_males <- input_list[input_list_aff_males_vector]
    
    carrierFemales <- switch (length(input_list_carrier_female),
                              one_ind_us(input_list_carrier_female, "0/1", filterTag, genomeVersion),
                              two_ind_us(input_list_carrier_female, "0/1", filterTag, genomeVersion),
                              three_ind_us(input_list_carrier_female, "0/1", filterTag, genomeVersion),
                              four_ind_us(input_list_carrier_female, "0/1", filterTag, genomeVersion),
                              five_ind_us(input_list_carrier_female, "0/1", filterTag, genomeVersion))
    
    affected <- switch (length(input_list_aff_males)+1,
                        stop("No affected males detected. Something might be wrong with the .ped file..."),
                        one_ind_us(input_list_aff_males, "1/1", filterTag, genomeVersion),
                        two_ind_us(input_list_aff_males, "1/1", filterTag, genomeVersion),
                        three_ind_us(input_list_aff_males, "1/1", filterTag, genomeVersion),
                        four_ind_us(input_list_aff_males, "1/1", filterTag, genomeVersion),
                        five_ind_us(input_list_aff_males, "1/1", filterTag, genomeVersion))
    
    if (length(input_list_unaff_males) < 1) {
      affected_filtered <- affected
    } else {
      affected_filtered <- pullout_shared_us(affected, input_list_unaff_males, filterTag, genomeVersion)
    }
  }   
  
  if (length(carrierFemales) < 1) {
    xl_all <- affected_filtered
  } else {
    realcommonXlinked <- sharedVariants(affected_filtered, carrierFemales)
    xl_all <- affected_filtered[realcommonXlinked]
  }
  
  xl_allchr <- matchChromosomeNames(xl_all, txdb)
  
  Xlist <- which(seqnames(xl_allchr) == "chrX")
  xlinked <- xl_allchr[Xlist]

  ##########################
  ##                      ##
  ##      ANNOTATION      ##
  ##                      ##
  ##########################
  
  xlinked_annotated <- annotationEngine(xlinked, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
                                        radicalAAchangeMatrix=radicalAAchangeMatrix,
                                        otherAnnotations=otherAnnotations,
                                        allTranscripts=allTranscripts, BPPARAM=BPPARAM)

  ##########################
  ##                      ##
  ## BUILD RESULTS OBJECT ##
  ##                      ##
  ##########################

  cnAF <- colnames(mcols(xlinked_annotated))
  cnAF <- cnAF[grep("AF", cnAF)]
  MAFpopMask <- rep(TRUE, length(cnAF))
  names(MAFpopMask) <- cnAF

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="X-linked", variants=xlinked_annotated,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any", aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
