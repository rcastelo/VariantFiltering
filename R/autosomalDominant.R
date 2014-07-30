setMethod("autosomalDominant", signature(param="VariantFilteringParam"),
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

  if (class(txdb) != "TxDb")
    stop("argument 'txdb' should be a 'TxDb' object (see GenomicFeatures package)\n")
 
  if (length(input_list) > 1) {
    multiSample <- FALSE
  } else if (length(input_list) == 1) {
    multiSample <- TRUE
  } else {
    stop("A minimum of 1 vcf file have to be provided")
  }
  
  pedf <- read.table(ped, header=F)
  
  unaff <- pedf[pedf[, 6] == 1, ]
  unaff_ind <- as.character(unaff[, 2])
  
  aff <- pedf[pedf[, 6] == 2, ]
  aff_ind <- as.character(aff[, 2])
  
  
  if (multiSample) {
    vcf_vcf1 <- readVcf(unlist(input_list), genomeVersion)
    
    unaffected <- switch (nrow(unaff),
                        one_ind_ms(vcf_vcf1, "0/0", unaff, filterTag),
                        two_ind_ms(vcf_vcf1, "0/0", unaff, filterTag),
                        three_ind_ms(vcf_vcf1, "0/0", unaff, filterTag),
                        four_ind_ms(vcf_vcf1, "0/0", unaff, filterTag),
                        five_ind_ms(vcf_vcf1, "0/0", unaff, filterTag))
    
    affected <- switch (nrow(aff)+1,
                        stop("No affected individuals detected. Something might be wrong with the .ped file..."),
                        one_ind_ms_2opt(vcf_vcf1, "0/1", "1/1", aff, filterTag),
                        two_ind_ms_2opt(vcf_vcf1, "0/1", "1/1", aff, filterTag),
                        three_ind_ms_2opt(vcf_vcf1, "0/1", "1/1", aff, filterTag),
                        four_ind_ms_2opt(vcf_vcf1, "0/1", "1/1", aff, filterTag),
                        five_ind_ms_2opt(vcf_vcf1, "0/1", "1/1", aff, filterTag))
    
    if (length(unaffected) < 1) {
      dominant <- affected
    } else {
      realcommondominant <- sharedVariants(affected, unaffected)
      dominant <- affected[realcommondominant]
    }
    
  } else {
    
    input_list_ind <- sapply(input_list, function(x) gsub("\\.vcf$", "", gsub("\\.vcf\\.gz$", "", x, ignore.case=TRUE), ignore.case=TRUE))
    
    input_list_unaff_vector <- c()
    input_list_aff_vector <- c()
    
    for (i in 1:length(input_list_ind)) {
      if (input_list_ind[i] %in% unaff[, 2]) {
        input_list_unaff_vector <- c(input_list_unaff_vector, i)
      } else if (input_list_ind[i]%in% aff[, 2]) {
        input_list_aff_vector <- c(input_list_aff_vector, i)
      }
    }
    
    input_list_unaff <- input_list[input_list_unaff_vector]
    input_list_aff <- input_list[input_list_aff_vector]
    
    affected <- switch (length(input_list_aff)+1,
                        stop("No affected individuals detected. Something might be wrong with the .ped file..."),
                        one_ind_us_2opt(input_list_aff, "0/1", "1/1", filterTag, genomeVersion),
                        two_ind_us_2opt(input_list_aff, "0/1", "1/1", filterTag, genomeVersion),
                        three_ind_us_2opt(input_list_aff, "0/1", "1/1", filterTag, genomeVersion),
                        four_ind_us_2opt(input_list_aff, "0/1", "1/1", filterTag, genomeVersion),
                        five_ind_us_2opt(input_list_aff, "0/1", "1/1", filterTag, genomeVersion))
    
    if (length(input_list_unaff) < 1) {
      dominant <- affected
    } else {
      dominant <- pullout_shared_us(affected, input_list_unaff, filterTag, genomeVersion)
    }
  }
  
  dominant <- matchChromosomeNames(dominant, txdb)
  
  ##########################
  ##                      ##
  ##      ANNOTATION      ##
  ##                      ##
  ##########################
  
  dominant_annotated <- annotationEngine(dominant, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
                                         radicalAAchangeMatrix=radicalAAchangeMatrix,
                                         otherAnnotations=otherAnnotations,
                                         allTranscripts=allTranscripts, BPPARAM=BPPARAM)

  ##########################
  ##                      ##
  ## BUILD RESULTS OBJECT ##
  ##                      ##
  ##########################

  cnAF <- colnames(mcols(dominant_annotated))
  cnAF <- cnAF[grep("AF", cnAF)]
  MAFpopMask <- rep(TRUE, length(cnAF))
  names(MAFpopMask) <- cnAF

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="autosomal dominant", variants=dominant_annotated,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any", aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
