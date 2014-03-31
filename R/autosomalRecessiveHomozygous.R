setMethod("autosomalRecessiveHomozygous", signature(param="VariantFilteringParam"),
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
    stop("A minimum of 1 vcf file has to be provided")
  }
  
  pedf <- read.table(ped, header=F)
  
  unaff <- pedf[pedf[, 6] == 1, ]
  unaff_ind <- as.character(unaff[, 2])

  aff <- pedf[pedf[, 6] == 2, ]
  aff_ind <- as.character(aff[, 2])
  
  vcf_vcf1 <- NULL
  if (multiSample) {
    vcf_vcf1 <- readVcf(unlist(input_list), genomeVersion)
    
    carriers <- switch (nrow(unaff),
                        one_ind_ms(vcf_vcf1, "0/1", unaff, filterTag),
                        two_ind_ms(vcf_vcf1, "0/1", unaff, filterTag),
                        three_ind_ms(vcf_vcf1, "0/1", unaff, filterTag),
                        four_ind_ms(vcf_vcf1, "0/1", unaff, filterTag),
                        five_ind_ms(vcf_vcf1, "0/1", unaff, filterTag))
    
    affected <- switch (nrow(aff)+1,
                        stop("No affected individuals detected. Something might be wrong with the .ped file..."),
                        one_ind_ms(vcf_vcf1, "1/1", aff, filterTag),
                        two_ind_ms(vcf_vcf1, "1/1", aff, filterTag),
                        three_ind_ms(vcf_vcf1, "1/1", aff, filterTag),
                        four_ind_ms(vcf_vcf1, "1/1", aff, filterTag),
                        five_ind_ms(vcf_vcf1, "1/1", aff, filterTag))
  } else {
   
    ## by now we have to disable this
    stop("More than one unique sample VCF files is currently unsupported. Please do the variant call with all the individuals altogether or proceed analysing each file separately")

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
    
    carriers <- switch(length(input_list_unaff),
                       one_ind_us(input_list_unaff, "0/1", filterTag, genomeVersion),
                       two_ind_us(input_list_unaff, "0/1", filterTag, genomeVersion),
                       three_ind_us(input_list_unaff, "0/1", filterTag, genomeVersion),
                       four_ind_us(input_list_unaff, "0/1", filterTag, genomeVersion),
                       five_ind_us(input_list_unaff, "0/1", filterTag, genomeVersion))
    
    affected <- switch(length(input_list_aff)+1,
                       stop("No affected individuals detected. Something might be wrong with the .ped file..."),
                       one_ind_us(input_list_aff, "1/1", filterTag, genomeVersion),
                       two_ind_us(input_list_aff, "1/1", filterTag, genomeVersion),
                       three_ind_us(input_list_aff, "1/1", filterTag, genomeVersion),
                       four_ind_us(input_list_aff, "1/1", filterTag, genomeVersion),
                       five_ind_us(input_list_aff, "1/1", filterTag, genomeVersion))
    
  }
  
  ### shared variants
  
  ### now lets select only those which fulfill the inheritance model
  
  recessive <- affected
  if (length(carriers) >= 1) {
    realcommonrecessive <- sharedVariants(affected, carriers)
    recessive <- affected[realcommonrecessive]
  }
  
  recessive <- matchChromosomeNames(recessive, txdb)

  ##########################
  ##                      ##
  ##      ANNOTATION      ##
  ##                      ##
  ##########################
  
  recessive_annotated <- annotationEngine(recessive, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
                                          radicalAAchangeMatrix=radicalAAchangeMatrix,
                                          otherAnnotations=otherAnnotations,
                                          allTranscripts=allTranscripts, BPPARAM=BPPARAM)

  ##########################
  ##                      ##
  ## BUILD RESULTS OBJECT ##
  ##                      ##
  ##########################

  cnAF <- colnames(mcols(recessive_annotated))
  cnAF <- cnAF[grep("AF", cnAF)]
  MAFpopMask <- rep(TRUE, length(cnAF))
  names(MAFpopMask) <- cnAF

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="autosomal recessive homozygous", variants=recessive_annotated,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any",
      aaChangeType="Any", MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
