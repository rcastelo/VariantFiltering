setMethod("deNovo", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  callobj <- match.call()
  callstr <- deparse(callobj)
  input_list <- as.list(path(param$vcfFiles))
  ped <- param$pedFilename
  sinfo <- param$seqInfos[[1]]
  orgdb <- param$orgdb
  txdb <- param$txdb
  snpdb <- param$snpdb
  radicalAAchangeMatrix <- param$radicalAAchangeMatrix
  allTranscripts <- param$allTranscripts
  otherAnnotations <- param$otherAnnotations
  filterTag <- param$filterTag
 
  genomeInfo <- sinfo
  
  if (!exists(as.character(substitute(BPPARAM))))
    stop(sprintf("Parallel back-end function %s given in argument 'BPPARAM' does not exist in the current workspace. Either you did not write correctly the function name or you did not load the package 'BiocParallel'.", as.character(substitute(BPPARAM))))
  
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
  unaff_dad <- unaff[unaff[, 5] == 1, ][, 2]
  unaff_mom <- unaff[unaff[, 5] == 2, ][, 2]
  
  aff <- pedf[pedf[, 6] == 2, ]
  aff_ind <- as.character(aff[, 2])
  
  if (!as.character(aff[1, 3]) == as.character(unaff_dad) && as.character(aff[1, 4]) == as.character(unaff_mom)) {
    stop("This analysis only admits both unaffected parents and an affected child. Please check the ped file provided")
  }
  
  if (multiSample) {
    message("Reading input VCF file into main memory.")
    vcf_vcf1 <- readVcf(unlist(input_list), genomeInfo)
    
    unaffected <- two_ind_ms(vcf_vcf1, "0/0", unaff, filterTag)
    affected <- vcf2GR_2options(vcf_vcf1, "0/1", "1/1", aff[, 2], filterTag)
    
    realdenovo <- sharedVariants(affected, unaffected)
    childDeNovo <- affected[realdenovo]
    
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
    
    message("Reading input VCF files into main memory.")
    vcf_carrier1 <- readVcf(unlist(input_list_father), filterTag, genomeInfo)
    vcf_carrier2 <- readVcf(unlist(input_list_mother), filterTag, genomeInfo)
    vcf_affected1 <- readVcf(unlist(input_list_aff), filterTag, genomeInfo)
    
   
    # we are interested in comparing all changes too see which of them are inherited, so we don't filter by any genotypic status 
    
    gr_carr1 <- vcf2GR_2options(vcf_carrier1, "0/1", "1/1", 1)
    gr_carr2 <- vcf2GR_2options(vcf_carrier2, "0/1", "1/1", 1)
    gr_aff1 <- vcf2GR_2options(vcf_affected1, "0/1", "1/1", 1)
            
    ### discard shared variants
   
    realcommongen11 <- sharedVariants(gr_aff1, gr_carr1)
    childNotParent1 <- gr_aff1[-realcommongen11]
    
    realcommongen12 <- sharedVariants(childNotParent1, gr_carr2)
    childDeNovo <- childNotParent1[-realcommongen12]
    
  } 
  
  denovo <- matchChromosomes(childDeNovo, txdb)
    
  ##########################
  ##                      ##
  ##      ANNOTATION      ##
  ##                      ##
  ##########################
  
  denovo_annotated <- annotationEngine(denovo, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
                                       radicalAAchangeMatrix=radicalAAchangeMatrix,
                                       otherAnnotations=otherAnnotations,
                                       allTranscripts=allTranscripts, BPPARAM=BPPARAM)

  ##########################
  ##                      ##
  ## BUILD RESULTS OBJECT ##
  ##                      ##
  ##########################

  cnAF <- colnames(mcols(denovo_annotated))
  cnAF <- cnAF[grep("AF", cnAF)]
  MAFpopMask <- rep(TRUE, length(cnAF))
  names(MAFpopMask) <- cnAF

  new("VariantFilteringResults", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="de novo", variants=denovo_annotated,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any", aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_)
})
