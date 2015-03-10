setMethod("allInheritanceModels", signature(param="VariantFilteringParam"),
          function(param, BPPARAM=bpparam()) {

  stop("this is being updated and cannot be used at the moment.")

  callobj <- match.call()
  callstr <- deparse(callobj)
  input_list <- as.list(path(param$vcfFiles))
##  ped <- param@pedFilename
  sinfo <- param$seqInfos[[1]]
  orgdb <- param$orgdb
  txdb <- param$txdb
  snpdb <- param$snpdb
  bsgenome <- param$bsgenome
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
    stop("Only one single VCF file is supported.")
  } else if (length(input_list) == 1) {
    multiSample <- TRUE
  } else {
    stop("A VCF file has to be provided")
  }
  
  ## pedf <- read.table(ped, header=F)
  
  ## unaff <- pedf[pedf[, 6] == 1, ]
  ## unaff_ind <- as.character(unaff[, 2])

  ## aff <- pedf[pedf[, 6] == 2, ]
  ## aff_ind <- as.character(aff[, 2])

  ## unaff_males <- pedf[pedf[, 6] == 1 & pedf[, 5] == 1, ]
  ## unaff_males_ind <- as.character(unaff_males[, 2])
  
  ## carrier_female <-  pedf[pedf[, 6] == 1 & pedf[, 5] == 2, ]
  ## carrier_female_ind <- as.character(carrier_female[, 2])

  ## unaff_dad <- unaff[unaff[, 5] == 1, ][, 2]
  ## unaff_mom <- unaff[unaff[, 5] == 2, ][, 2]
  
  message("Reading input VCF file into main memory.")
  vcf1 <- readVcf(unlist(input_list), genomeInfo)

  gr1 <- rowRanges(vcf1)
  gr1 <- .matchSeqinfo(gr1, txdb, bsgenome)

  ##########################
  ##                      ##
  ##      ANNOTATION      ##
  ##                      ##
  ##########################

  allinheritance_annotated <-  annotationEngine(gr1, orgdb=orgdb, txdb=txdb, snpdb=snpdb,
                                          radicalAAchangeMatrix=radicalAAchangeMatrix,
                                          otherAnnotations=otherAnnotations,
                                          allTranscripts=allTranscripts, BPPARAM=BPPARAM)

  
    for (i in seq_along(colnames(geno(vcf1)$GT))) {
      ind <- colnames(geno(vcf1)$GT)[i]
      homalt <- geno(vcf1)$GT[, ind] == "1/1" |
                geno(vcf1)$GT[, ind] == "1|1"
      homaltind <- DataFrame(homalt[names(allinheritance_annotated)])
      colnames(homaltind) <- paste0("homALT_", ind)

      mcols(allinheritance_annotated) <- cbind(mcols(allinheritance_annotated), homaltind)
  }

    for (i in seq_along(colnames(geno(vcf1)$GT))) {
      ind <- colnames(geno(vcf1)$GT)[i]
      het <- geno(vcf1)$GT[, ind] == "0/1" |
             geno(vcf1)$GT[, ind] == "0|1" |
             geno(vcf1)$GT[, ind] == "1/0" |
             geno(vcf1)$GT[, ind] == "1|0"
      hetind <- DataFrame(het[names(allinheritance_annotated)])
      colnames(hetind) <- paste0("het_", ind)

      mcols(allinheritance_annotated) <- cbind(mcols(allinheritance_annotated), hetind)
  }

    for (i in seq_along(colnames(geno(vcf1)$GT))) {
      ind <- colnames(geno(vcf1)$GT)[i]
      homref <- geno(vcf1)$GT[, ind] == "0/0" |
                geno(vcf1)$GT[, ind] == "0|0"
      homrefind <- DataFrame(homref[names(allinheritance_annotated)])
      colnames(homrefind) <- paste0("homREF_", ind)

      mcols(allinheritance_annotated) <- cbind(mcols(allinheritance_annotated), homrefind)
  }


    homaltcols <- grep("homALT_", colnames(mcols(allinheritance_annotated)))
    homaltcolsvars <- allinheritance_annotated[, homaltcols]
    df_homaltcolsvars <- as.data.frame(mcols(homaltcolsvars))
    indvarshomalt <- lapply(seq_along(homaltcolsvars), function(x) gsub("homALT_", "", colnames(df_homaltcolsvars)[which(  df_homaltcolsvars[x, ] == TRUE)]))

    DFindhomalt <- DataFrame(HomozygousAlt=CharacterList(indvarshomalt))
    mcols(allinheritance_annotated) <- cbind(mcols(allinheritance_annotated), DFindhomalt)


    hetcols <- grep("het", colnames(mcols(allinheritance_annotated)))
    hetcolsvars <- allinheritance_annotated[, hetcols]
    df_hetcolsvars <- as.data.frame(mcols(hetcolsvars))
    indvarshet <- lapply(seq_along(hetcolsvars), function(x) gsub("het_", "", colnames(df_hetcolsvars)[which(df_hetcolsvars[x, ] == TRUE)]))

    DFindhet <- DataFrame(Heterozygous=CharacterList(indvarshet))
    mcols(allinheritance_annotated) <- cbind(mcols(allinheritance_annotated), DFindhet)
    
    
    homrefcols <- grep("homREF_", colnames(mcols(allinheritance_annotated)))
    homrefcolsvars <- allinheritance_annotated[, homrefcols]
    df_homrefcolsvars <- as.data.frame(mcols(homrefcolsvars))
    indvarshomref <- lapply(seq_along(homrefcolsvars), function(x) gsub("homREF_", "", colnames(df_homrefcolsvars)[which(  df_homrefcolsvars[x, ] == TRUE)]))

    DFindhomref <- DataFrame(HomozygousRef=CharacterList(indvarshomref))
    mcols(allinheritance_annotated) <- cbind(mcols(allinheritance_annotated), DFindhomref)


  ##########################
  ##                      ##
  ## BUILD RESULTS OBJECT ##
  ##                      ##
  ##########################

  locMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$LOCATION)),
                                     levels(annotated_variants$LOCATION)))
  conMask <- do.call("names<-", list(rep(TRUE, nlevels(annotated_variants$CONSEQUENCE)),
                                     levels(annotated_variants$CONSEQUENCE)))
  MAFpopMask <- NA
  if ("MafDb" %in% sapply(param$otherAnnotations, class)) {
    ## assume AF columns are those containing AF[A-Z]+ and being of class 'numeric'
    cnAF <- colnames(mcols(annotated_variants))
    colsclasses <- sapply(mcols(annotated_variants), class)
    cnAF <- cnAF[intersect(grep("AF[A-Z]+", cnAF), grep("numeric", colsclasses))]
    MAFpopMask <- rep(TRUE, length(cnAF))
    names(MAFpopMask) <- cnAF
  }

    
  new("VariantFilteringResultsAIM", callObj=callobj, callStr=callstr, inputParameters=param,
      inheritanceModel="any", variants=allinheritance_annotated, bamViews=BamViews(),
      indselected=NA_character_, inheritancepattern="None", selectindexcase=NA_character_,
      selectcarrierrelative1=NA_character_, selectcarrierrelative2=NA_character_, selectaffectedrelative=NA_character_,
      selectcarrierallele1ch=NA_character_, selectcarrierallele2ch=NA_character_, selectaffrelative1ch=NA_character_,
      selectunaffectedrelative1ad=NA_character_, selectunaffectedrelative2ad=NA_character_, selectaffectedrelative1ad=NA_character_,
      selectcarrierrelativefemale1xl=NA_character_, selectaffectedrelativemale1xl=NA_character_, selectunaffectedmale1xl=NA_character_,
      selectparent1dn=NA_character_, selectparent2dn=NA_character_,
      dbSNPflag=NA_character_, OMIMflag=NA_character_, variantType="Any",
      locationMask=locMask, consequenceMask=conMask, aaChangeType="Any",
      MAFpopMask=MAFpopMask, naMAF=TRUE, maxMAF=1,
      minPhastCons=NA_real_, minPhylostratumIndex=NA_integer_,
      minCRYP5ss=NA_real_, minCRYP3ss=NA_real_, minCUFC=0)
})
