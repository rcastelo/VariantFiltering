
## function to score DNA-binding sites including variants. it produces scores for the binding site with
## the reference and alternative alleles of SNVs. it assumes that the input variantsVR is a VRanges object
.scoreBindingSiteVariants <- function(variantsVR, weightMatrices, bsgenome, BPPARAM=bpparam("SerialParam")) {

  ## adapt to sequence style and genome version from the input
  ## BSgenome object, thus assuming positions are based on the same
  ## genome even though might be differently specified (i.e., hg38 vs GRCh38)
  seqlevelsStyle(variantsVR) <- seqlevelsStyle(bsgenome)
  commonChr <- intersect(seqlevels(variantsVR), seqlevels(bsgenome))
  if (any(is.na(genome(variantsVR)))) {
    warning(sprintf("Assuming the genome build of the input variants is %s.", unique(genome(bsgenome)[commonChr])))
    genome(variantsVR) <- genome(bsgenome)
  } else if (any(genome(variantsVR)[commonChr] != genome(bsgenome)[commonChr])) {
    warning(sprintf("Assumming %s represent the same genome build.",
                    paste(c(unique(genome(variantsVR)[commonChr]), unique(genome(bsgenome)[commonChr])),
                          collapse=" and ")))
    genome(variantsVR) <- genome(bsgenome)
  }

  wmNames <- sapply(weightMatrices, wmName)
  mask <- nchar(wmNames) == 0 || is.na(wmNames)
  wmNames[mask] <- paste("WM", 1:sum(mask))

  ## create a DataFrame object for binding site score annotations
  dummyDF <- rep(list(rep(NA_real_, length(variantsVR))), length(wmNames)*3)
  names(dummyDF) <- paste(c("SCOREwmREF", "SCOREwmALT", "SCOREwmPOS"), ## "SITEwmREF", "CHRwmREF", "POSwmREF"),
                          rep(wmNames, each=3), sep="_")
  dummyDF <- do.call("DataFrame", dummyDF)
  ## dummyDF[[paste0("SITEwmREF_", wmNames)]] <- NA_character_
  ## dummyDF[[paste0("CHRwmREF_", wmNames)]] <- NA_character_
  ## dummyDF[[paste0("POSwmREF_", wmNames)]] <- NA_integer_

  snvMask <- rep(TRUE, length(variantsVR))
  if (!is.null(variantsVR$TYPE)) ## if variant type is annotated then use it to restrict to SNVs
    snvMask <- variantsVR$TYPE == "SNV"

  if (all(!snvMask)) ## if there are no SNVs, then there is nothing to score (by now)
    return(dummyDF)

  altAlleles <- alt(variantsVR[snvMask])
  strandAlleles <- rep("*", sum(snvMask))
  if (!is.null(variantsVR$LOCSTRAND[snvMask])) { ## if variant is annotated with strand
    altAlleles <- VariantFiltering:::.adjustForStrandSense(variantsVR[snvMask], altAlleles)
    strandAlleles <- variantsVR$LOCSTRAND[snvMask]
  }

  if (is.character(altAlleles))
    altAlleles <- DNAStringSetList(strsplit(as.character(altAlleles), split="", fixed=TRUE))
  stopifnot(class(altAlleles) == "DNAStringSetList") ## QC

  altAllelesReverse <- relist(reverseComplement(unlist(altAlleles)), altAlleles)

  for (i in seq(along=weightMatrices)) {
    message(sprintf("Scoring binding sites with weight matrix %s.", wmNames[i]))

    wmObj <- weightMatrices[[i]]
    wmColNames <- paste(c("SCOREwmREF", "SCOREwmALT", "SCOREwmPOS"), wmNames[i], sep="_")
    ## wmColNames <- paste(c("SCOREwmREF", "SCOREwmALT", "SCOREwmPOS", "SITEwmREF", "CHRwmREF", "POSwmREF"), wmNames[i], sep="_")

    locMask <- rep(TRUE, length(variantsVR))
    if (!is.null(variantsVR$LOCATION))
      locMask <- variantsVR$LOCATION %in% wmLocations(wmObj)

    if (any(snvMask & locMask))
      dummyDF[snvMask & locMask, wmColNames] <-
        .scoreBindingSiteVariantsThisWM(variantsVR[snvMask & locMask],
                                        altAlleles[locMask[snvMask]],
                                        altAllelesReverse[locMask[snvMask]],
                                        strandAlleles[locMask[snvMask]],
                                        wmObj, bsgenome, BPPARAM)
  }

  mcols(dummyDF) <- DataFrame(TAB=rep("BindingSites", ncol(dummyDF)))

  dummyDF
}

.scoreBindingSiteVariantsThisWM <- function(variantsVR, altAlleles, altAllelesReverse,
                                            strandAlleles, wmObj, bsgenome,
                                            BPPARAM=bpparam("SerialParam")) {

  widthWm <- width(wmObj)
  wmName <- wmName(wmObj)
 
  ## build GRanges to fetch genome sequence through getSeq()
  bsite_window <- GRanges(seqnames=seqnames(variantsVR),
                  ranges=ranges(variantsVR),
                  strand=strandAlleles)

  ## retrieve regions around the allele potentially involving binding sites
  wregion <- width(wmObj)*2-1
  bsite_window <- resize(bsite_window, width=wregion, fix="center")

  ## fetch sequences and generate corresponding ones replacing the alternative alleles
  bsite_window_strings <- getSeq(bsgenome, bsite_window)
  bsite_window_ALT_strings <- replaceAt(bsite_window_strings,
                                        IRanges(start=widthWm, end=widthWm),
                                        altAlleles)

  ## fetch sequences from the reverse strand whenever there is no strand annotated
  maskRev <- rep(TRUE, length(variantsVR))

  if (!is.null(variantsVR$LOCSTRAND) && !is.null(variantsVR$LOCATION))
    maskRev <- as.vector(variantsVR$LOCSTRAND == "*" |
                         variantsVR$LOCATION == "promoter")

  bsite_window_strings_rev <- bsite_window_ALT_strings_rev <- DNAStringSet()
  if (any(maskRev)) {
    bsite_window_strings_rev <- reverseComplement(bsite_window_strings[maskRev])
    bsite_window_ALT_strings_rev <- replaceAt(bsite_window_strings_rev,
                                              IRanges(start=widthWm, end=widthWm),
                                              altAllelesReverse[maskRev])
  }

  ## deal with binding sites that can only be scored strictly within
  ## the boundaries of the location of the weight matrix (e.g., fiveSpliceSite and threeSpliceSite)
  lefttrimmednt <- righttrimmednt <- rep(0, length(variantsVR))
  strictLocMask <- wmStrictLocations(wmObj)
  if (any(strictLocMask) && !is.null(variantsVR$LOCSTART) && !is.null(variantsVR$LOCEND)) {
    trimMask <- rep(TRUE, length(variantsVR))
    if (length(strictLocMask) > 1) { ## eventually not all locations should be strict
      names(strictLocMask) <- wmLocations(wmObj)
      stopifnot(all(variantsVR$LOCATION %in% names(strictLocMask))) ## QC
      trimMask <- strictLocMask[as.character(variantsVR$LOCATION)]
    }
    locgr <- GRanges(seqnames=seqnames(variantsVR),
                     ranges=IRanges(variantsVR$LOCSTART, variantsVR$LOCEND),
                     strand=variantsVR$LOCSTRAND)
    locgr <- pintersect(bsite_window, locgr)
    
    ## CURRENTLY LOCSTART AND LOCEND CONTAIN COORDINATES RELATIVE TO EACH ANNOTATED FEATURE
    ## EXCEPT IN fiveSpliceSite and threeSpliceSite THAT CONTAIN THE FEATURE COORDINATES
    ## UNTIL LOCSTART AND LOCEND CONTAIN START AND END OF ALL ANNTOTATED FEATURES LET'S JUST
    ## PUT BACK THE REGION COORDINATES WHEREVER THE INTERSECTION DOES NOT FIT THE WEIGHT MATRIX
    ## THIS WILL MAKE WORK STRICT LOCATIONS ONLY IN fiveSpliceSite AND threeSpliceSite BY NOW
    ranges(locgr[width(locgr) < widthWm]) <- ranges(bsite_window[width(locgr) < widthWm])

    lefttrimmednt[trimMask] <- start(locgr)[trimMask] - start(bsite_window)[trimMask]
    righttrimmednt[trimMask] <- end(bsite_window)[trimMask] - end(locgr)[trimMask]
    bsite_window_strings <- subseq(bsite_window_strings, start=lefttrimmednt+1,
                                   end=width(bsite_window_strings)-righttrimmednt)
    bsite_window_ALT_strings <- subseq(bsite_window_ALT_strings, start=lefttrimmednt+1,
                                       end=width(bsite_window_ALT_strings)-righttrimmednt)
    if (any(maskRev)) {
      bsite_window_strings_rev <- subseq(bsite_window_strings_rev, start=lefttrimmednt+1,
                                         end=width(bsite_window_strings_rev)-righttrimmednt)
      bsite_window_ALT_strings_rev <- subseq(bsite_window_ALT_strings_rev, start=lefttrimmednt+1,
                                             end=width(bsite_window_ALT_strings_rev)-righttrimmednt)
    }
  }
  totaltrimmednt <- lefttrimmednt + righttrimmednt
  squaredregions <- all(totaltrimmednt[1] == totaltrimmednt)

  ## score potential binding sites involving ALT alleles
  bsite_window_ALT_scores <- numeric(0)
  if (length(bsite_window_ALT_strings) > 0) { ## bpvec() returns an empty list() when input is empty
    if (BPPARAM$workers > 1)
      bsite_window_ALT_scores <- bpvec(X=bsite_window_ALT_strings,
                                       FUN=wmScore, object=wmFilename(wmObj), BPPARAM=BPPARAM)
    else
      bsite_window_ALT_scores <- wmScore(object=wmObj, dnaseqs=bsite_window_ALT_strings)
    bsite_window_ALT_scores <- unlist(bsite_window_ALT_scores, use.names=FALSE)
  }
  if (squaredregions)
    bsite_window_ALT_scores <- matrix(bsite_window_ALT_scores, ncol=widthWm-totaltrimmednt[1], byrow=TRUE)
  else
    bsite_window_ALT_scores <- split(bsite_window_ALT_scores, rep(1:length(variantsVR), times=widthWm-totaltrimmednt))

  ## score potential binding sites involving ALT alleles in those
  ## reverse complemented regions annotated without strand
  if (length(bsite_window_ALT_strings_rev) > 0) { ## bpvec() returns an empty list() when input is empty
    if (BPPARAM$workers > 1)
      sco_rev <- bpvec(X=bsite_window_ALT_strings_rev, FUN=wmScore, object=wmFilename(wmObj), BPPARAM=BPPARAM)
    else
      sco_rev <- wmScore(object=wmObj, dnaseqs=bsite_window_ALT_strings_rev)
    sco_rev <- unlist(sco_rev, use.names=FALSE)
    bsite_window_ALT_scores_rev <- vector(mode="list", length=length(variantsVR))
    if (squaredregions) {
      bsite_window_ALT_scores_rev <- matrix(NA, nrow=length(variantsVR),
                                            ncol=widthWm-totaltrimmednt[1])
      bsite_window_ALT_scores_rev[maskRev, ] <- matrix(sco_rev, ncol=widthWm-totaltrimmednt, byrow=TRUE)
      bsite_window_ALT_scores <- cbind(bsite_window_ALT_scores, bsite_window_ALT_scores_rev)
    } else {
      bsite_window_ALT_scores_rev[maskRev] <- split(sco_rev, rep(1:sum(maskRev), times=widthWm-totaltrimmednt[maskRev]))
      bsite_window_ALT_scores <- mapply(c, bsite_window_ALT_scores, bsite_window_ALT_scores_rev)
    }
    rm(sco_rev)
    rm(bsite_window_ALT_scores_rev)
  }

  ## fetch the maximum score and the binding site position leading to this maximum score
  bsite_max_ALT_scores <- matrix(NA_real_, nrow=length(variantsVR), ncol=2)
  if (squaredregions) {
    bsite_max_ALT_scores <- cbind(bsite_window_ALT_scores[, 1], 1)
    if (ncol(bsite_window_ALT_scores) > 1) {
      bsite_max_ALT_scores <- t(apply(bsite_window_ALT_scores, 1, function(x) {
                                        maxsco <- maxpos <- NA_real_
                                        if (any(!is.na(x))) {
                                          maxpos <- which.max(x)
                                          maxsco <- x[maxpos]
                                        }
                                        c(maxsco, maxpos)
                                      }))
    }
  } else {
    bsite_max_ALT_scores <- t(sapply(bsite_window_ALT_scores, function(x) {
                                       maxsco <- maxpos <- NA_real_
                                       if (any(!is.na(x))) {
                                         maxpos <- which.max(x)
                                         maxsco <- x[maxpos]
                                       } else if (length(x) == 1) ## if there's no window ensure ALT & REF are scored
                                         maxpos <- 1              ## thus assuming no window region => no search
                                       c(maxsco, maxpos)
                                     }))
  }

  ## score REF alleles at the binding site of the highest score with the ALT allele
  posHighestSitesALT <- bsite_max_ALT_scores[, 2]
  bsite_strings_at_ALT <- DNAStringSet()
  if (length(posHighestSitesALT) > 0 && any(!is.na(posHighestSitesALT))) {
    mask <- !is.na(posHighestSitesALT)
    bsite_strings_at_ALT <- bsite_window_strings[mask]
    bsite_strings_at_ALT[posHighestSitesALT[mask] > widthWm-totaltrimmednt[mask]] <-
      bsite_window_strings_rev[posHighestSitesALT[mask & maskRev] > widthWm-totaltrimmednt[mask & maskRev]]
    maskFwd <- posHighestSitesALT[mask] <= widthWm-totaltrimmednt[mask]
    bsite_strings_at_ALT[maskFwd] <- subseq(bsite_strings_at_ALT[maskFwd],
                                            start=posHighestSitesALT[mask][maskFwd], width=widthWm)
    bsite_strings_at_ALT[!maskFwd] <- subseq(bsite_strings_at_ALT[!maskFwd],
                                             start=posHighestSitesALT[mask][!maskFwd]-widthWm+totaltrimmednt[mask][!maskFwd],
                                             width=widthWm)
  }

  bsite_REF_scores <- numeric(0)
  if (length(bsite_strings_at_ALT) > 0) { ## bpvec() returns an empty list() when input is empty
    if (BPPARAM$workers > 1)
      bsite_REF_scores <- bpvec(X=bsite_strings_at_ALT, FUN=wmScore, object=wmFilename(wmObj), BPPARAM=BPPARAM)
    else
      bsite_REF_scores <- wmScore(object=wmObj, dnaseqs=bsite_strings_at_ALT)
    bsite_REF_scores <- unlist(bsite_REF_scores, use.names=FALSE)
  }

  ## store the position of the allele respect to the position of the first conserved nucleotide
  ## if there is no such nucleotide, then the position is with respect to the first nucleotide of
  ## the binding site
  wmRefPos <- conservedPositions(wmObj)
  if (length(wmRefPos) == 0)
    wmRefPos <- 1
  wmRefPos <- wmRefPos[1]
  relposbsite <- integer(0)
  if (wmRefPos > 1)
    relposbsite <- seq(-wmRefPos+1, -1, by=1)
  relposbsite <- c(relposbsite, seq(1, widthWm-length(relposbsite)))

  bsite_scores_pos <- rep(list(rep(NA_real_, length(variantsVR))), 3)
  names(bsite_scores_pos) <- paste(c("SCOREwmREF", "SCOREwmALT", "SCOREwmPOS"), wmName, sep="_")
  bsite_scores_pos <- do.call("DataFrame", bsite_scores_pos)

  bsite_scores_pos[[paste("SCOREwmALT", wmName, sep="_")]] <- round(bsite_max_ALT_scores[, 1], digits=2)

  ## positions corresponding to hits in the negative strand are
  ## coded as negative integers. when the weight matrix has at least
  ## one fully conserved position, the sign is reversed
  bsite_max_ALT_scores[, 2] <- ifelse(!is.na(bsite_max_ALT_scores[, 2]) & (bsite_max_ALT_scores[, 2] > widthWm-totaltrimmednt),
                                      -1*(start(variantsVR) - (start(bsite_window)+lefttrimmednt+bsite_max_ALT_scores[, 2]-widthWm-1) + 1),
                                      start(variantsVR) - (start(bsite_window)+lefttrimmednt+bsite_max_ALT_scores[, 2]-1) + 1)
  bsite_scores_pos[[paste("SCOREwmPOS", wmName, sep="_")]] <-
    sign(bsite_max_ALT_scores[, 2])*relposbsite[abs(bsite_max_ALT_scores[, 2])]
  mask <- !is.na(bsite_scores_pos[[paste("SCOREwmPOS", wmName, sep="_")]])
  bsite_scores_pos[[paste("SCOREwmREF", wmName, sep="_")]][mask] <- round(bsite_REF_scores, digits=2)

  ## bsite_scores_pos[[paste("SITEwmREF", wmName, sep="_")]] <- as.character(bsite_strings_at_ALT)
  ## bsite_scores_pos[[paste("CHRwmREF", wmName, sep="_")]] <- as.character(seqnames(variantsVR))
  ## bsite_scores_pos[[paste("POSwmREF", wmName, sep="_")]] <- start(variantsVR)

  bsite_scores_pos
}
