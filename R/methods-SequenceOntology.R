## form URLs as
## http://sequenceontology.org/browser/current_release/term/SO:0001578
## 
## SO:0001619 nc_transcript_variant (in a non-coding RNA gene)
## SO:0001578 stop_lost (elongated transcript)

## function to annotate a sequence ontology graph given as
## a 'graph' object in the second argument 'gSO' with
## variants given in the 'VRanges' first argument object 'variantsVR'
annotateSO <- function(variantsVR, gSO) {

  if (length(variantsVR) == 0)
    return(gSO)

  gSOattr <- names(nodeData(gSO, n=nodes(gSO)[1])[[1]])

  if (!"vcfIdx" %in% gSOattr) ## if the attribute 'vcfIdx' does not exist, create it
    nodeDataDefaults(gSO, "vcfIdx") <- integer(0)

  if (!"varIdx" %in% gSOattr) ## if the attribute 'varIdx' does not exist, create it
    nodeDataDefaults(gSO, "varIdx") <- integer(0)

  VRangesIdx <- 1:length(variantsVR)

  annotatedVCFIDX <- annotatedVARIDX <- integer(0)

  ## SO:0001819 synonymous_variant
  mask <- variantsVR$CONSEQUENCE %in% "synonymous"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001819", "vcfIdx") <- list("SO:0001819"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001819", "varIdx") <- list("SO:0001819"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001627 intron_variant
  mask <- variantsVR$LOCATION %in% "intron"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001627", "vcfIdx") <- list("SO:0001627"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001627", "varIdx") <- list("SO:0001627"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001587 stop_gained (premature stop codon)
  mask <- variantsVR$CONSEQUENCE %in% "nonsense"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001587", "vcfIdx") <- list("SO:0001587"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001587", "varIdx") <- list("SO:0001587"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0002012 start_lost
  mask <- start(variantsVR$CDSLOC) == 1 & as.character(variantsVR$REFCODON) %in% "ATG" &
          !as.character(variantsVR$VARCODON) %in% "ATG"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0002012", "vcfIdx") <- list("SO:0002012"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0002012", "varIdx") <- list("SO:0002012"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }
 
  ## SO:0001589 frameshift_variant
  mask <- variantsVR$CONSEQUENCE %in% "frameshift"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001589", "vcfIdx") <- list("SO:0001589"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001589", "varIdx") <- list("SO:0001589"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001583 missense_variant (amino acid change)
  mask <- variantsVR$CONSEQUENCE %in% "nonsynonymous"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001583", "vcfIdx") <- list("SO:0001583"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001583", "varIdx") <- list("SO:0001583"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001624 3_prime_UTR_variant
  mask <- variantsVR$LOCATION %in% "threeUTR"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001624", "vcfIdx") <- list("SO:0001624"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001624", "varIdx") <- list("SO:0001624"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }
  
  ## SO:0001623 5_prime_UTR_variant
  mask <- variantsVR$LOCATION %in% "fiveUTR"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001623", "vcfIdx") <- list("SO:0001623"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001623", "varIdx") <- list("SO:0001623"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001575 splice_donor_variant
  mask <- variantsVR$LOCATION %in% "fiveSpliceSite"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001575", "vcfIdx") <- list("SO:0001575"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001575", "varIdx") <- list("SO:0001575"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001574 splice_acceptor_variant
  mask <- variantsVR$LOCATION %in% "threeSpliceSite"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001574", "vcfIdx") <- list("SO:0001574"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001574", "varIdx") <- list("SO:0001574"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001628 intergenic
  mask <- variantsVR$LOCATION %in% "intergenic"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001628", "vcfIdx") <- list("SO:0001628"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001628", "varIdx") <- list("SO:0001628"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001820 inframe_indel
  mask <- variantsVR$TYPE %in% c("Insertion", "Deletion") & variantsVR$LOCATION %in% "coding" &
          width(variantsVR) == 3
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001820", "vcfIdx") <- list("SO:0001820"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001820", "varIdx") <- list("SO:0001820"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## SO:0001631 upstream_gene_variant (intergenic)
  mask <- variantsVR$LOCATION %in% "promoter"
  if (any(mask)) {
    vcfidx <- unique(variantsVR$VCFIDX[mask])
    nodeData(gSO, "SO:0001631", "vcfIdx") <- list("SO:0001631"=vcfidx)
    annotatedVCFIDX <- c(annotatedVCFIDX, vcfidx)
    varidx <- VRangesIdx[mask]
    nodeData(gSO, "SO:0001631", "varIdx") <- list("SO:0001631"=varidx)
    annotatedVARIDX <- c(annotatedVARIDX, varidx)
  }

  ## S0:0001060 sequence_variants (everything else that is not annotated to other SO term)
  missingVCFIDX <- setdiff(unique(variantsVR$VCFIDX), unique(annotatedVCFIDX))
  missingVARIDX <- setdiff(VRangesIdx, unique(annotatedVARIDX))
  if (length(missingVCFIDX) > 0)
    nodeData(gSO, "SO:0001060", "vcfIdx") <- list("SO:0001060"=missingVCFIDX)
  if (length(missingVARIDX) > 0)
    nodeData(gSO, "SO:0001060", "varIdx") <- list("SO:0001060"=missingVARIDX)

  gSO
}

## add SO filters and cutoffs metadata to the input variants VRanges object
addSOmetadata <- function(variantsVR) {
  sofilter <- function(x) {
                soterms <- VariantFiltering::cutoffs(x)$SOterms
                if (is.null(soterms) || all(is.na(soterms)) ||
                    length(soterms) == 0 || soterms[1] == "Any") ## no SO terms specified imply no restriction
                  soterms <- graph::nodes(VariantFiltering::sog(x))[sapply(graph::nodeData(VariantFiltering::sog(x), graph::nodes(VariantFiltering::sog(x)), "varIdx"), length) > 0]
                else {
                  soterms <- VariantFiltering:::.findSOIDs(VariantFiltering::cutoffs(x)$SOterms, VariantFiltering::sog(x))
                  ## restrict to the given SO terms and all their more
                  ## specific ancestors with at least one annotated variant
                  soterms <- unique(c(soterms, VariantFiltering:::.ancestorsSO(VariantFiltering::param(x), soterms)))
                  soterms <- soterms[sapply(graph::nodeData(VariantFiltering::sog(x), soterms, "varIdx"), length) > 0]
                }

                mask <- rep(TRUE, length(x))
                if (length(soterms) > 0) {
                  varidxinsog <- unique(unlist(graph::nodeData(VariantFiltering::sog(x), soterms, "varIdx"), use.names=FALSE))
                  mask <- rep(FALSE, length(x))
                  mask[varidxinsog] <- TRUE
                }

                mask
              }
  attr(sofilter, "description") <- "Sequence Ontology annotations"
  attr(sofilter, "TAB") <- "Transcript"
  environment(sofilter) <- baseenv()
  SOmetadata <- list(filters=list(SOterms=sofilter),
                     cutoffs=CutoffsList(SOterms="Any"))
  metadata(mcols(variantsVR))$filters <- c(metadata(mcols(variantsVR))$filters, SOmetadata$filters)
  metadata(mcols(variantsVR))$cutoffs <- c(metadata(mcols(variantsVR))$cutoffs, SOmetadata$cutoffs)

  variantsVR
}


## private functions

## builds a descendants matrix of an input acyclic digraph 'g'
## the resulting matrix has for every vertex in the rows a logical
## mask indicating what other vertices are reachable by directed paths
.buildDescendantsMatrix <- function(g) {
  stopifnot(edgemode(g) == "directed") ## QC

  dmat <- matrix(FALSE, nrow=numNodes(g), ncol=numNodes(g),
                 dimnames=list(nodes(g), nodes(g)))

  ## get the topological order of the vertices. assume that the sequence ontology
  ## graph is given such that the most general term SO:0001060 'sequence_variant'
  ## is the common sink of all terms. then traverse in reverse-topological order
  ## adding the children found
  to <- rev(tsort(g))
  edl <- edges(g)
  for (v in to) {
    if (length(edl[[v]]) > 0) {
      childdesc <- dmat[edl[[v]], , drop=FALSE]
      dmat[v, ] <- apply(childdesc, 2, function(x) any(x))
      dmat[v, edl[[v]]] <- TRUE
    }
  }

  dmat
}

## builds an ancestors matrix of an input acyclic digraph 'g'
## the resulting matrix has for every vertex in the rows a logical
## mask indicating what other vertices reach this vertex by directed paths
.buildAncestorsMatrix <- function(g) {
  stopifnot(edgemode(g) == "directed") ## QC

  amat <- matrix(FALSE, nrow=numNodes(g), ncol=numNodes(g),
                 dimnames=list(nodes(g), nodes(g)))

  ## get the topological order of the vertices. assume that the sequence ontology
  ## graph is given such that the most general term SO:0001060 'sequence_variant'
  ## is the common sink of all terms. then traverse in topological order
  ## adding the parents found
  to <- tsort(g)
  edl <- inEdges(g)
  for (v in to) {
    if (length(edl[[v]]) > 0) {
      parentsanc <- amat[edl[[v]], , drop=FALSE]
      amat[v, ] <- apply(parentsanc, 2, function(x) any(x))
      amat[v, edl[[v]]] <- TRUE
    }
  }

  amat
}

## returns SO IDs from the SO IDs or SO descriptions given in 'soterms'
## that are present in the SO graph 'gSO'
.findSOIDs <- function(soterms, gSO) {
  if (soterms[1] == "Any")
    return(nodes(gSO))

  mask <- soterms %in% nodes(gSO)
  if (any(!mask)) {
    solabels <- unlist(nodeData(gSO, nodes(gSO), "label"))
    mt <- match(soterms, solabels)
    soterms[!is.na(mt)] <- names(solabels[mt[!is.na(mt)]])
    mask <- soterms %in% nodes(gSO)
    if (any(!mask))
      warning(sprintf("%s are not valid sequence ontology terms.", paste(soterms[!mask], collapse=", ")))
  }
  soterms <- soterms[mask]
  soterms
}

## returns an edge list from a logical adjacency matrix
.aMat2edgeList <- function(amat, vtc) {
  stopifnot(!is.null(colnames(amat)) | !is.null(rownames(amat))) ## QC, amat should have row & column names

  wh <- which(amat, arr.ind=TRUE)
  wh[] <- colnames(amat)[wh]
  split(wh[, "col"], wh[, "row"])[vtc]
}

## returns the ancestor vertices from 'vtc' in the SO graph stored
## in the input VariantFilteringParam object
.ancestorsSO <- function(vfparam, vtc) {
  aedges <- .aMat2edgeList(soamat(vfparam), vtc)
  unique(unlist(aedges, use.names=FALSE))
}

## returns the descendant vertices from 'vtc' in the SO graph stored
## in the input VariantFilteringParam object
.descendantsSO <- function(vfparam, vtc) {
  dedges <- .aMat2edgeList(sodmat(vfparam), vtc)
  unique(unlist(dedges, use.names=FALSE))
}
