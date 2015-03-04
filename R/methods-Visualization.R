##
## visualization methods
##

## plot method
setMethod("plot", signature(x="VariantFilteringResults"),
          function(x, what, flank=20, ...) {
            require(Gviz)
            options(ucscChromosomeNames=FALSE)

            if (is.null(what) || missing(what))
              stop("Please set a value in the 'what' argument either as a variant identifier, a gene identifier, a chromosome name or a genomic position (as a GRanges object).")

            if (length(what) == 0)
              stop("The value in the 'what' argument has length zero.")

            fv <- filteredVariants(x)
            vars1 <- fv[[1]]
            seqlevelsStyle(vars1) <- seqlevelsStyle(param(x)$bsgenome) ## switch to the sequence style of the genome

            rng <- GRanges()
            if (class(what) == "character") {
              maskVarID <- vars1$VARID %in% what | vars1$dbSNP %in% what
              maskGene <- vars1$GENEID %in% what | vars1$GENE %in% what

              if (any(maskVarID | maskGene)) {
                rng <- range(vars1[maskVarID | maskGene])

                if (length(rng) == 0)
                  stop("None of the given values in the 'what' argument exists in the given 'VariantFilteringResults' object 'x'.")
                if (length(rng) > 1)
                  stop("The given values in the 'what' argument correspond to more than one chromosome.")
              } else {
                gr <- GRanges(seqnames=what[1], IRanges(1, 1))
                tryCatch({
                  seqlevelsStyle(gr) <- seqlevelsStyle(param(x)$bsgenome) ## switch to the sequence style of the genome
                }, error=function(err) {
                  stop(sprintf("The value %s given in the 'what' argument does not match any variant or gene identifer or chromosome name for %s.", what[1], organism(param(x)$bsgenome)))
                })

                chrlen <- seqlengths(param(x)$bsgenome)[as.character(seqnames(gr))]
                maskChrom <- as.vector(seqnames(vars1) %in% as.character(seqnames(gr)))
                if (any(maskChrom))
                  rng <- GRanges(seqnames=as.character(seqnames(gr)), IRanges(start=1, end=chrlen))
                else
                  stop(sprintf("There are no variants annotated to chromosome %s (%s %s).",
                               what[1], provider(param(x)$bsgenome), as.character(seqnames(gr))))
              }

            } else if (class(what) == "GRanges") {

              seqlevelsStyle(what) <- sealevelsStyle(param(x)$bsgenome) ## switch to the sequence style of the genome
              ov <- findOverlaps(what, vars1)
              if (length(ov) == 0)
                stop("The genomic ranges in the 'what' argument do not overlap any variant in the given 'VariantFilteringResults' object 'x'.")

              rng <- range(vars1[subjectHits(ov)])

              if (length(rng) > 1)
                stop("The given genomic ranges in the 'what' argument correspond to more than one chromosome.")

            } else
              stop("The value of the 'what' argument must be either a character vector or a 'GRanges' object.")

            if (class(what) != "GRanges" && flank > 0) {
              chrlen <- seqlengths(param(x)$bsgenome)[as.character(seqnames(rng))]
              if (start(rng) > flank && end(rng) < chrlen-flank)
                rng <- flank(rng, width=flank, both=TRUE)
              else {
                if (start(rng)-flank < 1)
                  start(rng) <- 1
                else
                  rng <- resize(rng, width=width(rng)+flank, fix="end")

                if (end(rng)+flank > chrlen)
                  end(rng) <- chrlen
                else
                  rng <- resize(rng, width=width(rng)+flank, fix="start")
              }

            }

            ov <- findOverlaps(rng, vars1)
            vars1 <- vars1[subjectHits(ov)]
            txdb <- param(x)$txdb
            txdbmdata <- metadata(param(x)$txdb)
            geneannot <- genes(txdb)
            geneannot$symbol <- select(param(x)$orgdb, keys=names(geneannot), columns="SYMBOL")$SYMBOL
            seqlevelsStyle(geneannot) <- seqlevelsStyle(param(x)$bsgenome) ## switch to the sequence style of the genome

            tracks <- GenomeAxisTrack()
            tracks <- c(tracks,
                        GeneRegionTrack(geneannot, genome=genome(param(x)$bsgenome)[1],
                                        chromosome=as.character(seqnames(vars1))[1], shape="arrow",
                                        name=paste(txdbmdata[txdbmdata$name %in% "Data source", "value"], "Genes"),
                                        transcriptAnnotation="symbol"))
            tracks <- c(tracks,
                        AnnotationTrack(start=start(vars1), width=width(vars1), chromosome=as.character(seqnames(vars1)),
                                        strand=strand(vars1), genome=genome(param(x)$bsgenome)[1], shape="box",
                                        id=sprintf("%s\n(%s)", vars1$VARID, vars1$DESC), name="Variants"))

            tracks <- c(tracks, SequenceTrack(param(x)$bsgenome, chromosome=as.character(seqnames(vars1))[1]))


            plotTracks(tracks, chromosome=as.vector(seqnames(vars1))[1],
                       from=start(rng), to=end(rng), featureAnnotation="id", fontcolor.feature="black", ...)
            
          })


strandedBamImport <- function (file, selection) {
    if (!file.exists(paste(file, "bai", sep = ".")))
        stop("Unable to find index for BAM file '", file, "'. You can
build an index using the following command:\n\t",
             "library(Rsamtools)\n\tindexBam(\"", file, "\")")
    sinfo <- scanBamHeader(file)[[1]]
    res <- if (!as.character(seqnames(selection)[1]) %in%
               names(sinfo$targets)) {
        mcols(selection) <- DataFrame(score = 0)
        selection
    }else {
        param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                              which = selection, flag =
scanBamFlag(isUnmappedQuery = FALSE))
        x <- scanBam(file, param = param)[[1]]
        gr <- GRanges(strand=x[["strand"]], ranges=IRanges(x[["pos"]],
width = x[["qwidth"]]), seqnames=seqnames(selection)[1])
        grs <- split(gr, strand(gr))
        cov <- lapply(grs[c("+", "-")], function(y) coverage(ranges(y),
width=end(selection)))
        pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y),
end(y))))))
        if(length(pos)==0){
            mcols(selection) <- DataFrame(plus=0, minus=0)
            selection
        }else{
            GRanges(seqnames = seqnames(selection)[1],
ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)),
                    plus=as.numeric(cov[["+"]][head(pos, -1)]),
                    minus=-as.numeric(cov[["-"]][head(pos, -1)]))
        }
    }
    return(res)
}

