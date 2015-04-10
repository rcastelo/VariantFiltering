##
## visualization methods
##

## plot method
setMethod("plot", signature(x="VariantFilteringResults"),
          function(x, what, sampleName, flankingNt=20, showAlnNtCutoff=200, isPaired=FALSE, ...) {

            options(ucscChromosomeNames=FALSE)

            if (is.null(what) || missing(what))
              stop("Please set a value in the 'what' argument either as a variant identifier, a gene identifier, a chromosome name or a genomic position (as a GRanges object).")

            if (length(what) == 0)
              stop("The value in the 'what' argument has length zero.")

            bsgenome <- param(x)$bsgenome

            ## need to adapt to the seqlevelsStyle from the BAM file, when available
            leadseqlevelsstyle <- seqlevelsStyle(bsgenome)
            if (!missing(sampleName)) {
              if (!is.character(sampleName))
                stop("The argument 'sampleName' should be a chracter string.")

              if (length(sampleName) > 1)
                sampleName <- sampleName[1]

              if (!sampleName %in% samples(x))
                stop(sprintf("%s does not form part of the samples in the input object.", sampleName))

              if (!sampleName %in% rownames(bamSamples(bamFiles(x))))
                stop(sprintF("There is no BAM file associated to the sample %s. Use 'bamFiles()' to add it.", sampleName))

              whbam <- which(rownames(bamSamples(bamFiles(x))) %in% sampleName)
              hdr <- scanBamHeader(bamPaths(bamFiles(x))[whbam])
              leadseqlevelsstyle <- seqlevelsStyle(names(hdr[[1]]$targets))
              seqlevelsStyle(bsgenome) <- leadseqlevelsstyle
            }

            fv <- filteredVariants(x)
            vars1 <- fv[[1]]
            seqlevelsStyle(vars1) <- leadseqlevelsstyle

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
                  seqlevelsStyle(gr) <- leadseqlevelsstyle
                }, error=function(err) {
                  stop(sprintf("The value %s given in the 'what' argument does not match any variant or gene identifer or chromosome name for %s.", what[1], organism(bsgenome)))
                })

                chrlen <- seqlengths(bsgenome)[as.character(seqnames(gr))]
                maskChrom <- as.vector(seqnames(vars1) %in% as.character(seqnames(gr)))
                if (any(maskChrom))
                  rng <- GRanges(seqnames=as.character(seqnames(gr)), IRanges(start=1, end=chrlen))
                else
                  stop(sprintf("There are no variants annotated to chromosome %s (%s %s).",
                               what[1], provider(bsgenome), as.character(seqnames(gr))))
              }

            } else if (class(what) == "GRanges") {

              seqlevelsStyle(what) <- leadseqlevelsstyle
              ov <- findOverlaps(what, vars1)
              if (length(ov) == 0)
                stop("The genomic ranges in the 'what' argument do not overlap any variant in the given 'VariantFilteringResults' object 'x'.")

              rng <- range(vars1[subjectHits(ov)])

              if (length(rng) > 1)
                stop("The given genomic ranges in the 'what' argument correspond to more than one chromosome.")

            } else
              stop("The value of the 'what' argument must be either a character vector or a 'GRanges' object.")

            if (class(what) != "GRanges" && flankingNt > 0) {
              chrlen <- seqlengths(bsgenome)[as.character(seqnames(rng))]
              if (start(rng)-flankingNt < 1)
                start(rng) <- 1
              else
                rng <- resize(rng, width=width(rng)+flankingNt, fix="end")

              if (end(rng)+flankingNt > chrlen)
                end(rng) <- chrlen
              else
                rng <- resize(rng, width=width(rng)+flankingNt, fix="start")
            }

            ## fetch all variants that can be anntotated in the plotting region
            ov <- findOverlaps(rng, vars1)
            vars1 <- vars1[subjectHits(ov)]
            chr2display <- as.character(seqnames(vars1))[1]

            ## prepare annotations for display, take care of restoring later the
            ## seqlevelsStyle of the TxDb object
            txdb <- param(x)$txdb
            seqlevelsstyletxdb <- seqlevelsStyle(txdb)
            seqlevelsStyle(txdb) <- leadseqlevelsstyle
            seqlevels(txdb, force=TRUE) <- chr2display
            txdbmdata <- metadata(txdb)

            ## geneannot <- exonsBy(txdb, by="gene")
            ## flatgeneannot <- unlist(geneannot)
            ## flatgeneannot$symbol <- select(param(x)$orgdb, keys=names(flatgeneannot), columns="SYMBOL")$SYMBOL
            ## seqlevelsStyle(flatgeneannot) <- leadseqlevelsstyle
            ## flatgeneannot <- flatgeneannot[seqnames(flatgeneannot) %in% as.character(seqnames(vars1))[1]]
            ## geneannot <- split(flatgeneannot, names(flatgeneannot))

            ## seqlevelsStyle(geneannot) <- leadseqlevelsstyle

            tracks <- GenomeAxisTrack()
            tracks <- c(tracks,
                        GeneRegionTrack(txdb,
                                        name=paste(txdbmdata[txdbmdata$name %in% "Data source", "value"], "Genes")))
            restoreSeqlevels(txdb)
            seqlevelsStyle(txdb) <- seqlevelsstyletxdb

            ## tracks <- c(tracks,
            ##             GeneRegionTrack(geneannot, genome=genome(bsgenome)[1],
            ##                             chromosome=as.character(seqnames(vars1))[1], shape="arrow",
            ##                             name=paste(txdbmdata[txdbmdata$name %in% "Data source", "value"], "Genes"),
            ##                             transcriptAnnotation="symbol"))

            strack <- SequenceTrack(bsgenome, chromosome=as.character(seqnames(vars1))[1])

            if (width(rng) < showAlnNtCutoff && !missing(sampleName)) {
              whbam <- which(rownames(bamSamples(bamFiles(x))) %in% sampleName)
              tracks <- c(tracks,
                          AlignmentsTrack(range=bamPaths(bamFiles(x))[whbam], isPaired=isPaired,
                                          referenceSequence=strack, showMismatches=TRUE))
            } else {
              idlabel <- ifelse(is.na(vars1$HGVSc), sprintf("%s(%s)", vars1$VARID, vars1$HGVSg),
                                sprintf("%s\n(%s)", vars1$VARID, vars1$HGVSc))
              tracks <- c(tracks,
                          AnnotationTrack(start=start(vars1), width=width(vars1),
                                          chromosome=as.character(seqnames(vars1)),
                                          strand=strand(vars1), genome=genome(bsgenome)[1],
                                          shape="box", id=idlabel, name="Variants"))
            }

            tracks <- c(tracks, strack)

            plotTracks(tracks, chromosome=as.character(seqnames(vars1))[1],
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

