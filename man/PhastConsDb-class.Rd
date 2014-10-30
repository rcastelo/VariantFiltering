\name{PhastConsDb-class}
\docType{class}
\alias{PhastConsDb-class}
\alias{PhastConsDb}

\alias{PhastConsDb}
\alias{annotateVariants,PhastConsDb-method}
\alias{organism,PhastConsDb-method}
\alias{provider,PhastConsDb-method}
\alias{providerVersion,PhastConsDb-method}
\alias{referenceGenome,PhastConsDb-method}
\alias{scores,PhastConsDb,GRanges-method}
\alias{scores}
\alias{seqinfo,PhastConsDb-method}
\alias{seqlengths,PhastConsDb-method}
\alias{seqnames,PhastConsDb-method}
\alias{show,PhastConsDb-method}

\title{PhastConsDb class}

\description{
Class for annotation packages storing UCSC phastCons conservation scores.
}
\usage{
\S4method{annotateVariants}{PhastConsDb}(annObj, variantsGR, BPPARAM)
\S4method{scores}{PhastConsDb,GRanges}(object, gpos,
                                       summaryFun="mean",
                                       coercionFun="as.numeric",
                                       caching=TRUE)
}
\arguments{
  \item{annObj}{A \code{PhastConsDb} object.}
  \item{variantsGR}{A \code{GRanges} object with the variants to annotate.}
  \item{BPPARAM}{An object of class \code{BiocParallelParam} specifying parameters related to
    the parallel execution of this function. See function \code{\link[BiocParallel]{bpparam}()}
    from the \code{BiocParallel} package.}
  \item{object}{A \code{PhastConsDb} object.}
  \item{gpos}{A \code{GRanges} object with positions from where to retrieve phastCons scores.}
  \item{summaryFun}{Function to summarize phastCons scores when more than one position is retrieved.
                    By default is set to the arithmetic mean.}
  \item{coercionFun}{Function to coerce the stored phastCons scores, before the summary function is
                     applied. By default phastCons scores are coerced to real values.}
  \item{caching}{Flag setting whether phastCons scores per chromosome should be kept cached in memory
                 (\code{TRUE}, default) or not (\code{FALSE}). The latter option minimizes the memory
                 footprint but slows down the performance when the \code{scores()} method is called
                 multiple times.}
}
\details{
The \code{PhastConsDb} class and associated methods serve the purpose of creating annotation packages that store
phastCons nucleotide-level conservation scores from the UCSC Genome Browser. One such annotation package is:

\tabular{ll}{
  \code{phastCons100way.UCSC.hg19} \tab Nucleotide-level phastCons conservation scores from the UCSC Genome Browser
        downloaded in November 2013 from \url{http://genome.ucsc.edu} and calculated from 100 vertebrate species. \cr
}

}
\value{.}
\author{R. Castelo}
\seealso{
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
}
\examples{
if (require(phastCons100way.UCSC.hg19)) {
  library(GenomicRanges)

  phastCons100way.UCSC.hg19
  scores(phastCons100way.UCSC.hg19,
         GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
}
}
\source{
  \url{http://genome.ucsc.edu}
}
\keyword{datasets}
