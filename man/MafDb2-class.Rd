\name{MafDb2-class}
\docType{class}
\alias{MafDb2-class}
\alias{MafDb2}

\alias{annotateVariants,MafDb2-method}
\alias{mafByOverlaps,MafDb2-method}
%%\alias{mafByOverlaps}
\alias{mafById,MafDb2-method}
%%\alias{mafById}
\alias{populations,MafDb2-method}
%%\alias{populations}
\alias{organism,MafDb2-method}
\alias{provider,MafDb2-method}
\alias{providerVersion,MafDb2-method}
\alias{referenceGenome,MafDb2-method}
\alias{seqinfo,MafDb2-method}
\alias{seqlengths,MafDb2-method}
\alias{seqnames,MafDb2-method}
\alias{seqlevelsStyle,MafDb2-method}
\alias{show,MafDb2-method}
\alias{$,MafDb2-method}

\title{MafDb2 class}

\description{
Class for annotation packages storing minor allele frequency data.
}
\usage{
\S4method{mafByOverlaps}{MafDb2}(x, ranges, pop, caching)
\S4method{mafById}{MafDb2}(x, ids, pop, caching)
\S4method{populations}{MafDb2}(x)
}
\arguments{
  \item{x}{A \code{MafDb2} object.}
  \item{ranges}{Either a \code{GRanges} object, a \code{GPos} object or a character string vector with the format \code{"CHR:START[-END]"}.}
  \item{ids}{A character string vector with variant identifiers annotated by the MAF data source, typically dbSNP 'rs' identifiers. Note that the mapping of these identifiers to genomic positions and MAF values might be a subset of the most up to date dbSNP 'rs' identifier assignment to variants. To access the latter, please use the \code{snpsById()} method from the \code{BSgenome} package with the desired \code{SNPlocs.*} package.}
  \item{pop}{A character string vector with the populations for which we want to retrieve MAF values.}
  \item{caching}{logical; TRUE (default) indicates that the function stores into main memeory the MAF data as it gets loaded from disk, improving performance; FALSE forces this function to load MAF data from disk each time, decreasing performance and memory requirements.}
}
\details{
The \code{MafDb2} class and associated methods serve the purpose of creating
annotation packages that store minor allele frequency data that can be queried by genomic position. Two such annotation packages are:

\tabular{ll}{
  \code{MafDb.ExAC.r0.3.1.snvs.hs37d5} \tab MAF values from the ExAC consortium. \cr
  \code{MafDb.ExAC.r0.3.1.nonTCGA.snvs.hs37d5} \tab MAF values from the ExAC consortium. \cr
}

This object class tries to reduce the disk space required to store allele frequencies (AFs)
for millions of SNPs by coding AF float values, which range between 0 and 1, into a single-byte
\code{raw} object type. To achieve this, the original AF values are rounded to one significant
digit for AF < 0.01 and two significant digits for AF >= 0.01. When a variant has multiple
alternate alleles, only the larger MAF value is stored.
}
\author{R. Castelo}
\examples{

  ## lookup allele frequencies for rs1129038, a SNP associated to blue and brown eye colors
  ## as reported by Eiberg et al. Blue eye color in humans may be caused by a perfectly associated
  ## founder mutation in a regulatory element located within the HERC2 gene inhibiting OCA2 expression.
  ## Human Genetics, 123(2):177-87, 2008 [http://www.ncbi.nlm.nih.gov/pubmed/18172690]

  if (require(MafDb.ExAC.r0.3.1.snvs.hs37d5)) {
    mafdb <- MafDb.ExAC.r0.3.1.snvs.hs37d5
    mafdb

    ## specialized interface
    populations(mafdb)

    rng <- GRanges("15", IRanges(28356859, 28356859))
    mafByOverlaps(mafdb, rng)
    mafByOverlaps(mafdb, "15:28356859-28356859")
    mafByOverlaps(mafdb, "15:28356859")
    mafById(mafdb, "rs1129038")
  }

}

\source{
  \url{http://exac.broadinstitute.org}
}
\keyword{datasets}
