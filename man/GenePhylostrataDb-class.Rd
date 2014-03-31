\name{GenePhylostrataDb-class}
\docType{class}
\alias{GenePhylostrataDb-class}
\alias{GenePhylostrataDb}

\alias{annotateVariants,GenePhylostrataDb-method}
\alias{genePhylostrata,GenePhylostrataDb-method}
\alias{genePhylostratum,GenePhylostrataDb,missing-method}
\alias{genePhylostratum,GenePhylostrataDb,character-method}
\alias{GenePhylostrataDb}
\alias{genePhylostrata}
\alias{genePhylostratum}
\alias{organism,GenePhylostrataDb-method}
\alias{humanGenesPhylostrata}
\alias{show,GenePhylostrataDb-method}

\title{PhastConsDb class}

\description{
Class for storing gene-level conservation information in the form of levels of phylogenetic strata; see Neme and Tautz (2013).
}
\details{
The \code{GenePhylostrataDb} class and associated methods serve the purpose of storing and manipulating gene-level
conservation information in the form of levels of phylogenetic strata (Neme and Tautz, 2013).
One such objects is created at loading time by the \code{VariantFiltering} package with the constructor function
\code{GenePhylostrataDb()}, and it is called \code{humanGenesPhylostrata}.

}
\references{
  Neme, R. and Tautz, D. Phylogenetic patterns of emergence of new genes support a model of frequent de novo evolution. BMC Genomics, 14:117, 2013
}
\value{.}
\author{R. Castelo}
\seealso{
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
}
\examples{
humanGenesPhylostrata
}
\source{
  \url{http://genomebiology.com/content/supplementary/1471-2164-14-117-s1.xlsx}
}
\keyword{datasets}
