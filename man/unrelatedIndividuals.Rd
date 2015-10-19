\name{unrelatedIndividuals}
\alias{unrelatedIndividuals}
\alias{unrelatedIndividuals,VariantFilteringParam-method}

\title{
Analysis designed to be applied over a group of unrelated individuals
}
\description{
This function is designed to create an object to deepen into the
variants presented by a group of unrelated individuals
}
\usage{
\S4method{unrelatedIndividuals}{VariantFilteringParam}(param,
                                                       svparam=ScanVcfParam(),
                                                       BPPARAM=bpparam("SerialParam"))
}
\arguments{
  \item{param}{
  A \code{\link{VariantFilteringParam}} object built from a multisample VCF file.
  }
  \item{svparam}{
  An instance of a \code{ScanVcfParam} object to enable analyzing a subset of variants and samples. This object
  is passed internally to a call to the \code{\link[VariantAnnotation]{readVcf}()} function in the
  \code{VariantAnnotation} package, see its help page for a complete description of this functionality.
  }
  \item{BPPARAM}{
  An object of class \code{\link[BiocParallel]{BiocParallelParam}} specifiying parameters related to
  the parallel execution of some of the tasks and calculations within this function.
  See function \code{\link[BiocParallel]{bpparam}()} from the \code{BiocParallel} package.
  }
}
\details{
This function requires as an input a \code{\link{VariantFilteringParam}} class object built from an input multisample VCF file.
}
\value{
An object of class \code{\link{VariantFilteringResults}} including functional annotations on all variants.
}
\references{
Elurbe D.M., Mila, M., Castelo, R. The VariantFiltering package, in preparation.
}
\author{
Dei M. Elurbe and R. Castelo
}
\seealso{
\code{\link{autosomalRecessiveHomozygous}}
\code{\link{autosomalRecessiveHeterozygous}}
\code{\link{autosomalDominant}}
\code{\link{deNovo}}
\code{\link{xLinked}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
param <- VariantFilteringParam(vcfFilenames=CEUvcf)
uInd <- unrelatedIndividuals(param)
uInd
}
}
\keyword{utilities}
