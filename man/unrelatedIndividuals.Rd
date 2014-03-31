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
                                                       BPPARAM=bpparam())
}
\arguments{
  \item{param}{
  An VariantFilteringParam object containing a VCF file.
  }
  \item{BPPARAM}{
  An object of class 'BiocParallelParam' specifiying parameters related to
  the parallel execution of some of the tasks and calculations within this function.
  See function \code{\link[BiocParallel]{bpparam}()} from the \code{BiocParallel} package.
  }
}
\details{
This function requires as an input a VariantFilteringParam class object which contains the vcf file ready for the analysis.
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
\code{\link{allInheritanceModels}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.gz")
param <- VariantFilteringParam(vcfFileNames=CEUvcf)
uInd <- unrelatedIndividuals(param)
uInd
}
}
\keyword{utilities}
