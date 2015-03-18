\name{allInheritanceModels}
\alias{allInheritanceModels}
\alias{allInheritanceModels,VariantFilteringParam-method}

\title{
Analysis designed to be applied over a group of related individuals
}
\description{
This method filters variants from a group of related individuals
annotating compatible inheritance models of segregation.
}
\usage{
\S4method{allInheritanceModels}{VariantFilteringParam}(param,
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
\code{\link{autosomalRecessiveHomozygous}}
\code{\link{autosomalRecessiveHeterozygous}}
\code{\link{autosomalDominant}}
\code{\link{deNovo}}
\code{\link{xLinked}}
\code{\link{allInheritanceModels}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
param <- VariantFilteringParam(vcfFilenames=CEUvcf)
aim <- allInheritanceModels(param)
aim
}
}
\keyword{utilities}
