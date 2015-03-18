\name{xLinked}
\alias{xLinked}
\alias{xLinked,VariantFilteringParam-method}

\title{
X-Linked inheritance analysis
}
\description{
This function identifies variants that appear only in the X chromosome of the unaffected females as heterozygous, don't appear in the unaffected males analyzed and finally are present (as homozygous) in the affected male(s).
}

\usage{
\S4method{xLinked}{VariantFilteringParam}(param,
                                          BPPARAM=bpparam())
}
\arguments{
  \item{param}{
  An VariantFilteringParam object containing VCF file(s). From 1 to 5
  independent files for affected individuals, 0 to 5 more for the
  carrier females and 0 to 5 additional individuals as unaffected males
  (up to 15 individuals). If the VCF is a multi-sample, it can contain
  the same amount of individuals divided likewise.
  }
  \item{BPPARAM}{
  An object of class 'BiocParallelParam' specifiying parameters related to
  the parallel execution of some of the tasks and calculations within this function.
  See function \code{\link[BiocParallel]{bpparam}()} from the \code{BiocParallel} package.
  }
}
\details{
This function requires as an input a VariantFilteringParam class object which contains the vcf files ready for the analysis, along with a ped file which specifies the characteristics of each individual present in the analysis.
}
\value{
An object of class \code{\link{VariantFilteringResults}} including functional annotations on all selected variants.
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
\code{\link{allInheritanceModels}}
\code{\link{unrelatedIndividuals}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

## This actually won't run b/c in this trio de descendant is a female
CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped")
param <- VariantFilteringParam(vcfFilenames=CEUvcf,
                     pedFilename=CEUped)
xlid <- xLinked(param)
xlid
}
}
\keyword{utilities}
