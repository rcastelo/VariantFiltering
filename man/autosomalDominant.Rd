\name{autosomalDominant}
\alias{autosomalDominant}
\alias{autosomalDominant,VariantFilteringParam-method}

\title{
Autosomal dominant inheritance analysis
}
\description{
This function identifies variants present in all the affected individual(s) discarding the ones that also occur in at least one of the unaffected subjects.
}
\usage{
\S4method{autosomalDominant}{VariantFilteringParam}(param,
                                                    BPPARAM=bpparam())
}
\arguments{
  \item{param}{
  An VariantFilteringParam object containing VCF file(s). From 1 to 5 independent files for affected individuals and 0 to 5 more for unaffected ones (up to 10 individuals). If the VCF is a multi-sample, it can contain the same amount of individuals divided likewise.
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
\code{\link{xLinked}}
\code{\link{deNovo}}
\code{\link{allInheritanceModels}}
\code{\link{unrelatedIndividuals}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.gz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped")
param <- VariantFilteringParam(vcfFileNames=CEUvcf,
                     pedFileName=CEUped)
aDo <- autosomalDominant(param)
aDo
}
}
\keyword{utilities}
