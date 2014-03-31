\name{deNovo}
\alias{deNovo}
\alias{deNovo,VariantFilteringParam-method}

\title{
De Novo variants analysis
}
\description{
This function has been created in order to search for \emph{de novo} variants in one individual, discarding the ones shared with his/her parents.
}
\usage{
\S4method{deNovo}{VariantFilteringParam}(param,
                                         BPPARAM=bpparam())
}
\arguments{
  \item{param}{
  An VariantFilteringParam object containing VCF file(s). 1 independent files for affected individuals and 2 more for unaffected ones (both parents). If the VCF is a multi-sample, it can contain the same amount of individuals divided likewise.  }
  \item{BPPARAM}{
  An object of class 'BiocParallelParam' specifiying parameters related to
  the parallel execution of some of the tasks and calculations within this function.
  See function \code{\link[BiocParallel]{bpparam}()} from the \code{BiocParallel} package.
  }
}
\details{
This function requires as an input a VariantFilteringParam class object which contains the vcf files ready for the analysis, along with a ped file which specifies the characteristics of each individual present in the analysis.

Vcf from both parents are required and only one child at time can be analyzed.
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
\code{\link{autosomalDominant}}
\code{\link{autosomalRecessiveHeterozygous}}
\code{\link{xLinked}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.gz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped")
param <- VariantFilteringParam(vcfFileNames=CEUvcf,
                     pedFileName=CEUped)
deNo <- deNovo(param)
deNo
}
}
\keyword{utilities}
