\name{autosomalRecessiveHomozygous}
\alias{autosomalRecessiveHomozygous}
\alias{autosomalRecessiveHomozygous,VariantFilteringParam-method}

\title{
Autosomal recessive inheritance analysis: Homozygous
}
\description{
This function works analyzing the variants of the unaffected individuals storing the common heterozygous ones and comparing them with the common homozygous variants between the affected subjects.
}
\usage{
\S4method{autosomalRecessiveHomozygous}{VariantFilteringParam}(param,
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
This function requires as an input an \code{VariantFilteringParam} class object which contains the VCF files ready for the analysis, along with a ped file which specifies the characteristics of each individual present in the analysis.
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
\code{\link{autosomalRecessiveHeterozygous}}
\code{\link{autosomalDominant}}
\code{\link{xLinked}}
\code{\link{deNovo}}
\code{\link{allInheritanceModels}}
\code{\link{unrelatedIndividuals}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{
library(VariantFiltering)

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"),
                    "CEUtrio.vcf.gz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"),
                    "CEUtrio.ped")
param <- VariantFilteringParam(vcfFileNames=CEUvcf, pedFileName=CEUped)
reHo <- autosomalRecessiveHomozygous(param)
reHo
}
}
\keyword{utilities}
