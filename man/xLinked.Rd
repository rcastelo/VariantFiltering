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
                                          svparam=ScanVcfParam(),
                                          BPPARAM=bpparam("SerialParam"))
}
\arguments{
  \item{param}{
  A \code{\link{VariantFilteringParam}} object built from a multisample VCF file with at least one affected individual and zero or more unaffected ones, and from a PED file specyfing the family relationships among individuals as well as their gender and phenotype status (affected or unaffected).
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
This function requires as an input a \code{\link{VariantFilteringParam}} class object built from an input multisample VCF file, along with a PED file.
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
