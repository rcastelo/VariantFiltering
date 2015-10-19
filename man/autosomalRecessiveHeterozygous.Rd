\name{autosomalRecessiveHeterozygous}
\alias{autosomalRecessiveHeterozygous}
\alias{autosomalRecessiveHeterozygous,VariantFilteringParam-method}
\title{
Autosomal recessive inheritance analysis: Heterozygous
}
\description{
This function aims to analyze the variants of the unaffected individuals, storing and grouping the heterozygous ones by gene. The affected individuals ought present two or more different heterozygous changes in the gene, and at least one of them shall come from each parent.
}

\usage{
\S4method{autosomalRecessiveHeterozygous}{VariantFilteringParam}(param,
                                                                 svparam=ScanVcfParam(),
                                                                 BPPARAM=bpparam("SerialParam"))
}
\arguments{
  \item{param}{
  A \code{\link{VariantFilteringParam}} object built from a multisample VCF file with at least one affected individual and two or more unaffected ones, and from a PED file specyfing the family relationships among individuals as well as their gender and phenotype status (affected or unaffected).
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
\code{\link{autosomalDominant}}
\code{\link{xLinked}}
\code{\link{deNovo}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped")
param <- VariantFilteringParam(vcfFilenames=CEUvcf,
                     pedFilename=CEUped)
reHet <- autosomalRecessiveHeterozygous(param)
reHet
}
}
\keyword{utilities}
