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
                                                                 BPPARAM=bpparam())
}
\arguments{
  \item{param}{
  An VariantFilteringParam object containing VCF file(s). From 1 to 5 independent files for affected individuals and 2 more for the carriers of each allele (required). If the VCF is a multi-sample, it can contain the same amount of individuals divided likewise.
  }
  \item{BPPARAM}{
  An object of class 'BiocParallelParam' specifiying parameters related to
  the parallel execution of some of the tasks and calculations within this function.
  See function \code{\link[BiocParallel]{bpparam}()} from the \code{BiocParallel} package.
  }
}
\details{
This function requires as an input a VariantFilteringParam class object which contains the vcf files ready for the analysis, along with a ped file which specifies the characteristics of each individual present in the analysis.

For this function, it is not possible to work with all the transcripts and it is limited to variants located into coding regions.
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
\code{\link{allInheritanceModels}}
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
