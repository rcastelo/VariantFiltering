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
                                                    svparam=ScanVcfParam(),
                                                    use=c("everything", "complete.obs", "all.obs"),
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
  \item{use}{character string specifying the policy to apply on missing genotypes when comparing them. This policy can be either '"everything"' (default), '"complete.obs"' or '"all.obs"'. The default policy ('"everything"') will propagate NA truth values using the behavior of the R logical operators, with the exception that when the final truth value associated with a variant is \code{NA}, the variant is ultimately discarded.
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
An object of class \code{\link{VariantFilteringResults}} including functional annotations on variants segregating according to an autosomal dominant inheritance model.
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
\code{\link{unrelatedIndividuals}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped")
param <- VariantFilteringParam(vcfFilename=CEUvcf,
                     pedFilename=CEUped)
aDo <- autosomalDominant(param)
aDo
}
}
\keyword{utilities}
