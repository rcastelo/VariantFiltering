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
  An VariantFilteringParam object containing VCF file(s). From 1 to 5 independent files for affected individuals and 0 to 5 more for unaffected ones (up to 10 individuals). If the VCF is a multi-sample, it can contain the same amount of individuals divided likewise.
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
\code{\link{autosomalRecessiveHeterozygous}}
\code{\link{autosomalDominant}}
\code{\link{xLinked}}
\code{\link{deNovo}}
\code{\link{unrelatedIndividuals}}
\code{\link{VariantFilteringResults}}
}
\examples{
\dontrun{
library(VariantFiltering)

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"),
                    "CEUtrio.vcf.bgz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"),
                    "CEUtrio.ped")
param <- VariantFilteringParam(vcfFilename=CEUvcf, pedFilename=CEUped)
reHo <- autosomalRecessiveHomozygous(param)
reHo
}
}
\keyword{utilities}
