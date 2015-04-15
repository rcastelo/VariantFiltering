\name{VariantFiltering-package}
\alias{VariantFiltering-package}
\alias{VariantFiltering}
\docType{package}
\title{
Filtering of coding and non-coding genetic variants
}
\description{
The \code{VariantFiltering} package filters coding and non-coding genetic variants using different
criteria such as an inheritance model (autosomal recessive -both, homozygous and heterozygous-,
autosomal dominant, X-linked and \emph{de novo}), amino acid change consequence, minimum allele
frequencies, cryptic splice site potential, conservation, etc.
}
\section{Functions}{
  \itemize{
    \item \code{\link{autosomalRecessiveHomozygous}} identify homozygous variants in the affected individual(s) while the unaffected ones present these same variants but in heterozygous state. Autosomal recessive inheritance pattern.
    \item \code{\link{autosomalRecessiveHeterozygous}} identify variants grouped by genes with two (or more) heterogeneous alleles (at least one on each allele, i.e. coming from each parent). Autosomal recessive inheritance pattern.
    \item \code{\link{autosomalDominant}} identify variants present in all the affected individual(s) discarding the ones that also occur in at least one of the unaffected subjects. Autosomal dominant inheritance pattern.
    \item \code{\link{xLinked}} identify variants that appear only in the X chromosome of the unaffected females as heterozygous, don't appear in the unaffected males analyzed and finally are present (as homozygous) in the affected male(s). X-linked inheritance pattern.
    \item \code{\link{deNovo}} identify variants in the affected
    individual that have not been inherited.
    \item\code{\link{unrelatedIndividuals}} annotate variants without filtering by any inheritance pattern.
  }
}
\author{
Dei M. Elurbe and Robert Castelo.

Maintainer: Robert Castelo <robert.castelo@upf.edu>
}

\references{
Elurbe D.M., Mila, M., Castelo, R. VariantFiltering: filtering of coding and non-coding genetic variants, in preparation.
}

\keyword{package}
