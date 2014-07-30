\name{WeightMatrix-class}
\docType{class}
\alias{WeightMatrix-class}
\alias{WeightMatrix}

\alias{width,WeightMatrix-method}
\alias{conservedPositions,WeightMatrix-method}
\alias{show,WeightMatrix-method}
\alias{wmScore,WeightMatrix,DNAStringSet-method}
\alias{wmScore,WeightMatrix,character-method}

\alias{conservedPositions}
\alias{readWm}
\alias{wmAcceptorSites}
\alias{wmDonorSites}
\alias{wmScore}

\title{Weight matrix class}

\description{
Class for storing weight matrices that VariantFiltering uses to score potential cryptic splice sites.
}
\usage{
\S4method{width}{WeightMatrix}(x)
\S4method{conservedPositions}{WeightMatrix}(x)
\S4method{show}{WeightMatrix}(object)
\S4method{wmScore}{WeightMatrix,DNAStringSet}(object, dnaseqs)
\S4method{wmScore}{WeightMatrix,character}(object, dnaseqs)
}
\arguments{
  \item{x}{A \code{WeightMatrix} object.}
  \item{object}{A \code{WeightMatrix} object.}
  \item{dnaseqs}{Either a vector of character strings a \code{DNAStringSet} object, both of which store
        nucleotide sequences to be scored using the input \code{WeightMatrix} object.}
}
\details{
The \code{WeightMatrix} class and associated methods serve the purpose of enabling the \code{VariantFiltering} package
to score synonymous and intronic genetic variants for potential cryptic splice sites. The class and the methods,
however, are exposed to the end user since they could be useful for other analysis purposes.

The \code{VariantFiltering} package contains two weight matrices, one for 5'ss and another for 3'ss, which have been built
by a statistical method that accounts for dependencies between the splice site positions, minimizing the rate of
false positive predictions. The method concretely builds these models by inclusion-driven learning of Bayesian
networks and further details can be found in the paper of Castelo and Guigo (2004).

The function \code{readWm()} reads a weight matrix stored in a text file in a particular format and returns
a \code{WeightMatrix} object. See the \code{.ibn} files located in the \code{extdata} folder of the \code{VariantFiltering}
package, as an example of this format.

The method \code{wmScore()} scores one or more sequences of nucleotides using the input \code{WeightMatrix} object.
If the sequences are longer than the width of the weight matrix, this function will score every possible site
within those sequences. It returns a vector of with the calculated scores. When the scores cannot be calculated
because of a conserved position that does not occur in the sequence (i.e., absence of a GT dinucleotide with the
5'ss weight matrix), it returns \code{NA} as corresponding score value.

The method \code{width()} takes a \code{WeightMatrix} object as input and returns the number of positions of the
weight matrix.

The method \code{conservedPositions()} takes a \code{WeightMatrix} object as input and returns the number of
fully conserved positions in the weight matrix.
}
\value{.}
\references{
Castelo, R and Guigo, R. Splice site identification by idlBNs. Bioinformatics, 20(1):i69-i76, 2004.
}
\author{R. Castelo}
\examples{
wm <- readWm(file.path(system.file("extdata", package="VariantFiltering"), "hsap.donors.hcmc10_15_1.ibn"))
width(wm)
conservedPositions(wm)
wmScore(wm, "CAGGTAGGA")
wmScore(wm, "CAGGAAGGA")
wmScore(wm, "CAGGTCCTG")
wmScore(wm, "CAGGTCGTGGAG")
}
