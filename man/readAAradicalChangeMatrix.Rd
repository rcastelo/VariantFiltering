\name{readAAradicalChangeMatrix}
\alias{readAAradicalChangeMatrix}

\title{
Read matrix of amino acid radical changes
}
\description{
Function to read and parse a tab-separated file of amino acid properties into a matrix of logical values indicating whether the change of one amino acid by another can be considered radical or conservative according to the chemical properties specified in the input file.
}

\usage{
readAAradicalChangeMatrix(file)
}
\arguments{
  \item{file}{
    A file containing a classification of amino acids with respect to one or more chemical properties. Its particular format
    should match the one from the file called \code{AA_chemical_properties_HanadaGojoboriLi2006.tsv} found in the
    \code{extdata} folder of this package. This file is based on Table 1 from Hanada et al. (2006).
  }
}
\details{
The input file should contain one or more columns each of them forming a logical mask specifing sets of amino acids
sharing some chemical property.
}
\value{
  An squared symmetric matrix with as many rows and columns as amino acids and whose cells contain
  logical values. These values are set to TRUE whenever the amino acid change implied by row and column
  is considered radical and FALSE when considered conservative. Amino acid changes
  within a chemical property are defined as conservative and radical otherwise.
}
\references{
  Hanada, K., Gojobori, T. and Li, W. Radical amino acid change versus positive selection in the evolution of viral envelope proteins.
  Gene, 385:83-88, 2006.
}
\author{
R. Castelo
}
\seealso{
\code{\link{VariantFilteringParam}}
}
\examples{
aamat <- readAAradicalChangeMatrix(file.path(system.file("extdata", package="VariantFiltering"),
                                   "AA_chemical_properties_HanadaGojoboriLi2006.tsv"))
aamat[1:5, 1:5]
}
\keyword{utilities}
