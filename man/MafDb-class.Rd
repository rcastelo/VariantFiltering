\name{MafDb-class}
\docType{class}
\alias{MafDb-class}
\alias{MafDb}

\alias{annotateVariants}
\alias{annotateVariants,MafDb-method}
\alias{annotateVariants,OrgDb-method}
\alias{annotateVariants,PolyPhenDb-method}
\alias{annotateVariants,PROVEANDb-method}
\alias{annotateVariants,SNPlocs-method}
\alias{annotateVariants,XtraSNPlocs-method}
\alias{annotateVariants,TxDb-method}
\alias{snpid2maf,MafDb-method}
\alias{snpid2maf}
\alias{knownVariantsMAFcols,MafDb-method}
\alias{knownVariantsMAFcols}
\alias{keytypes,MafDb-method}
\alias{keytypes}
\alias{keys,MafDb-method}
\alias{keys}
\alias{columns,MafDb-method}
\alias{columns}
\alias{select,MafDb-method}
\alias{select}
\alias{codeAF2RAW}
\alias{decodeRAW2AF}

\title{MafDb class}

\description{
Class for annotation packages storing minor allele frequency data.
}
\usage{
\S4method{snpid2maf}{MafDb}(mafdb, varID)
\S4method{knownVariantsMAFcols}{MafDb}(mafdb)
\S4method{keytypes}{MafDb}(x)
\S4method{keys}{MafDb}(x, keytype)
\S4method{columns}{MafDb}(x)
\S4method{select}{MafDb}(x, keys, columns, keytype)
\S4method{annotateVariants}{MafDb}(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam"))
}
\arguments{
  \item{mafdb}{A \code{MafDb} object.}
  \item{x}{A \code{MafDb} object.}
  \item{varID}{A variant identifier, typically a \code{rsxxxx} dbSNP identifier.}
  \item{keytype}{the keytype that matches the keys used. For \code{MafDb} objects
    there is at the moment only one type of key which is the variant identifier
    provided by the original data manufacturer.}
  \item{keys}{the keys to select records from the database. All possible keys
    are turned by using the \code{keys} method.}
  \item{columns}{the columns or kinds of things that can be retrieved from the
    database. As with \code{keys}, all possible columns are returned by using the
    \code{columns} method.}
  \item{annObj}{A \code{MafDb} object.}
  \item{variantsVR}{A \code{VRanges} object with the variants to annotate.}
  \item{param}{A \code{\link{VariantFilteringParam-class}} object.}
  \item{BPPARAM}{An object of class \code{BiocParallelParam} specifying parameters related to
    the parallel execution of this function. See function \code{\link[BiocParallel]{bpparam}()}
    from the \code{BiocParallel} package.}
}
\details{
The \code{MafDb} class and associated methods serve the purpose of creating
annotation packages that store minor allele frequency data. Two such
annotation packages are:

\tabular{ll}{
  \code{MafDb.1Kgenomes.phase1.hs37d5} \tab MAF values from the 1000 Genomes Project Phase 1. \cr

  \code{MafDb.1Kgenomes.phase3.hs37d5} \tab MAF values from the 1000 Genomes Project Phase 3. \cr
}

This object class tries to reduce the disk space required to store allele frequencies (AFs)
for millions of SNPs by coding AF float values, which range between 0 and 1, into a single-byte
\code{raw} object type. To achieve this, the original AF values are rounded and coded as follows:

\itemize{
  \item AF >= 0.01 & AF <= 1 values are rounded to 2 digits, where values 0.01, ..., 0.99, 1, are coded as
        raw byte values 1 to 100.
  \item AF >= 0.001 & AF < 0.01 values are rounded to 3 digits, where values 0.001, ..., 0.009 are coded
        as raw byte values 101 to 109.
  \item AF >= 0.0001 & AF < 0.001 are rounded to 4 digits, where values 0.0001, ..., 0.0009 are coded
        as raw byte values 111 to 119.
  \item AF >= 0.00001 & AF < 0.0001 are rounded to 5 digits, where values 0.00001, ..., 0.00009 are coded
        as raw byte values 121 to 129.
  \item AF < 0.00001 are rounded to 6 digits, where values 0, 0.000001, ..., 0.000009 are coded
        as raw byte values 130 to 139.
  \item AF \code{NA} values are coded to raw byte value of 255. Note that by default \code{NA}
        values are coded by the raw byte value 0 but this corresponds by default to the null
        string when raw byte values are coerced into char, which is problematic when storing
        this data as \code{CHAR} values in a SQLite database. This precludes using this original coding
        of \code{NA} values.
}

A further compression of these data is performed in the cases of variants with mutiple
alternative alleles. In those cases, instead of storing the AF of each alternate allele only
the maximum AF value is stored.
}
\author{R. Castelo}
\seealso{
  \code{\link[MafDb.1Kgenomes.phase1.hs37d5]{MafDb.1Kgenomes.phase1.hs37d5}}
  \code{\link[MafDb.1Kgenomes.phase3.hs37d5]{MafDb.1Kgenomes.phase3.hs37d5}}
}
\examples{

  ## lookup allele frequencies for rs1129038, a SNP associated to blue and brown eye colors
  ## as reported by Eiberg et al. Blue eye color in humans may be caused by a perfectly associated
  ## founder mutation in a regulatory element located within the HERC2 gene inhibiting OCA2 expression.
  ## Human Genetics, 123(2):177-87, 2008 [http://www.ncbi.nlm.nih.gov/pubmed/18172690]

  if (require(MafDb.1Kgenomes.phase3.hs37d5)) {
    mafdb <- MafDb.1Kgenomes.phase3.hs37d5
    mafdb

    ## specialized interface
    knownVariantsMAFcols(mafdb)
    snpid2maf(mafdb, "rs1129038")

    ## standard AnnotationDbi interface
    keytypes(mafdb)
    columns(mafdb)
    select(mafdb, keys="rs1129038", columns=c("varID", "AF"))
  }

}
\source{
  \url{ftp://ftp.1000genomes.ebi.ac.uk}
}
\keyword{datasets}
