\name{MafDb2-class}
\docType{class}
\alias{MafDb2-class}
\alias{MafDb2}

\alias{annotateVariants,MafDb2-method}
\alias{mafByOverlaps,MafDb2-method}
\alias{mafByOverlaps}
\alias{populations,MafDb2-method}
\alias{populations}
%%\alias{keytypes,MafDb2-method}
%%\alias{keytypes}
%%\alias{keys,MafDb2-method}
%%\alias{keys}
%%\alias{columns,MafDb2-method}
%%\alias{columns}
%%\alias{select,MafDb2-method}
%%\alias{select}
%%\alias{codeAF2RAW}
%%\alias{decodeRAW2AF}
\alias{organism,MafDb2-method}
\alias{provider,MafDb2-method}
\alias{providerVersion,MafDb2-method}
\alias{referenceGenome,MafDb2-method}
\alias{seqinfo,MafDb2-method}
\alias{seqlengths,MafDb2-method}
\alias{seqnames,MafDb2-method}
\alias{seqlevelsStyle,MafDb2-method}
\alias{show,MafDb2-method}

%%\S4method{keytypes}{MafDb2}(x)
%%\S4method{keys}{MafDb2}(x, keytype)
%%\S4method{columns}{MafDb2}(x)
%%\S4method{select}{MafDb2}(x, keys, columns, keytype)
%%\S4method{annotateVariants}{MafDb2}(annObj, variantsVR, param, BPPARAM=bpparam("SerialParam"))

\title{MafDb2 class}

\description{
Class for annotation packages storing minor allele frequency data.
}
\usage{
\S4method{mafByOverlaps}{MafDb2}(x, ranges, pop, caching)
\S4method{populations}{MafDb2}(x)
}
\arguments{
  \item{x}{A \code{MafDb2} object.}
  \item{ranges}{Either a \code{GRanges} object or a character string vector with the format \code{"CHR:START[-END]"}.}
  \item{pop}{A character string with the populations for which we want to retrieve MAF values.}
  \item{caching}{logical; TRUE (default) indicates that the function stores into main memeory the MAF data as it gets loaded from disk, improving performance; FALSE forces this function to load MAF data from disk each time, decreasing performance and memory requirements.}
%%  \item{varID}{A variant identifier, typically a \code{rsxxxx} dbSNP identifier.}
%%  \item{keytype}{the keytype that matches the keys used. For \code{MafDb2} objects
%%    there is at the moment only one type of key which is the variant identifier
%%    provided by the original data manufacturer.}
%%  \item{keys}{the keys to select records from the database. All possible keys
%%    are turned by using the \code{keys} method.}
%%  \item{columns}{the columns or kinds of things that can be retrieved from the
%%    database. As with \code{keys}, all possible columns are returned by using the
%%    \code{columns} method.}
%%  \item{annObj}{A \code{MafDb2} object.}
%%  \item{variantsVR}{A \code{VRanges} object with the variants to annotate.}
%%  \item{param}{A \code{\link{VariantFilteringParam-class}} object.}
%%  \item{BPPARAM}{An object of class \code{BiocParallelParam} specifying parameters related to
%%    the parallel execution of this function. See function \code{\link[BiocParallel]{bpparam}()}
%%    from the \code{BiocParallel} package.}
}
\details{
The \code{MafDb2} class and associated methods serve the purpose of creating
annotation packages that store minor allele frequency data. Two such
annotation packages are:

\tabular{ll}{
  \code{MafDb2.ALL.wgs.phase1.release.v3.20101123} \tab MAF values from the 1000 Genomes project
        downloaded in October 2013 from
        \url{ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz}. \cr

  \code{MafDb2.ESP6500SI.V2.SSA137} \tab MAF values from 6500 ESP exomes downloaded in March 2015 from
        \url{http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz}. \cr
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
%%\seealso{
%%  \code{\link[MafDb2.ExAC.r0.3.1.snvs.hs37d5]{MafDb2.ExAC.r0.3.1.snvs.hs37d5}}
%%}
\examples{

  ## lookup allele frequencies for rs1129038, a SNP associated to blue and brown eye colors
  ## as reported by Eiberg et al. Blue eye color in humans may be caused by a perfectly associated
  ## founder mutation in a regulatory element located within the HERC2 gene inhibiting OCA2 expression.
  ## Human Genetics, 123(2):177-87, 2008 [http://www.ncbi.nlm.nih.gov/pubmed/18172690]

\dontrun{

  if (require(MafDb.ExAC.r0.3.1.snvs.hs37d5)) {
    MafDb.ExAC.r0.3.1.snvs.hs37d5

    ## specialized interface
    populations(MafDb.ExAC.r0.3.1.snvs.hs37d5)
  }

}

}
\source{
  \url{http://exac.broadinstitute.org}
}
\keyword{datasets}
