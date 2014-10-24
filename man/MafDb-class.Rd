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
\alias{annotateVariants,TxDb-method}
\alias{fetchKnownVariantsByID,MafDb-method}
\alias{fetchKnownVariantsByID}
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
Class for annotation packages storing minimum allele frequency data.
}
\usage{
\S4method{fetchKnownVariantsByID}{MafDb}(mafdb, varID)
\S4method{knownVariantsMAFcols}{MafDb}(mafdb)
\S4method{keytypes}{MafDb}(x)
\S4method{keys}{MafDb}(x, keytype)
\S4method{columns}{MafDb}(x)
\S4method{select}{MafDb}(x, keys, columns, keytype)
\S4method{annotateVariants}{MafDb}(annObj, variantsGR)
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
  \item{variantsGR}{A \code{GRanges} object with the variants to annotate.}
}
\details{
The \code{MafDb} class and associated methods serve the purpose of creating
annotation packages that store minimum allele frequency data. Two such
annotation packages are:

\tabular{ll}{
  \code{MafDb.ALL.wgs.phase1.release.v3.20101123} \tab MAF values from the 1000 Genomes project
        downloaded in November 2013 from
        \url{ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz}. \cr

  \code{MafDb.ESP6500SI.V2.SSA137.dbSNP138} \tab MAF values from 6500 ESP exomes downloaded in November 2013 from
        \url{http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf.tar.gz}. \cr
}

This object class tries to reduce the disk space required to store allele frequencies (AFs)
for millions of SNPs by coding AF float values, which range between 0 and 1, into a single-byte
\code{raw} object type. To achieve this, the original AF values are rounded and coded as follows:

\itemize{
  \item AF < 0.001 values are rounded to 4 digits, where values 0, 0.0001, ..., 0.0009 are coded
        as raw byte values 1 to 10.
  \item AF < 0.01 values are rouneded to 3 digits, where values 0.001, ... 0.009, are coded as
        raw byte values 11 to 19.
  \item AF >= 0.01 values are rounded to 2 digits and coded as raw byte values 20 to 119.
  \item AF \code{NA} values are coded to raw byte value of 255. Note that by default \code{NA}
        values are coded by the raw byte value 0 but this corresponds by default to the null
        string when raw byte values are coerced into char. This precludes using this original coding
        of \code{NA} values.
}

A further compression of these data is performed in the cases of variants with mutiple
alternative alleles. In those cases, instead of storing the AF of each alternate allele only
the maximum AF value is stored and as alternate allele the string concatenation of all alleles,
separated by a slash character, is stored.
}
\value{.}
\author{R. Castelo}
\seealso{
  \code{\link[MafDb.ALL.wgs.phase1.release.v3.20101123]{MafDb.ALL.wgs.phase1.release.v3.20101123}}
  \code{\link[MafDb.ESP6500SI.V2.SSA137.dbSNP138]{MafDb.ESP6500SI.V2.SSA137.dbSNP138}}
}
\examples{
  if (require(MafDb.ESP6500SI.V2.SSA137.dbSNP138)) {
    MafDb.ESP6500SI.V2.SSA137.dbSNP138

    ## specialized interface
    knownVariantsMAFcols(MafDb.ESP6500SI.V2.SSA137.dbSNP138)
    fetchKnownVariantsByID(MafDb.ESP6500SI.V2.SSA137.dbSNP138, "rs199529001")

    ## standard AnnotationDbi interface
    keytypes(MafDb.ESP6500SI.V2.SSA137.dbSNP138)
    columns(MafDb.ESP6500SI.V2.SSA137.dbSNP138)
    select(MafDb.ESP6500SI.V2.SSA137.dbSNP138,
           keys="rs199529001", columns=c("varID", "chrom", "AF"))
  }
}
\source{
  \url{ftp://ftp.1000genomes.ebi.ac.uk}

  \url{http://evs.gs.washington.edu/EVS}
}
\keyword{datasets}
