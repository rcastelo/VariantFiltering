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
\alias{annotateVariants,TranscriptDb-method}
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
The \code{MafDb} class and associated methods serve the purpose of creating annotation packages that store
minimum allele frequency data. Two such annotation packages are:

\tabular{ll}{
  \code{MafDb.ALL.wgs.phase1.release.v3.20101123} \tab MAF values from the 1000 Genomes project
        downloaded in November 2013 from
        \url{ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz}. \cr

  \code{MafDb.ESP6500SI.V2.SSA137.dbSNP138} \tab MAF values from 6500 ESP exomes downloaded in November 2013 from
        \url{http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf.tar.gz}. \cr
}

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
