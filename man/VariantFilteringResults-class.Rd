\name{VariantFilteringResults-class}
\docType{class}

% Classes:
\alias{class:VariantFilteringResults}
\alias{VariantFilteringResults}
\alias{VariantFilteringResults-class}

\alias{class:VariantFilteringResultsAIM}
\alias{VariantFilteringResultsAIM}
\alias{VariantFilteringResultsAIM-class}

\alias{class:VariantFilteringResultsUI}
\alias{VariantFilteringResultsUI}
\alias{VariantFilteringResultsUI-class}

% Methods:
\alias{show,CompressedVRangesList-method}
\alias{show,VariantFilteringResults-method}
\alias{show,VariantFilteringResultsAIM-method}
\alias{show,VariantFilteringResultsUI-method}
\alias{param}
\alias{param,VariantFilteringResults-method}
\alias{inheritanceModel}
\alias{inheritanceModel,VariantFilteringResults-method}
\alias{samples}
\alias{samples<-}
\alias{samples,VariantFilteringResults-method}
\alias{samples<-,VariantFilteringResults-method}
\alias{samples<-,VariantFilteringResults,character-method}
\alias{resetSamples}
\alias{resetSamples,VariantFilteringResults-method}
\alias{bamFiles}
\alias{bamFiles<-}
\alias{bamFiles,VariantFilteringResults-method}
\alias{bamFiles<-,VariantFilteringResults-method}
\alias{bamFiles<-,VariantFilteringResults,BamViews-method}
\alias{selectIndividual}
\alias{selectIndividual<-}
\alias{selectIndividual,VariantFilteringResultsAIM-method}
\alias{selectIndividual<-,VariantFilteringResultsAIM-method}
\alias{selectIndividual<-,VariantFilteringResultsAIM,character-method}
\alias{selectIndividual,VariantFilteringResultsUI-method}
\alias{selectIndividual<-,VariantFilteringResultsUI-method}
\alias{selectIndividual<-,VariantFilteringResultsUI,character-method}
\alias{dbSNPpresent}
\alias{dbSNPpresent<-}
\alias{dbSNPpresent,VariantFilteringResults-method}
\alias{dbSNPpresent<-,VariantFilteringResults-method}
\alias{dbSNPpresent<-,VariantFilteringResults,ANY-method}
\alias{variantType}
\alias{variantType<-}
\alias{variantType,VariantFilteringResults-method}
\alias{variantType<-,VariantFilteringResults,character-method}
\alias{variantLocation}
\alias{variantLocation<-}
\alias{variantLocation,VariantFilteringResults-method}
\alias{variantLocation<-,VariantFilteringResults-method}
\alias{variantLocation<-,VariantFilteringResults,logical-method}
\alias{variantLocation<-,VariantFilteringResults,ANY,logical-method}
\alias{variantConsequence}
\alias{variantConsequence<-}
\alias{variantConsequence,VariantFilteringResults-method}
\alias{variantConsequence<-,VariantFilteringResults-method}
\alias{variantConsequence<-,VariantFilteringResults,logical-method}
\alias{variantConsequence<-,VariantFilteringResults,ANY,logical-method}
\alias{aaChangeType}
\alias{aaChangeType<-}
\alias{aaChangeType,VariantFilteringResults-method}
\alias{aaChangeType<-,VariantFilteringResults,character-method}
\alias{OMIMpresent}
\alias{OMIMpresent<-}
\alias{OMIMpresent,VariantFilteringResults-method}
\alias{OMIMpresent<-,VariantFilteringResults-method}
\alias{OMIMpresent<-,VariantFilteringResults,character-method}
\alias{naMAF}
\alias{naMAF<-}
\alias{naMAF,VariantFilteringResults-method}
\alias{naMAF<-,VariantFilteringResults,logical-method}
\alias{maxMAF}
\alias{maxMAF<-}
\alias{maxMAF,VariantFilteringResults-method}
\alias{maxMAF<-,VariantFilteringResults,numeric-method}
\alias{MAFpop,VariantFilteringResults-method}
\alias{MAFpop}
\alias{MAFpop<-}
\alias{MAFpop<-,VariantFilteringResults,logical-method}
\alias{MAFpop<-,VariantFilteringResults-method}
\alias{MAFpop<-,VariantFilteringResults,ANY,logical-method}
\alias{minPhastCons}
\alias{minPhastCons<-}
\alias{minPhastCons,VariantFilteringResults-method}
\alias{minPhastCons<-,VariantFilteringResults-method}
\alias{minPhastCons<-,VariantFilteringResults,ANY-method}
\alias{minPhylostratum}
\alias{minPhylostratum<-}
\alias{minPhylostratum,VariantFilteringResults-method}
\alias{minPhylostratum<-,VariantFilteringResults-method}
\alias{minPhylostratum<-,VariantFilteringResults,ANY-method}
\alias{minCRYP5ss}
\alias{minCRYP5ss<-}
\alias{minCRYP5ss,VariantFilteringResults-method}
\alias{minCRYP5ss<-,VariantFilteringResults,ANY-method}
\alias{minCRYP5ss<-,VariantFilteringResults-method}
\alias{minCRYP3ss}
\alias{minCRYP3ss<-}
\alias{minCRYP3ss,VariantFilteringResults-method}
\alias{minCRYP3ss<-,VariantFilteringResults,ANY-method}
\alias{minCRYP3ss<-,VariantFilteringResults-method}
\alias{minCUFC}
\alias{minCUFC<-}
\alias{minCUFC,VariantFilteringResults-method}
\alias{minCUFC<-,VariantFilteringResults,numeric-method}
\alias{allVariants,VariantFilteringResults-method}
\alias{allVariants}
\alias{filteredVariants,VariantFilteringResults-method}
\alias{filteredVariants,VariantFilteringResultsAIM-method}
\alias{filteredVariants,VariantFilteringResultsUI-method}
\alias{filteredVariants}
\alias{reportVariants,VariantFilteringResults-method}
\alias{reportVariants,VariantFilteringResultsAIM-method}
\alias{reportVariants,VariantFilteringResultsUI-method}
\alias{reportVariants}
\alias{plot,VariantFilteringResults,ANY-method}

\title{VariantFiltering results class}

\description{
The class \code{VariantFilteringResults} is defined to ease manipulating the results from the calls to functions that filter input genetic variants according to a desired segregating inheritance model (xLinked(), recessiveHomozygous(), etc).

The classes \code{VariantFilteringResultsAIM} and
\code{VariantFilteringResultsUI} are both inherited from the
\code{VariantFilteringResults} one but created through the
\code{allInheritanceModels()} and \code{unrelatedIndividuals()}
functions respectively. 

}
\details{
The class \code{VariantFilteringResults} serves as a purpose of manipulating the results and applying filters on functional annotations of the variants. This manipulation takes place by means of accessor functions enumerated below.
}
\section{Variant Annotation data}{
An \code{VariantFilteringResults} object contains the following annotation data for each variant:
\describe{
  \item{GENE}{
  Gene name given by HGNC.
  }
  \item{CHR}{
  Chromosome name where the variant is located.
  }
  \item{POS}{
  Nucleotide position of the variant within the chromosome.
  }
  \item{LOCATION}{
  Region where the variant is located (coding, intronic, splice site, promoter...)
  }
  \item{TYPE}{
  Type of variant among single nucleotide variant (SNV), insertion,
  deletion, multinucleotide variant (MNV) and deletion followed by an
  insertion (Delins). These types are determined using functions \code{isSNV()},
  \code{isInsertion()}, \code{isDeletion()}, \code{isSubstitution()} and
  \code{isDelins()} from the \code{\link[VariantAnnotation]{VariantAnnotation}}
  package.
  }
  \item{dbSNP}{
  dbSNP identifier derived by position from the annotation packages given
  in the \code{snpdb} argument of the \code{\link{VariantFilteringParam()}}
  function.
  }
  \item{CONSEQUENCE}{
  Consequence of the variant when located in the coding region (synonymous,
  nonsynonymous, etc.)
  }
  \item{TXNAME}{
  Transcript name extracted from the \code{TxDb} annotation package given
  by the \code{txdb} argument of the \code{\link{VariantFilteringParam()}}
  function.
  }
  \item{DESC}{
  HGVS description of the variant.
  }
  \item{PROT}{
  Reference and alternative amino acid in the amino acid chain.
  }
  \item{PolyPhen2}{
  PolyPhen2 prediction for the variant, if corresponding annotation package
  is provided to the \code{otherAnnotations} argument of the
  \code{\link{VariantFilteringParam()}} function.
  }
  \item{SIFT}{
  SIFT/PROVEAN prediction for the variant, if corresponding annotation package
  is provided to the \code{otherAnnotations} argument of the
  \code{\link{VariantFilteringParam()}} function.
  }
  \item{OMIM}{
  OMIM identifier of the gene associated to the variant.
  }
  \item{(POP)_AFKG}{
  Minor Allele Frequencies (MAF) from the 1000 Genomes project, where \code{(POP)}
  is the acronym for the population from where the MAF values are derived.
  When \code{(POP)} is empty, then the MAF value refers to the global MAF for
  that annotation source.
  }
  \item{(POP)_AFESP}{
  Minor Allele Frequencies (MAF) from the NHLBI Exome Sequencing project, where \code{(POP)}
  is the acronym for the population from where the MAF values are derived.
  When \code{(POP)} is empty, then the MAF value refers to the global MAF for
  that annotation source.
  }
  \item{CRYP5ssREF}{
  Score for the cryptic 5'ss for the REF allele respect to the ALT allele.
  }
  \item{CRYP5ssALT}{
  Maximum score for a potential cryptic 5'ss created by the ALT allele.
  }
  \item{CRYP5ssPOS}{
  Position of the allele respect to the position of the dinucleotide \code{GT},
  considering those as positions 1 and 2.
  }
  \item{CRYP3ssREF}{
  Score for the cryptic 3'ss for the REF allele respect to the ALT allele.
  }
  \item{CRYP3ssALT}{
  Maximum score for a potential cryptic 3'ss created by the ALT allele.
  }
  \item{CRYP3ssPOS}{
  Position of the allele respect to the position of the dinucleotide \code{AG},
  considering those as positions 1 and 2.
  }
}
}
\section{Methods}{
In the calls below, \code{x} is a \code{VariantFilteringResults} object.
\describe{
  \item{}{
    \code{param(x)}: return the \code{VariantFilteringParam} input parameter object employed in the call that produced the \code{VariantFilteringResults} object \code{x}.
  }
  \item{}{
    \code{inheritanceModel(x)}: return the model of inheritance employed in the call that produced the \code{VariantFilteringResults} object \code{x}.
  }
  \item{}{
    \code{dbSNPpresent(x)}: flag whether to filter variants present or absent from dbSNP (\code{NA} -do not filter-, \code{"Yes"}, \code{"No"}).
  }
  \item{}{
    \code{variantType(x)}: filter by type of variant (\code{"Any"}, \code{"SNV"}, \code{"Indel"}, \code{"MNV"}).
  }
  \item{}{
    \code{variantLocation(x)}: filter by variant location (\code{"coding"}, \code{"intron"}, \code{"threeUTR"}, \code{"fiveUTR"}, \code{"intergenic"}, \code{"spliceSite"}, \code{"promoter"}).
  }
  \item{}{
    \code{variantConsequence(x)}: filter by variant consequence (\code{"snynonymous"}, \code{"nonsynonymous"}, \code{"frameshift"}, \code{"nonsense"}, \code{"not translated"}).
  }
  \item{}{
    \code{aaChangeType(x)}: filter by type of change of amino acid (\code{"Any"}, \code{"Radical"}, \code{"Conservative"}).
  }
  \item{}{
    \code{OMIMpresent(x)}: flag whether to filter variants whose associated genes are present or absent from OMIM (\code{NA} -do not filter-, \code{"Yes"}, \code{"No"}).
  }
  \item{}{
    \code{naMAF(x)}: flag whether NA maximum MAF values should be included in the filtered variants.
  }
  \item{}{
    \code{maxMAF(x)}: maximum MAF value that a variant may meet among the selected populations.
  }
  \item{}{
    \code{minPhastCons(x)}: minimum phastCons score for nucleotide conservation (\code{NA} -do not filter-, [0-1]).
  }
  \item{}{
    \code{minPhylostratum(x)}: minimum phylostratum for gene conservation (\code{NA} -do not filter-, [1-20]).
  }
  \item{}{
    \code{MAFpop(x)}: selection of populations to use when filtering by maximum MAF value.
  }
  \item{}{
    \code{minCRYP5ss(x)}: minimum weight matrix score on a cryptic 5'ss. \code{NA} indicates this filter is not applied.
  }
  \item{}{
    \code{minCRYP3ss(x)}: minimum weight matrix score on a cryptic 3'ss. \code{NA} indicates this filter is not applied.
  }
  \item{}{
    \code{minCUFC(x)}: minimum absolute codon-usage log2 fold-change.
  }
  \item{}{
    \code{allVariants(x)}: get a \code{GRanges} object with all variants without applying any filter.
  }
  \item{}{
    \code{filteredVariants(x)}: get a \code{GRanges} object with the variants obtained after applying all the filters.
  }
  \item{}{
    \code{selectIndividual(x)}: selection of individuals for further
    analysis. \code{NA} indicates that all individuals present in the
    VCF file are selected.  
  }
  \item{}{
    \code{reportVariants(x, type=c("shiny", "csv", "tsv"), file=NULL)}:
    Builds a report from an \code{VariantFilteringResult} object. Using the
    \code{type} argument, the report can take the form of a flat file in CSV
    or TSV format or a web \code{shiny} app (default) that enables applying
    functional annotation filters in an interactive manner.

    When the \code{shiny} app is closed this method returns a
    \code{VariantFilteringResult} object with the corresponding filters
    switched on or off according to how the app has been interactively used.
  }
  \item{}{
    \code{plot(x, what, flank, ...)}: Method to plot variants using the
    \code{Gviz} package. The argument \code{what} can be either a character
    vector specifying gene or variant identifiers or a chromosome name, or
    a \code{GRanges} object specifying a genomic region. The argument
    \code{flank} is a number of nucleotides to extend the plotting region
    derived from the argument \code{what}. Further arguments in \code{...}
    are passed to the \code{Gviz} function \code{plotTracks()}. See the
    vignette of \code{Gviz} for what can these extra arguments be.
  }
 }
}

\author{D.M. Elurbe and R. Castelo}

\examples{
\dontrun{
library(VariantFiltering)

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"),
                    "CEUtrio.vcf.gz")
CEUped <- file.path(system.file("extdata", package="VariantFiltering"),
                    "CEUtrio.ped")
param <- VariantFilteringParam(vcfFileNames=CEUvcf, pedFileName=CEUped)
reHo <- autosomalRecessiveHomozygous(param)
naMAF(reHo) <- FALSE
maxMAF(reHo) <- 0.05
reHo
head(filteredVariants(reHo))
reportVariants(reHo, type="csv", file="reHo.csv")
}
}

\keyword{classes,methods}
