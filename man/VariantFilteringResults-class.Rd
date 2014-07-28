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
\alias{show,VariantFilteringResults-method}
\alias{show,VariantFilteringResultsAIM-method}
\alias{show,VariantFilteringResultsUI-method}
\alias{param}
\alias{param,VariantFilteringResults-method}
\alias{inheritanceModel}
\alias{inheritanceModel,VariantFilteringResults-method}
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
  Gene name by the HGNC
  }
  \item{CHR}{
  Chromosome number, following the UCSC standards
  }
  \item{POS}{
  Location of the variant inside de chromosome
  }
  \item{LOCATION}{
  Region where the variant is located (coding, splice site, promoter...)
  }
  \item{TYPE}{
  Type of variant (SNV, InDel...)
  }
  \item{dbSNP}{
  dbSNP ID
  }
  \item{CONSEQUENCE}{
  Consequence in the translation of the sequence
  }
  \item{TXNAME}{
  Transcript name extracted from the \code{TxDb} object defined at \code{txdb}
  }
  \item{CDS}{
  Reference and variant nucleotides in the cDNA
  }
  \item{PROT}{
  Reference and variant amino acid in the amino acid chain
  }
  \item{PolyPhen2}{
  PolyPhen2 prediction for the variant
  }
  \item{SIFT}{
  SIFT prediction for the variant
  }
  \item{OMIM}{
  OMIM ID entry for the gene
  }
  \item{AFKG}{
  Minor Allele Frequency from the 1000 Genomes project, all populations
  }
  \item{AMR_AFKG}{
  Minor Allele Frequency from the 1000 Genomes project, Ad Mixed American
  }
  \item{ASN_AFKG}{
  Minor Allele Frequency from the 1000 Genomes project, East Asian
  }
  \item{AFR_AFKG}{
  Minor Allele Frequency from the 1000 Genomes project, African
  }  
  \item{EUR_AFKG}{
  Minor Allele Frequency from the 1000 Genomes project, European
  }    
  \item{AFESP}{
  Minor Allele Frequency from the NHLBI Exome Sequencing Project, all populations
  }
  \item{EA_AFESP}{
  Minor Allele Frequency from the NHLBI Exome Sequencing Project, European American
  }
  \item{AA_AFESP}{
  Minor Allele Frequency from the NHLBI Exome Sequencing Project, African American
  }
  \item{CRYP5ssREF}{
  Score for the cryptic 5'ss for the REF allele respect to the ALT allele
  }
  \item{CRYP5ssALT}{
  Maximum score for the cryptic 5'ss taking into account all the possible positions within the window
  }
  \item{CRYP5ssPOS}{
  Position of the allele respect to the position of the dinucleotide GT, considering those as positions 1 and 2
  }
  \item{CRYP3ssREF}{
  Score for the cryptic 3'ss for the REF allele respect to the ALT allele
  }
  \item{CRYP3ssALT}{
  Maximum score for the cryptic 3'ss taking into account all the possible positions within the window
  }
  \item{CRYP3ssPOS}{
  Position of the allele respect to the position of the dinucleotide AG, considering those as positions 1 and 2
  }
}
}
\section{Accessor methods}{
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
