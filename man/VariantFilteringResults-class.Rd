\name{VariantFilteringResults-class}
\docType{class}

% Class
\alias{class:VariantFilteringResults}
\alias{VariantFilteringResults}
\alias{VariantFilteringResults-class}
\alias{class:CutoffsList}
\alias{CutoffsList}
\alias{CutoffsList-class}

% Methods:
\alias{show,CompressedVRangesList-method}
\alias{show,VariantFilteringResults-method}
\alias{show,VariantFilteringResultsUI-method}
\alias{summary,VariantFilteringResults-method}
\alias{param}
\alias{param,VariantFilteringResults-method}
\alias{length}
\alias{length,VariantFilteringResults-method}
\alias{inheritanceModel}
\alias{inheritanceModel,VariantFilteringResults-method}
\alias{sog}
\alias{sog,VariantFilteringResults-method}
\alias{filters}
\alias{filters<-}
\alias{filters,VariantFilteringResults-method}
\alias{filters<-,VariantFilteringResults-method}
\alias{filtersMetadata}
\alias{filtersMetadata,VariantFilteringParam-method}
\alias{filtersMetadata,VariantFilteringResults-method}
\alias{cutoffs}
\alias{cutoffs<-}
\alias{cutoffs,VariantFilteringResults-method}
\alias{cutoffs<-,VariantFilteringResults-method}
\alias{cutoffs<-,VariantFilteringResults,CutoffsList-method}
\alias{cutoffs<-,VariantFilteringResults,logical-method}
\alias{change}
\alias{change<-}
\alias{change<-,CutoffsList}
\alias{change<-,CutoffsList,ANY,integer-method}
\alias{change<-,CutoffsList,ANY,numeric-method}
\alias{change<-,CutoffsList,ANY,character-method}
\alias{change<-,CutoffsList,ANY,logical-method}
\alias{sortings}
\alias{sortings<-}
\alias{sortings,VariantFilteringResults-method}
\alias{sortings<-,VariantFilteringResults,CutoffsList-method}
\alias{softFilterMatrix,VariantFilteringResults-method}
\alias{softFilterMatrix<-,VariantFilteringResults-method}
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
\alias{annoGroups}
\alias{annoGroups,VariantFilteringResults-method}
\alias{referenceGenome,VariantFilteringResults-method}
\alias{allVariants,VariantFilteringResults-method}
\alias{allVariants}
\alias{filteredVariants,VariantFilteringResults-method}
\alias{filteredVariants}
\alias{reportVariants,VariantFilteringResults-method}
\alias{reportVariants}
\alias{browseVariants}
\alias{plot,VariantFilteringResults,ANY-method}

\title{VariantFilteringResults class}

\description{
The \code{VariantFilteringResults} class is used to store the kind of object obtained as a result of an analysis using the functions \code{\link{unrelatedIndividuals}()}, \code{\link{autosomalRecessiveHomozygous}()}, \code{\link{autosomalRecessiveHeterozygous}()}, \code{\link{autosomalDominant}()}, \code{\link{deNovo}()} and \code{\link{xLinked}()}. Its purpose is to ease the task of filtering and prioritizing the variants annotated by those functions.
}
\details{
Variants are stored within a \code{VariantFilteringResults} object using a \code{\link[VariantAnnotation]{VRanges}} object, which also holds the variant annotations in its metadata columns. VariantFiltering adds the following core set of annotations.
\describe{
  \item{LOCATION}{
  Region where the variant is located (coding, intronic, splice site, promoter, ...) as given by the function \code{\link[VariantAnnotation]{locateVariants}()} from the \code{VariantAnnotation} package.
  }
  \item{LOCSTART}{
  Start position of the variant within the region defined by the \code{LOCATION} annotation.
  }
  \item{GENEID}{
  Gene identifier derived with the transcript-centric annotation package given in the \code{txdb} argument
  of the \code{\link{VariantFilteringParam}()} function, typically an Entrez Gene identifier.
  }
  \item{GENE}{
  Gene name given by HGNC derived with the gene-centric annotation package given in the \code{orgdb} argument
  of the \code{\link{VariantFilteringParam}()} function.
  }
  \item{TYPE}{
  Type of variant, either a single nucleotide variant (SNV), an insertion,
  a deletion, a multinucleotide variant (MNV) or a deletion followed by an
  insertion (Delins). These types are determined using functions
  \code{\link[VariantAnnotation]{isSNV}()},
  \code{\link[VariantAnnotation]{isInsertion}()},
  \code{\link[VariantAnnotation]{isDeletion}()},
  \code{\link[VariantAnnotation]{isSubstitution}()} and
  \code{\link[VariantAnnotation]{isDelins}()} from the
  \code{VariantAnnotation} package.
  }
  \item{dbSNP}{
  dbSNP identifier derived by position from the annotation packages given
  in the \code{snpdb} argument of the \code{\link{VariantFilteringParam}()}
  function.
  }
  \item{cDNALOC}{
  Location of the variant along the processed transcript, when the variant
  belongs to an exonic region.
  }
  \item{CONSEQUENCE}{
  Consequence of the variant when located in the coding region (synonymous,
  nonsynonymous, missense, nonsense o frameshift) as given by the function
  \code{\link[VariantAnnotation]{predictCoding}()} from the \code{VariantAnnotation} package.
  }
  \item{TXNAME}{
  Transcript name extracted from the \code{TxDb} annotation package given
  by the \code{txdb} argument of the \code{\link{VariantFilteringParam}()}
  function.
  }
  \item{HGVSg}{
  HGVS description of the variant at genomic level.
  }
  \item{HGVSc}{
  HGVS description of the variant at coding level.
  }
  \item{HGVSp}{
  HGVS description of the variant at protein level.
  }
  \item{OMIM}{
  OMIM identifier of the gene associated to the variant derived with the gene-centric
  annotation package given in the \code{orgdb} argument
  of the \code{\link{VariantFilteringParam}()} function.
  }
  \item{AAchangeType}{
  In the case of coding variants, whether the amino acid change is conservative or
  radical according to the matrix of amino acid biochemical properties given in the
  argument \code{radicalAAchangeFilename} of the \code{\link{VariantFilteringParam}()} function.
  }
  \item{SCORE5ssREF}{
  Score for the cryptic 5'ss for the REF allele respect to the ALT allele.
  }
  \item{SCORE5ssALT}{
  Maximum score for a potential cryptic 5'ss created by the ALT allele.
  }
  \item{SCORE5ssPOS}{
  Position of the allele respect to the position of the dinucleotide \code{GT},
  considering those as positions 1 and 2.
  }
  \item{SCORE3ssREF}{
  Score for the cryptic 3'ss for the REF allele respect to the ALT allele.
  }
  \item{SCORE3ssALT}{
  Maximum score for a potential cryptic 3'ss created by the ALT allele.
  }
  \item{SCORE3ssPOS}{
  Position of the allele respect to the position of the dinucleotide \code{AG},
  considering those as positions 1 and 2.
  }
}
}
\section{Accessors}{
  A \code{VariantFilteringResults} has the following set of accessor methods.

  \describe{
    \item{\code{length(x)}:}{ total number of variants stored internatlly within the
      \code{VRanges} object. Note that this number will be typically larger than the number
      of variantes in the input VCF object because each of them is copied for each combination
      of alternate allele, annotated region and sample.
    }
    \item{\code{param(x)}:}{ returns the \code{VariantFilteringParam} input parameter
      object employed in the call that produced the \code{VariantFilteringResults} object \code{x}.
    }
    \item{\code{inheritanceModel(x)}:}{ returns the model of inheritance employed in the
      call that produced the \code{VariantFilteringResults} object \code{x}.
    }
    \item{\code{samples(object)}:}{ active samples from which the current filtered variants were derived. If the
      \code{x} was obtained with \code{\link{unrelatedIndividuals}()}, then the replace method
      \code{samples(object)<-} can be used to restrict the subset of active samples. In every other case
      (\code{\link{autosomalDominant}()}, etc. ) active samples cannot be changed.
    }
    \item{\code{resetSamples(object)}:}{ set back as active samples the initial set of samples specified
      in the input parameter object.
    }
    \item{\code{sog(x)}:}{ Sequence Ontology (SO) graph (actually, an acyclic digraph)
      returned as a \code{\link[graph]{graphNEL}} object, whose vertices are SO terms,
      edges represent ontology relationships and vertex attributes \code{vcfIdx} and
      \code{varIdx} contain what variants are annotated to each SO term. These annotations
      can be directly retrieved from the SO graph with the \code{\link[graph]{nodeData}()}
      function from the \code{graph} package. The \code{summary()} function described
      in this manual page allows one to tally the number of variants in each SO term throughout
      the entire SO hierarchy.
    }
    \item{\code{bamFiles(x)}:}{ access and update the \code{\link[Rsamtools]{BamViews}} object containing
      references to BAM files from which the input VCF files were derived. Initially this is empty.
    }
    \item{\code{allVariants(x, groupBy="sample")}:}{ returns a \code{\link[VariantAnnotation]{VRangesList}}
      object with all variants grouped by default by sample. Using the argument \code{groupBy}
      we can specify any metadata column to be used to group variants. If the value given to
      \code{groupBy} does not correspond to any such columns, a
      \code{\link[VariantAnnotation]{VRanges}} object with all variants together is returned.
    }
    \item{\code{filteredVariants(x, groupBy="sample")}:}{ it works like \code{allVariants(x)}
      but instead of returning all variants, it returns only those who pass the active
      filters; see \code{filters()} and \code{cutoffs()} below.
    }
  }
}

\section{Filters and cutoffs}{

  The variants contained in a \code{VariantFilteringResults} object can be filtered using
  the \code{\link[S4Vectors]{FilterRules}} mechanism, defined in the \code{S4Vectors} package,
  by using the functions \code{filters()} and \code{cutoffs()} described below. There are
  additional functions, also described in this section, to facilitate this task on the set
  of core annotations provided by \code{VariantFiltering}.

  \describe{
    \item{\code{filters(x)}:}{ get the current \code{\link[S4Vectors]{FilterRules}} object that defines
      the available set of filter criteria that one can use to filter the variants contained in
      \code{x}. This can also be used as a replacement function \code{filters(x)<-} to update
      this set of filters. The actual filtering is done when calling the function
      \code{filteredVariants()}.
    }
    \item{\code{filtersMetadata(x)}:}{ metadata about the available filters.
    }
    \item{\code{cutoffs(x)}:}{ get cutoffs from the available filters.
    }
    \item{\code{change(x, cutoff)<-}:}{ change cutoffs from the available filters. Here, argument \code{x}
      is a \code{CutoffsList} object given by the method \code{cutoffs()}, and argument \code{cutoff}
      is a character string with the name of the cutoff.
    }
    \item{\code{softFilterMatrix(x)}:}{ get and update the variant by filter matrix; see
      \code{\link[VariantAnnotation]{softFilterMatrix}()} in the \code{VariantAnnotation} package.
    }
    \item{\code{dbSNPpresent(x)}:}{ flag whether to filter variants present or absent from dbSNP (\code{NA} -do not filter-, \code{"Yes"}, \code{"No"}).
    }
    \item{\code{variantType(x)}:}{ filter by type of variant ( \code{"SNV"}, \code{"Insertion"}, \code{"Deletion"}, \code{"MNV"}, \code{"Delins"}).
    }
    \item{\code{variantLocation(x)}:}{ filter by variant location (\code{"coding"}, \code{"intron"}, \code{"threeUTR"}, \code{"fiveUTR"}, \code{"intergenic"}, \code{"spliceSite"}, \code{"promoter"}).
    }
    \item{\code{variantConsequence(x)}:}{ filter by variant consequence (\code{"snynonymous"}, \code{"nonsynonymous"}, \code{"frameshift"}, \code{"nonsense"}, \code{"not translated"}).
    }
    \item{\code{aaChangeType(x)}:}{ filter by type of change of amino acid (\code{"Radical"}, \code{"Conservative"}).
    }
    \item{\code{OMIMpresent(x)}:}{ flag whether to filter variants whose associated genes are present or absent from OMIM (\code{NA} -do not filter-, \code{"Yes"}, \code{"No"}).
    }
    \item{\code{naMAF(x)}:}{ flag whether NA maximum MAF values should be included in the filtered variants.
    }
    \item{\code{maxMAF(x)}:}{ maximum MAF value that a variant may meet among the selected populations.
    }
    \item{\code{minPhastCons(x)}:}{ minimum phastCons score for nucleotide conservation (\code{NA} -do not filter-, [0-1]).
    }
    \item{\code{minPhylostratum(x)}:}{ minimum phylostratum for gene conservation (\code{NA} -do not filter-, [1-20]).
    }
    \item{\code{MAFpop(x)}:}{ selection of populations to use when filtering by maximum MAF value.
    }
    \item{\code{minScore5ss(x)}:}{ minimum weight matrix score on a cryptic 5'ss. \code{NA} indicates this filter is not applied.
    }
    \item{\code{minScore3ss(x)}:}{ minimum weight matrix score on a cryptic 3'ss. \code{NA} indicates this filter is not applied.
    }
    \item{\code{minCUFC(x)}:}{ minimum absolute codon-usage log2 fold-change.
    }
  }
}

\section{Summarization, visualization and reporting}{

  The following functions help in summarizing, visualizing and reporting the
  fiiltered variants.

  \describe{
    \item{\code{summary(object, method=c("SO", "SOfull", "bioc"))}:}{ tally the current
      filtered set of variants to features. By default, features are Sequence
      Ontology (SO) terms to which variants are annotated by \code{VariantFiltering}.
      The \code{method} argument allows the user to change this default setting to
      tallying throughout the entire SO hierarchy. Both options, \code{SO} and
      \code{SOfull} can be used in combination with the cutoff \code{SOterms}; see
      the vignette. The option \code{method="bioc"} considers as features the
      regions and consequences annotated by functions
      \code{\link[VariantAnnotation]{locateVariants}()} and
      \code{\link[VariantAnnotation]{predictCoding}()} from the \code{VariantAnnotation}
      package. The result is returned as a \code{data.frame} object.
    }
    \item{\code{plot(x, what, sampleName, flankingNt=20, showAlnNtCutoff=200, isPaired=FALSE, ...)}:}{
      Plot variants using the \code{Gviz} package. The argument \code{what} can be
      either a character vector specifying gene or variant identifiers or a
      chromosome name, or a \code{GRanges} object specifying a genomic region. The
      argument \code{sampleName} is optional and allows the user to plot the aligned
      reads and coverage from a specific sample, located in the plotted region, when
      the corresponding BAM file has been linked to the object with \code{bamFiles()}.
      The argument \code{flankingNt} is a number of nucleotides to extend the plotting
      region derived from the argument \code{what}. The argument \code{showAlnNtCutoff}
      is the region size cutoff below which it will be attempted to plot the aligned reads.
      The argument \code{isPaired} is passed directly to the \code{Gviz} function
      \code{\link[Gviz]{AlignmentsTrack}()} which streams over the BAM file to plot the reads
      and sets whether the BAM file contains single (default) or paired-end reads.
      Further arguments in \code{...} are passed to the \code{Gviz} function
      \code{\link[Gviz]{plotTracks}()} and can be used to fine-tune the final plot; see
      the vignette of \code{Gviz} to find out what these arguments are.
    }
    \item{\code{reportVariants(x, type=c("shiny", "csv", "tsv"), file=NULL)}:}{
      Builds a report from the \code{VariantFilteringResult} object \code{x}. Using
      the \code{type} argument, the report can take the form of a flat file in CSV
      or TSV format or a web \code{shiny} app (default) that enables applying
      functional annotation filters in an interactive manner.

      When the \code{shiny} app is closed this method returns a
      \code{VariantFilteringResult} object with the corresponding filters
      switched on or off according to how the app has been interactively used.
    }

  }
}

\author{R. Castelo}

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
