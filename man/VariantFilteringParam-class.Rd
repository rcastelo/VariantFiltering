\name{VariantFilteringParam-class}
\docType{class}

\alias{sequence_variant.gSOXP}

%% Classes:
\alias{class:VariantFilteringParam}
\alias{VariantFilteringParam}
\alias{VariantFilteringParam-class}

%% Methods:
\alias{show,VariantFilteringParam-method}
\alias{$,VariantFilteringParam-method}
\alias{names,VariantFilteringParam-method}
\alias{filters,VariantFilteringParam-method}
\alias{cutoffs,VariantFilteringParam-method}
\alias{sog,VariantFilteringParam-method}

%% Functions:
\alias{spliceSiteMatricesHuman}
\alias{variantLocations}

\title{VariantFiltering parameter class}

\description{
The class \code{VariantFilteringParam} is defined to ease configuring the call
to the functions that filter input genetic variants according to a desired
segregating inheritance model (\code{xLinked()}, \code{autosomalRecessiveHomozygous()}, etc).
}
\usage{
VariantFilteringParam(vcfFilenames, pedFilename=NA_character_,
                      bsgenome="BSgenome.Hsapiens.1000genomes.hs37d5",
                      orgdb="org.Hs.eg.db",
                      txdb="TxDb.Hsapiens.UCSC.hg19.knownGene",
                      snpdb="SNPlocs.Hsapiens.dbSNP144.GRCh37",
                      weightMatricesFilenames=NA,
                      weightMatricesLocations=rep(list(variantLocations()), length(weightMatricesFilenames)),
                      weightMatricesStrictLocations=rep(list(FALSE), length(weightMatricesFilenames)),
                      radicalAAchangeFilename=file.path(system.file("extdata",
                                                                    package="VariantFiltering"),
                                                        "AA_chemical_properties_HanadaGojoboriLi2006.tsv"),
                      codonusageFilename=file.path(system.file("extdata",
                                                               package="VariantFiltering"),
                                                   "humanCodonUsage.txt"),
                      geneticCode=getGeneticCode("SGC0"),
                      allTranscripts=FALSE,
                      regionAnnotations=list(CodingVariants(), IntronVariants(),
                                             FiveSpliceSiteVariants(), ThreeSpliceSiteVariants(),
                                             PromoterVariants(), FiveUTRVariants(), ThreeUTRVariants()),
                      otherAnnotations=c("MafDb.1Kgenomes.phase3.hs37d5",
                                         "MafDb.ExAC.r0.3.1.snvs.hs37d5",
                                         "PolyPhen.Hsapiens.dbSNP131",
                                         "SIFT.Hsapiens.dbSNP137",
                                         "phastCons100way.UCSC.hg19",
                                         "humanGenesPhylostrata"),
                      geneKeytype=NA_character_,
                      yieldSize=NA_integer_)
\S4method{show}{VariantFilteringParam}(object)
\S4method{$}{VariantFilteringParam}(x, name)
\S4method{names}{VariantFilteringParam}(x)
}
\arguments{
  \item{vcfFilenames}{
  Character string of the input VCF file name.
  }
  \item{pedFilename}{
  Character string of the pedigree file name in PED format.
  }
  \item{bsgenome}{
  Character string of a genome annotation package (\code{BSgenome.Hsapiens.1000genomes.hs37d5} by default).
  }
  \item{orgdb}{
  Character string of a gene-centric annotation package (\code{org.Hs.eg.db} by default).
  }
  \item{txdb}{
  Character string of a transcript-centric annotation package (\code{TxDb.Hsapiens.UCSC.hg19.knownGene} by default).
  The package \code{GenomicFeatures} provides infraestructure to build such annotation packages from different sources such as online UCSC tracks, Biomart tables, or \code{GFF} files.
  }
  \item{snpdb}{
  Character string of a SNP-centric annotation package (\code{SNPlocs.Hsapiens.dbSNP.20120608} by default).
  }
  \item{weightMatricesFilenames}{
  Character string of filenames of position weight matrices for binding site recognition.
  The default \code{NA} value indicates that no binding sites will be scored. To use this feature to score, for instance,
  splice sites in human, assign to this argument the function \code{spliceSiteMatricesHuman()}. See the files (\code{hsap.donors.hcmc10_15_1.ibn} and \code{hsap.acceptors.hcmc10_15_1.ibn}) returned by this function for details on their format.
  }
  \item{weightMatricesLocations}{
  Keywords of the annotated locations to variants under which a weight matrix will be used for scoring binding sites. This argument
  is only used when \code{weightMatricesFilenames!=NA} and, in such case, then more than one matrix is provided, this argument should be
  a list of character vectors with as many elements as matrices given in \code{weightMatricesFilenames}. The possible values can be obtained
  by typing \code{variantLocations()}.
  }
  \item{weightMatricesStrictLocations}{
  Logical vector flagging whether a weight matrix should be scoring binding sites strictly within the boundaries of the given locations.
  This argument is only used when \code{weightMatricesFilenames!=NA} and, in such case, then more than one matrix is provided, this
  argument should be a list of logical vectors with as many elements as matrices given in \code{weightMatricesFilenames}.
  }
  \item{radicalAAchangeFilename}{
  Name of a tab-separated text file containing chemical properties of amino acids. These properties are interpreted such that amino acid changes within a property are considered "conservative" and between properties are considered "radical". See the default file (\code{AA_chemical_properties_HanadaGojoboriLi2006.tsv}) for details on its format.
  }
  \item{codonusageFilename}{
  Name of a text file containing the codon usage.
  }
  \item{geneticCode}{
  Name character vector of length 64 describing the genetic code. The default value is \code{getGeneticCode("SGC0")}, the standard genetic code. An alternative genetic code, for instance, is \code{getGeneticCode("SGC1")}, the vertebrate mitochondrial genetic code. See \code{\link[Biostrings]{getGeneticCode}} in the Biostrings package for further details.
  }
  \item{allTranscripts}{
  Logical. This option allows the user to choose between working with all the transcripts affected by the variant (\code{allTranscripts=TRUE}) or with only one transcript per variant.
  }
  \item{regionAnnotations}{
  List of \code{\link[VariantAnnotation]{VariantType-class}} objects defining what regions to annotate.
  }
  \item{otherAnnotations}{
  Character vector of names of annotation packages or annotation objects.
  }
  \item{geneKeytype}{
  Character vector of the type of key gene identifier provided by the transcript-centric (TxDb) annotation package to interrogate the organism-centric (OrgDb) annotation package. The default value (\code{NA_character_} indicates that it will be assumed to be an Entrez identifier unless the values in the \code{GENEID} column returned by the TxDb package start with \code{ENSG} and then it will be assumed that they are Ensembl gene identifiers, or with one of \code{NM_, NP_, NR_, XM_, XP_, XR_ or YP_} and then it will be assumed that they are RefSeq gene identifiers.
  }
  \item{yieldSize}{
  Number of variants to yield each time the input VCF file is read. This argument is passed to the \code{\link[Rsamtools]{TabixFile}} function when opening the input VCF file and it allows to iterate through the variants in chunks of the given size to limit the memory requirements. Its default value (\code{NA_integer_}) implies that the whole input VCF file will be read into main memory.
  }
  \item{object}{
  A VariantFilteringParam object created through \code{VariantFilteringParam()}.
  }
  \item{x}{
  A VariantFilteringParam object created through \code{VariantFilteringParam()}.
  }
  \item{name}{
  Slot name of a VariantFilteringParam object. Use \code{names()} to find out what these slots are.
  }
}
\details{
The class \code{VariantFilteringParam} serves as a purpose of simplifying the call to the inheritance model function and its subsequent annotation and filtering steps. It also groups all the parameters that the user can customize (i.e newer versions of the annotation packages, when available). 

The method \code{VariantFilteringParam()} creates an \code{VariantFilteringParam} object used as an input argument to other functions such as \code{autosomalRecessiveHomozygous()}, etc.

The method \code{names()} allows one to see the names of the slots from a \code{VariantFilteringParam} object. Using the \code{$} operator, one can retrieve the values of these slots in an analogous way to a \code{list}.
}

\value{
An \code{VariantFilteringParam} object is returned by the method \code{VariantFilteringParam}.
}

\author{D.M. Elurbe, P. Puigdevall and R. Castelo}

\examples{
vfpar <- VariantFilteringParam(system.file("extdata", "CEUtrio.vcf.bgz", package="VariantFiltering"),
                               system.file("extdata", "CEUtrio.ped", package="VariantFiltering"),
                               snpdb=list(), otherAnnotations=character(0))
vfpar
names(vfpar)
vfpar$vcfFiles
}

\keyword{classes,methods}
