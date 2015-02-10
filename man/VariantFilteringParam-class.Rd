\name{VariantFilteringParam-class}
\docType{class}

%% Classes:
\alias{class:VariantFilteringParam}
\alias{VariantFilteringParam}
\alias{VariantFilteringParam-class}

%% Methods:
\alias{VariantFilteringParam,character-method}
\alias{show,VariantFilteringParam-method}
\alias{$,VariantFilteringParam-method}
\alias{names,VariantFilteringParam-method}

\title{VariantFiltering parameter class}

\description{
The class \code{VariantFilteringParam} is defined to ease configuring the call to the functions that filter input genetic variants according to a desired segregating inheritance model (\code{xLinked()}, \code{autosomalRecessiveHomozygous()}, etc).
}
\usage{
\S4method{VariantFilteringParam}{character}(vcfFilenames, pedFilename=NA_character_,
                                            bsgenome="BSgenome.Hsapiens.UCSC.hg19",
                                            orgdb="org.Hs.eg.db",
                                            txdb="TxDb.Hsapiens.UCSC.hg19.knownGene",
                                            snpdb="SNPlocs.Hsapiens.dbSNP.20120608",
                                            spliceSiteMatricesFilenames=c(file.path(system.file("extdata", package="VariantFiltering"), "hsap.donors.hcmc10_15_1.ibn"),
                                                                          file.path(system.file("extdata", package = "VariantFiltering"), "hsap.acceptors.hcmc10_15_1.ibn")),
                                            radicalAAchangeFilename=file.path(system.file("extdata", package="VariantFiltering"),
                                                                              "AA_chemical_properties_HanadaGojoboriLi2006.tsv"),
                                            codonusageFilename=file.path(system.file("extdata", package="VariantFiltering"),"humanCodonUsage.txt"),
                                            geneticCode=getGeneticCode("SGC0"),
                                            allTranscripts=FALSE,
                                            otherAnnotations=c("MafDb.ESP6500SI.V2.SSA137.dbSNP138",
                                                               "MafDb.ALL.wgs.phase1.release.v3.20101123",
                                                               "PolyPhen.Hsapiens.dbSNP131",
                                                               "SIFT.Hsapiens.dbSNP137",
                                                               "phastCons100way.UCSC.hg19",
                                                               "humanGenesPhylostrata"),
                                            filterTag=NA_character_,
                                            geneKeytype=NA_character_,
                                            yieldSize=NA_integer_)
\S4method{show}{VariantFilteringParam}(object)
\S4method{$}{VariantFilteringParam}(x, name)
\S4method{names}{VariantFilteringParam}(x)
}
\arguments{
  \item{vcfFilenames}{
  Character vector of paths to VCF files.
  }
  \item{pedFilename}{
  Character string of the pedigree file name in PED format.
  }
  \item{bsgenome}{
  Character string of a genome annotation package (\code{BSgenome.Hsapiens.UCSC.hg19} by defautl).
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
  %\item{readDepthGTfield}{
  %Field name in the input VCF file that contains the total read depth at each position and sample. Default value is \code{DP}.
  %}
  %\item{allelicDepthGTfield}{
  %Field name in the input VCF file that contains the read depth of each allele per sample. Default value is \code{AD}.
  %}
  \item{spliceSiteMatricesFilenames}{
  Names of text files containing position weight matrices for splice site recognition. See the default files (\code{hsap.donors.hcmc10_15_1.ibn} and \code{hsap.acceptors.hcmc10_15_1.ibn}) for details on their format.
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
  \item{otherAnnotations}{
  Character vector of names of annotation packages or annotation objects.
  }
  \item{filterTag}{
  Character vector of tags used to select only those variants which present the corresponding string (or strings) within the VCF's FILTER column. \code{NA_character_} (default) means that it will consider all the variants present in the VCF file to perform the analysis. An example of possible value is PASS, a commonly used flag to indicate that the call fullfils all the quality parameters, so only variants flagged this way will be used in the analysis.
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
p <- VariantFilteringParam(list.files(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.gz$", full.names=TRUE), list.files(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped", full.names=TRUE))
p
names(p)
p$vcfFiles
}

\keyword{classes,methods}
