## based on the SNPlocs class defined in the BSgenome package
## maybe this class could be generalized to one whose data correspond
## to numerical values associated to nucleotide positions and
## sotored in a RleList object
setClass("MafDb",
         representation=representation(provider="character", ## e.g. UCSC
                                       provider_version="character", ## creation date in compact format
                                       download_url="character", ## download URL of the ExAC data
                                       download_date="character", ## date on which data was downloaded
                                       reference_genome="GenomeDescription", ## extracted from BSgenome.*

                                       ## package name and absolute path to local directory where to find
                                       ## the serialized objects containing the data
                                       data_pkgname="character",
                                       data_dirpath="character",
                                       data_serialized_objnames="character",
                                       data_tag="character",
                                       data_pops="character",

                                       .data_cache="environment"))


## based on the SNPlocs class defined in the BSgenome package
## maybe this class could be generalized to one whose data correspond
## to numerical values associated to nucleotide positions and
## sotored in a RleList object
setClass("MafDb2",
         representation=representation(provider="character", ## e.g. UCSC
                                       provider_version="character", ## creation date in compact format
                                       download_url="character", ## download URL of the ExAC data
                                       download_date="character", ## date on which data was downloaded
                                       reference_genome="GenomeDescription", ## extracted from BSgenome.*

                                       ## package name and absolute path to local directory where to find
                                       ## the serialized objects containing the data
                                       data_pkgname="character",
                                       data_dirpath="character",
                                       data_serialized_objnames="character",
                                       data_pops="character",

                                       .data_cache="environment"))

## based on the SNPlocs class defined in the BSgenome package
## maybe this class could be generalized to one whose data correspond
## to numerical values associated to nucleotide positions and
## sotored in a RleList object
setClass("PhastConsDb",
         representation=representation(provider="character", ## e.g. UCSC
                                       provider_version="character", ## creation date in compact format
                                       download_url="character", ## download URL of all phastCons data
                                       download_date="character", ## date on which data was downloaded
                                       reference_genome="GenomeDescription", ## extracted from BSgenome.*

                                       ## package name and absolute path to local directory where to find
                                       ## the serialized objects containing the phastCons scores
                                       data_pkgname="character",
                                       data_dirpath="character",
                                       data_serialized_objnames="character",

                                       .data_cache="environment"))

setClass("GenePhylostrataDb",
         representation=representation(organism="character",
                                       download_url="character",
                                       download_date="character",
                                       source_pub="character",

                                       ## package name and absolute path to local directory where to find
                                       ## the serializd object containing the phastCons scores
                                       data_pkgname="character",
                                       data_dirpath="character",
                                       data_serialized_objnames="character",

                                       ## values defining the strata
                                       strataNames="character",
                                       strataTaxID="character",

                                       .data_cache="environment"))

## define a class in the VariantAnnotation VariantType hierarchy to
## annotate separately 5' and 3' splice sites and be able to set user-defined
## boundaries beyond the default canonical di-nucleotides
setClass("FiveSpliceSiteVariants", contains="SpliceSiteVariants",
         representation(minIntronLength="integer",
                        upstream="integer",
                        downstream="integer"))

setClass("ThreeSpliceSiteVariants", contains="SpliceSiteVariants",
         representation(minIntronLength="integer",
                        upstream="integer",
                        downstream="integer"))


setClass("VariantFilteringParam",
         representation=representation(callObj="call",
                                       callStr="character",
                                       vcfFiles="TabixFileList",
                                       seqInfos="list",
                                       sampleNames="character",
                                       pedFilename="character",
                                       bsgenome="BSgenome",
                                       orgdb="OrgDb",
                                       txdb="TxDb",
                                       snpdb="list",
                                       gSO="graphNEL",
                                       gSOdmat="matrix",
                                       gSOamat="matrix",
                                       weightMatrices="list",
                                       radicalAAchangeFilename="character",
                                       radicalAAchangeMatrix="matrix",
                                       geneticCode="character",
                                       regionAnnotations="SimpleList",
                                       otherAnnotations="list",
                                       allTranscripts="logical",
                                       codonusageFilename="character",
                                       codonusageTable="numeric",
                                       geneKeytype="character",
                                       filters="FilterRules",
                                       filterDescriptions="DataFrame",
                                       qualityFilterNames="character",
                                       cutoffs="list",
                                       yieldSize="integer"))

setClass("VariantFilteringResults",
         representation=representation(callObj="call",
                                       callStr="character",
                                       ## inputParameters="VariantFilteringParam",
                                       inputParameters="ANY",
                                       activeSamples="character",
                                       inheritanceModel="character",
                                       variants="VRanges",
                                       bamViews="BamViews",
                                       gSO="graphNEL",
                                       filters="FilterRules",
                                       cutoffs="list",
                                       dbSNPflag="character",
                                       OMIMflag="character",
                                       variantTypeMask="logical",
                                       locationMask="logical",
                                       consequenceMask="logical",
                                       aaChangeType="character",
                                       MAFpopMask="logical",
                                       naMAF="logical",
                                       maxMAF="numeric",
                                       minPhastCons="numeric",
                                       minPhylostratumIndex="integer",
                                       minScore5ss="numeric",
                                       minScore3ss="numeric",
                                       minCUFC="numeric"))

setClass("WeightMatrix",
         representation(wm="externalptr",
                        filename="character",
                        locations="character",
                        strictLocations="logical"))
