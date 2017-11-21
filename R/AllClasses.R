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
                                       bsgenome="character",
                                       orgdb="character",
                                       txdb="character",
                                       snpdb="character",
                                       gSO="graphNEL",
                                       gSOdmat="matrix",
                                       gSOamat="matrix",
                                       weightMatrices="list",
                                       radicalAAchangeFilename="character",
                                       radicalAAchangeMatrix="matrix",
                                       geneticCode="character",
                                       regionAnnotations="SimpleList",
                                       otherAnnotations="character",
                                       otherAnnotationsClass="character",
                                       allTranscripts="logical",
                                       codonusageFilename="character",
                                       codonusageTable="numeric",
                                       geneKeytype="character",
                                       filters="FilterRules",
                                       filtersMetadata="DataFrame",
                                       qualityFilterNames="character",
                                       yieldSize="integer"))

setClass("CutoffsList", contains="SimpleList",
         prototype=prototype(elementType="ANY"))

setClass("VariantFilteringResults",
         representation=representation(callObj="call",
                                       callStr="character",
                                       genomeDescription="GenomeDescription",
                                       inputParameters="VariantFilteringParam",
                                       activeSamples="character",
                                       inheritanceModel="character",
                                       variants="VRanges",
                                       bamViews="BamViews",
                                       gSO="graphNEL",
                                       filters="FilterRules",
                                       filtersMetadata="DataFrame",
                                       cutoffs="CutoffsList",
                                       sortings="CutoffsList",
                                       annoGroups="list",
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

