## for MafDb class
setGeneric("dbConn", function(x) standardGeneric("dbConn"))
setGeneric("knownVariantsMAFcols", function(mafdb, ...) standardGeneric("knownVariantsMAFcols"))
setGeneric("fetchKnownVariantsByID", function(mafdb, ...) standardGeneric("fetchKnownVariantsByID"))

## for PhastConsDb class
setGeneric("scores", function(object, gpos, ...) standardGeneric("scores"))

## for GenePhylostrataDb class
setGeneric("genePhylostratum", function(object, ids) standardGeneric("genePhylostratum"))
setGeneric("genePhylostrata", function(object) standardGeneric("genePhylostrata"))

## for VariantFilteringResults class
setGeneric("VariantFilteringResults", function(variants, ...) standardGeneric("VariantFilteringResults"))
setGeneric("inheritanceModel", function(x) standardGeneric("inheritanceModel"))
setGeneric("param", function(x) standardGeneric("param"))
setGeneric("allVariants", function(x, ...) standardGeneric("allVariants"))
setGeneric("filteredVariants", function(x, ...) standardGeneric("filteredVariants"))
setGeneric("dbSNPpresent", function(x) standardGeneric("dbSNPpresent"))
setGeneric("dbSNPpresent<-", function(x, value) standardGeneric("dbSNPpresent<-"))
setGeneric("OMIMpresent", function(x) standardGeneric("OMIMpresent"))
setGeneric("OMIMpresent<-", function(x, value) standardGeneric("OMIMpresent<-"))
setGeneric("variantType", function(x) standardGeneric("variantType"))
setGeneric("variantType<-", function(x, value) standardGeneric("variantType<-"))
setGeneric("aaChangeType", function(x) standardGeneric("aaChangeType"))
setGeneric("aaChangeType<-", function(x, value) standardGeneric("aaChangeType<-"))
setGeneric("naMAF", function(x) standardGeneric("naMAF"))
setGeneric("naMAF<-", function(x, value) standardGeneric("naMAF<-"))
setGeneric("MAFpop", function(x) standardGeneric("MAFpop"))
setGeneric("MAFpop<-", function(x, popkey=NA, value) standardGeneric("MAFpop<-"))
setGeneric("maxMAF", function(x) standardGeneric("maxMAF"))
setGeneric("maxMAF<-", function(x, value) standardGeneric("maxMAF<-"))
setGeneric("minPhastCons", function(x) standardGeneric("minPhastCons"))
setGeneric("minPhastCons<-", function(x, value) standardGeneric("minPhastCons<-"))
setGeneric("minPhylostratum", function(x) standardGeneric("minPhylostratum"))
setGeneric("minPhylostratum<-", function(x, value) standardGeneric("minPhylostratum<-"))
setGeneric("minCRYP5ss", function(x) standardGeneric("minCRYP5ss"))
setGeneric("minCRYP5ss<-", function(x, value) standardGeneric("minCRYP5ss<-"))
setGeneric("minCRYP3ss", function(x) standardGeneric("minCRYP3ss"))
setGeneric("minCRYP3ss<-", function(x, value) standardGeneric("minCRYP3ss<-"))

## for VariantFilteringResultsAIM class
setGeneric("selectInheritancePattern", function(x) standardGeneric("selectInheritancePattern"))
setGeneric("selectInheritancePattern<-", function(x, value) standardGeneric("selectInheritancePattern<-"))
setGeneric("selectIndexCase", function(x) standardGeneric("selectIndexCase"))
setGeneric("selectIndexCase<-", function(x, value) standardGeneric("selectIndexCase<-"))

setGeneric("selectCarrierRelative1", function(x) standardGeneric("selectCarrierRelative1"))
setGeneric("selectCarrierRelative1<-", function(x, value) standardGeneric("selectCarrierRelative1<-"))
setGeneric("selectCarrierRelative2", function(x) standardGeneric("selectCarrierRelative2"))
setGeneric("selectCarrierRelative2<-", function(x, value) standardGeneric("selectCarrierRelative2<-"))
setGeneric("selectAffectedRelative", function(x) standardGeneric("selectAffectedRelative"))
setGeneric("selectAffectedRelative<-", function(x, value) standardGeneric("selectAffectedRelative<-"))

setGeneric("selectCarrierAllele1CH", function(x) standardGeneric("selectCarrierAllele1CH"))
setGeneric("selectCarrierAllele1CH<-", function(x, value) standardGeneric("selectCarrierAllele1CH<-"))
setGeneric("selectCarrierAllele2CH", function(x) standardGeneric("selectCarrierAllele2CH"))
setGeneric("selectCarrierAllele2CH<-", function(x, value) standardGeneric("selectCarrierAllele2CH<-"))
setGeneric("selectAffRelative1CH", function(x) standardGeneric("selectAffRelative1CH"))
setGeneric("selectAffRelative1CH<-", function(x, value) standardGeneric("selectAffRelative1CH<-"))

setGeneric("selectUnaffectedRelative1AD", function(x) standardGeneric("selectUnaffectedRelative1AD"))
setGeneric("selectUnaffectedRelative1AD<-", function(x, value) standardGeneric("selectUnaffectedRelative1AD<-"))
setGeneric("selectUnaffectedRelative2AD", function(x) standardGeneric("selectUnaffectedRelative2AD"))
setGeneric("selectUnaffectedRelative2AD<-", function(x, value) standardGeneric("selectUnaffectedRelative2AD<-"))
setGeneric("selectAffectedRelative1AD", function(x) standardGeneric("selectAffectedRelative1AD"))
setGeneric("selectAffectedRelative1AD<-", function(x, value) standardGeneric("selectAffectedRelative1AD<-"))

setGeneric("selectCarrRelFem1XL", function(x) standardGeneric("selectCarrRelFem1XL"))
setGeneric("selectCarrRelFem1XL<-", function(x, value) standardGeneric("selectCarrRelFem1XL<-"))
setGeneric("selectAffRelMale1XL", function(x) standardGeneric("selectAffRelMale1XL"))
setGeneric("selectAffRelMale1XL<-", function(x, value) standardGeneric("selectAffRelMale1XL<-"))
setGeneric("selectUnaffMale1XL", function(x) standardGeneric("selectUnaffMale1XL"))
setGeneric("selectUnaffMale1XL<-", function(x, value) standardGeneric("selectUnaffMale1XL<-"))

setGeneric("selectParent1DN", function(x) standardGeneric("selectParent1DN"))
setGeneric("selectParent1DN<-", function(x, value) standardGeneric("selectParent1DN<-"))
setGeneric("selectParent2DN", function(x) standardGeneric("selectParent2DN"))
setGeneric("selectParent2DN<-", function(x, value) standardGeneric("selectParent2DN<-"))

## for VariantFilteringResultsUI class
setGeneric("selectIndividual", function(x) standardGeneric("selectIndividual"))
setGeneric("selectIndividual<-", function(x, value) standardGeneric("selectIndividual<-"))
setGeneric("selectGene", function(x) standardGeneric("selectGene"))
setGeneric("selectGene<-", function(x, value) standardGeneric("selectGene<-"))

## for all three VariantFilteringResults
setGeneric("reportVariants", function(vfResultsObj, ...) standardGeneric("reportVariants"))

## for the rest ..
setGeneric("VariantFilteringParam", function(vcfFilenames, ...) standardGeneric("VariantFilteringParam"))
setGeneric("xLinked", function(param, ...) standardGeneric("xLinked"))
setGeneric("autosomalRecessiveHomozygous", function(param, ...) standardGeneric("autosomalRecessiveHomozygous"))
setGeneric("autosomalDominant", function(param, ...) standardGeneric("autosomalDominant"))
setGeneric("deNovo", function(param, ...) standardGeneric("deNovo"))
setGeneric("autosomalRecessiveHeterozygous", function(param, ...) standardGeneric("autosomalRecessiveHeterozygous"))
setGeneric("unrelatedIndividuals", function(param, ...) standardGeneric("unrelatedIndividuals"))
setGeneric("allInheritanceModels", function(param, ...) standardGeneric("allInheritanceModels"))
setGeneric("annotateVariants", function(annObj, ...) standardGeneric("annotateVariants"))
setGeneric("wmScore", function(object, dnaseqs, ...) standardGeneric("wmScore"))
setGeneric("conservedPositions", function(x, ...) standardGeneric("conservedPositions"))
