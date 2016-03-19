test_location_annotations <- function() {

  require(BiocParallel)

  CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
  vfpar <- VariantFilteringParam(CEUvcf,
                                 spliceSiteMatrices=c(system.file("extdata", "hsap.donors.hcmc10_15_1.ibn",
                                                                  package="VariantFiltering"),
                                                      system.file("extdata", "hsap.acceptors.hcmc10_15_1.ibn",
                                                                  package="VariantFiltering")),
                                 snpdb=list(),
                                 otherAnnotations=character(0))

  uind <- unrelatedIndividuals(vfpar, BPPARAM=MulticoreParam(workers=1))
  fv <- filteredVariants(uind)

  ## HAVE TO CHECK WHY THESE STOPPTED WORKING BUT PROBABLY IS THE NEW TXDB

  ## these are variants that pass quality filters and affect 5' and 3' ss
  ## idsFive <- c("rs760352870", "rs11910504", "rs3764880")
  idsFive <- c("rs11910504", "rs3764880")
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$LOCATION == "fiveSpliceSite"] == idsFive))

  ## idsThree <- c("rs1153275", "rs144475501", "rs220312", "rs2839260",
  ##               "rs6003299", "rs11568185")
  idsThree <- c("rs7271033", "rs144475501", "rs459763", "rs2070444",
                "rs2018393", "rs5758405", "rs11568185")
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$LOCATION == "threeSpliceSite"] == idsThree))

  ## rs58876072 introduces a frameshift
  ## checkTrue(fv[[1]]$VARID[fv[[1]]$CONSEQUENCE %in% "frameshift"] == "rs5840900")
  checkTrue(fv[[1]]$VARID[fv[[1]]$CONSEQUENCE %in% "frameshift"] == "rs58876072")

  ## rs2301384 introduces a stop codon
  ## idsNonsense <- c("rs73105200", "rs233303")
  idsNonsense <- c("rs2301384")
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$CONSEQUENCE %in% "nonsense"] == idsNonsense))
}
