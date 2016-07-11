test_location_annotations <- function() {
  CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
  vfpar <- VariantFilteringParam(CEUvcf,
                                 weightMatricesFilenames=spliceSiteMatricesHuman(),
                                 snpdb=list(),
                                 otherAnnotations=character(0))

  uind <- unrelatedIndividuals(vfpar)
  fv <- filteredVariants(uind)

  ## these are variants that pass quality filters and affect 5' and 3' ss
  idsFive <- c("rs760352870", "rs11910504", "rs3764880")
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$LOCATION == "fiveSpliceSite"] == idsFive))

  idsThree <- c("rs1153275", "rs144475501", "rs220312", "rs2839260",
                "rs6003299", "rs11568185")
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$LOCATION == "threeSpliceSite"] == idsThree))

  ## rs5840900 introduces a frameshift
  checkTrue(fv[[1]]$VARID[fv[[1]]$CONSEQUENCE %in% "frameshift"] == "rs5840900")

  ## rs73105200 and rs233303 introduce a stop codon
  idsNonsense <- c("rs73105200", "rs233303")
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$CONSEQUENCE %in% "nonsense"] == idsNonsense))
}
