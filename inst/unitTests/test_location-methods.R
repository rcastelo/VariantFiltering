test_location_annotations <- function() {
  CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
  vfpar <- VariantFilteringParam(CEUvcf,
                                 weightMatricesFilenames=spliceSiteMatricesHuman(),
                                 snpdb=character(0),
                                 otherAnnotations=character(0))

  uind <- unrelatedIndividuals(vfpar)
  fv <- filteredVariants(uind)

  ## these are variants that pass quality filters and affect 5' and 3' ss
  ## idsFive <- c("rs760352870", "rs11910504", "rs3764880") TxDb.Hsapiens.UCSC.hg19.knownGene rel 3.21
  idsFive <- c("rs760352870", "rs34310744", "rs11910504", "rs3764880") ## rel 3.22
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$LOCATION == "fiveSpliceSite"] == idsFive))

  ## rel 3.21
  ## idsThree <- c("rs1153275", "rs144475501", "rs220312", "rs2839260",
  ##               "rs6003299", "rs11568185")
  idsThree <- c("rs1015518", "rs3761264", "rs1153275", "rs144475501", "rs220312",
		"rs1475891", "rs2839260", "rs11568185") ## rel 3.22
  checkTrue(all(unique(fv[[1]]$VARID[fv[[1]]$LOCATION == "threeSpliceSite"]) == idsThree))

  ## frameshifts
  ## idsFrameshift <- "rs5840900" ## rel 3.21
  idsFrameshift <- c("rs5840900", "rs140511")
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$CONSEQUENCE %in% "frameshift"] == idsFrameshift))

  ## rs73105200 and rs233303 introduce a stop codon
  ## idsNonsense <- c("rs73105200", "rs233303") rel 3.21
  idsNonsense <- "rs73105200"
  checkTrue(all(fv[[1]]$VARID[fv[[1]]$CONSEQUENCE %in% "nonsense"] == idsNonsense))
}
