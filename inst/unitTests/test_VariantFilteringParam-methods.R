test_VariantFilteringParam <- function() {

  CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
  CEUped <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped")
  vfpar <- VariantFilteringParam(CEUvcf, CEUped, snpdb=character(0),
                                 otherAnnotations=character(0))

  checkTrue(class(get(vfpar$bsgenome)) == "BSgenome")
  checkTrue(class(get(vfpar$orgdb)) == "OrgDb")
  checkTrue(class(get(vfpar$txdb)) == "TxDb")
}
