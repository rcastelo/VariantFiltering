test_VariantFilteringParam <- function() {

  CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.vcf.bgz")
  CEUped <- file.path(system.file("extdata", package="VariantFiltering"), "CEUtrio.ped")
  vfpar <- VariantFilteringParam(CEUvcf, CEUped, snpdb=list(),
                                 otherAnnotations=character(0))

  checkTrue(class(vfpar$bsgenome) == "BSgenome")
  checkTrue(class(vfpar$orgdb) == "OrgDb")
  checkTrue(class(vfpar$txdb) == "TxDb")
}
