### R code from vignette source 'usingVariantFiltering.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: <setup
###################################################
pdf.options(useDingbats=FALSE)
options(width=80)


###################################################
### code chunk number 3: usingVariantFiltering.Rnw:87-94
###################################################
library(VariantFiltering)

CEUvcf <- file.path(system.file("extdata", package="VariantFiltering"),
                    "CEUtrio.vcf.bgz")
vfpar <- VariantFilteringParam(CEUvcf)
class(vfpar)
vfpar


###################################################
### code chunk number 4: usingVariantFiltering.Rnw:101-102
###################################################
args(VariantFilteringParam)


###################################################
### code chunk number 5: usingVariantFiltering.Rnw:115-120
###################################################
if (VariantFiltering:::.getOS() == "windows") {
  uind <- VariantFiltering:::.loadPrecomputedVariantFilteringResults()
} else
  uind <- unrelatedIndividuals(vfpar) ## ~75 sec.



###################################################
### code chunk number 6: usingVariantFiltering.Rnw:121-122 (eval = FALSE)
###################################################
## uind <- unrelatedIndividuals(vfpar)


###################################################
### code chunk number 7: usingVariantFiltering.Rnw:124-126
###################################################
class(uind)
uind


###################################################
### code chunk number 8: usingVariantFiltering.Rnw:137-138
###################################################
summary(uind)


###################################################
### code chunk number 9: usingVariantFiltering.Rnw:152-153
###################################################
summary(uind, method="bioc")


###################################################
### code chunk number 10: usingVariantFiltering.Rnw:159-161
###################################################
uindSO <- summary(uind, method="SOfull")
uindSO


###################################################
### code chunk number 11: usingVariantFiltering.Rnw:168-169
###################################################
uindSO[uindSO$Level == 6, ]


###################################################
### code chunk number 12: usingVariantFiltering.Rnw:174-175
###################################################
allVariants(uind)


###################################################
### code chunk number 13: displayedVariant
###################################################
path2bams <- file.path(system.file("extdata", package="VariantFiltering"),
                       paste0(samples(uind), ".subset.bam"))
bv <- BamViews(bamPaths=path2bams,
               bamSamples=DataFrame(row.names=samples(uind)))
bamFiles(uind) <- bv
bamFiles(uind)
plot(uind, what="rs6130959", sampleName="NA12892")


###################################################
### code chunk number 14: usingVariantFiltering.Rnw:211-216
###################################################
samples(uind)
samples(uind) <- c("NA12891", "NA12892")
uind
uindSO2sam <- summary(uind, method="SOfull")
uindSO2sam[uindSO2sam$Level == 6, ]


###################################################
### code chunk number 15: usingVariantFiltering.Rnw:222-224
###################################################
uind <- resetSamples(uind)
uind


###################################################
### code chunk number 16: usingVariantFiltering.Rnw:235-236
###################################################
filters(uind)


###################################################
### code chunk number 17: usingVariantFiltering.Rnw:242-243
###################################################
active(filters(uind))


###################################################
### code chunk number 18: usingVariantFiltering.Rnw:248-251
###################################################
active(filters(uind)) <- TRUE
active(filters(uind))
summary(uind)


###################################################
### code chunk number 19: usingVariantFiltering.Rnw:256-259
###################################################
active(filters(uind)) <- FALSE
active(filters(uind))["dbSNP"] <- TRUE
summary(uind)


###################################################
### code chunk number 20: usingVariantFiltering.Rnw:268-272
###################################################
change(cutoffs(uind), "SOterms") <- "UTR_variant"
active(filters(uind))["SOterms"] <- TRUE
summary(uind)
summary(uind, method="SOfull")


###################################################
### code chunk number 21: usingVariantFiltering.Rnw:284-295
###################################################
startLost <- function(x) {
  mask <- start(allVariants(x, groupBy="nothing")$CDSLOC) == 1 &
          as.character(allVariants(x, groupBy="nothing")$REFCODON) == "ATG" &
          as.character(allVariants(x, groupBy="nothing")$VARCODON) != "ATG"
  mask
}
filters(uind)$startLost <- startLost
active(filters(uind)) <- FALSE
active(filters(uind))["startLost"] <- TRUE
active(filters(uind))
summary(uind)


###################################################
### code chunk number 22: usingVariantFiltering.Rnw:307-318
###################################################
active(filters(uind)) <- FALSE
active(filters(uind))["maxMAF"] <- TRUE
cutoffs(uind)$maxMAF
cutoffs(uind)$maxMAF$popmask
cutoffs(uind)$maxMAF$maxvalue
change(cutoffs(uind)$maxMAF, "popmask") <- FALSE
cutoffs(uind)$maxMAF$popmask
change(cutoffs(uind)$maxMAF, "popmask") <- c(AFKGp1=TRUE)
cutoffs(uind)$maxMAF$popmask
change(cutoffs(uind)$maxMAF, "maxvalue") <- 0.01
summary(uind)


###################################################
### code chunk number 23: usingVariantFiltering.Rnw:324-325
###################################################
filteredVariants(uind)


###################################################
### code chunk number 24: usingVariantFiltering.Rnw:336-339 (eval = FALSE)
###################################################
## CEUped <- file.path(system.file("extdata", package="VariantFiltering"),
##                     "CEUtrio.ped")
## param <- VariantFilteringParam(vcfFilenames=CEUvcf, pedFilename=CEUped)


###################################################
### code chunk number 25: reportvariants (eval = FALSE)
###################################################
## reportVariants(uind, type="csv", file="uind.csv")


###################################################
### code chunk number 26: info
###################################################
toLatex(sessionInfo())


