## create a pedigree with one affected chiled (A) and two
## unaffected parents (B and C)
pedDf <- data.frame(FamilyID="1", IndividualID=LETTERS[1:3],
                    FatherID=factor(c("B", NA, NA), levels=LETTERS[1:3]),
                    MotherID=factor(c("C", NA, NA), levels=LETTERS[1:3]),
                    Gender=c(1, 1, 2), Phenotype=c(2, 1, 1))

## create a VRanges object with two synthetic variants on three individuals
vr <- VRanges(seqnames=rep("chr1", 6),
              ranges=IRanges(c(rep(10, 3), rep(20, 3)),
                             c(rep(10, 3), rep(20, 3))),
              ref=c(rep("C", 3), rep("A", 3)), alt=c(rep("T", 3), rep("G", 3)),
              refDepth=c(rep(10, 3), rep(5, 3)), altDepth=c(rep(5, 3), rep(10, 3)),
              totalDepth=rep(100, 6), sampleNames=rep(LETTERS[1:3], 2),
              softFilterMatrix=FilterMatrix(matrix=cbind(LowQual=rep(TRUE, 6)),
                                            filterRules=FilterRules(LowQual=TRUE)),
              GT=NULL)

test_deNovo_inheritance <- function() {

  ## first variant segregates as a de novo trait and second one does not
  vr$GT <- c("0/1", "0/0", "0/0", "0/1", "0/0", "1/1")

  ## mask variants segregating as a de novo trait
  mask <- VariantFiltering:::.deNovoMask(vObj=vr, pedDf=pedDf)

  ## the first variant should segregate as de novo while the second should not
  checkTrue(all(mask[1:3]) && all(!mask[4:6]))
}

test_autosomalRecessiveHomozygous_inheritance <- function() {

  require(BSgenome.Hsapiens.UCSC.hg19)

  ## first variant segregates as an autosomal recessive homozygous trait and second one does not
  vr$GT <- c("1/1", "0/1", "0/1", "0/1", "0/1", "0/1")

  ## mask variants segregating as a de novo trait
  mask <- VariantFiltering:::.autosomalRecessiveHomozygousMask(vObj=vr,
                                                               bsgenome=Hsapiens,
                                                               pedDf=pedDf)

  ## the first variant should segregate as de novo while the second should not
  checkTrue(all(mask[1:3]) && all(!mask[4:6]))
}

test_autosomalDominant_inheritance <- function() {

  require(BSgenome.Hsapiens.UCSC.hg19)

  ## first variant segregates as an autosomal recessive homozygous trait and second one does not
  pedDf$Phenotype[pedDf$IndividualID == "B"] <- 2
  vr$GT <- c("0/1", "1/1", "0/0", "0/1", "0/0", "0/0")

  ## mask variants segregating as a de novo trait
  mask <- VariantFiltering:::.autosomalDominantMask(vObj=vr,
                                                    bsgenome=Hsapiens,
                                                    pedDf=pedDf)

  ## the first variant should segregate as de novo while the second should not
  checkTrue(all(mask[1:3]) && all(!mask[4:6]))
}
