checkTrue(require(BSgenome.Hsapiens.1000genomes.hs37d5))
Hsapiens <- BSgenome.Hsapiens.1000genomes.hs37d5

## create a pedigree with one affected male child (A) and two
## unaffected parents (B -male- and C -female-)
pedDf <- data.frame(FamilyID="1", IndividualID=LETTERS[1:3],
                    FatherID=factor(c("B", NA, NA), levels=LETTERS[1:3]),
                    MotherID=factor(c("C", NA, NA), levels=LETTERS[1:3]),
                    Sex=c(1, 1, 2), Phenotype=c(2, 1, 1))

## create a VRanges object with two synthetic variants on three individuals
vr <- VRanges(seqnames=rep("1", 6),
              ranges=IRanges(c(rep(10, 3), rep(20, 3)),
                             c(rep(10, 3), rep(20, 3))),
              ref=c(rep("C", 3), rep("A", 3)), alt=c(rep("T", 3), rep("G", 3)),
              refDepth=c(rep(10, 3), rep(5, 3)), altDepth=c(rep(5, 3), rep(10, 3)),
              totalDepth=rep(100, 6), sampleNames=rep(LETTERS[1:3], 2),
              softFilterMatrix=FilterMatrix(matrix=cbind(LowQual=rep(TRUE, 6)),
                                            filterRules=FilterRules(LowQual=TRUE)),
              GT=NULL,
              seqinfo=seqinfo(Hsapiens))

test_deNovo_inheritance <- function() {

  thisvr <- vr

  ## first variant segregates as a de novo trait and second one does not
  thisvr$GT <- c("0/1", "0/0", "0/0", "0/1", "0/0", "1/1")

  ## mask variants segregating as a de novo trait
  mask <- VariantFiltering:::.deNovoMask(vObj=thisvr, pedDf=pedDf)

  ## the first variant should segregate as de novo while the second should not
  checkTrue(all(mask[1:3]) && all(!mask[4:6]))
}

test_autosomalRecessiveHomozygous_inheritance <- function() {

  thisvr <- vr

  ## first variant segregates as an autosomal recessive homozygous trait and second one does not
  thisvr$GT <- c("1/1", "0/1", "0/1", "0/1", "0/1", "0/1")

  ## mask variants segregating as an autosomal recessive homozygous trait
  mask <- VariantFiltering:::.autosomalRecessiveHomozygousMask(vObj=thisvr,
                                                               bsgenome=Hsapiens,
                                                               pedDf=pedDf)

  ## the first variant should segregate as autosomal recessive homozygous while the second should not
  checkTrue(all(mask[1:3]) && all(!mask[4:6]))
}

test_autosomalDominant_inheritance <- function() {

  thisvr <- vr

  ## first variant segregates as an autosomal dominant trait and second one does not
  pedDf$Phenotype[pedDf$IndividualID == "B"] <- 2
  thisvr$GT <- c("0/1", "1/1", "0/0", "0/1", "0/0", "0/0")

  ## mask variants segregating as an autosomal dominant trait
  mask <- VariantFiltering:::.autosomalDominantMask(vObj=thisvr,
                                                    bsgenome=Hsapiens,
                                                    pedDf=pedDf)

  ## the first variant should segregate as autosomal dominant while the second should not
  checkTrue(all(mask[1:3]) && all(!mask[4:6]))
  pedDf$Phenotype[pedDf$IndividualID == "B"] <- 1
}

test_xLinked_inheritance <- function() {

  thisvr <- vr

  ## first variant segregates as an X-linked trait and second one does not
  thisvr$GT <- c("1/1", "0/0", "0/1", "1/1", "0/1", "0/1")
  seqnames(thisvr) <- factor("X", levels=seqlevels(thisvr))

  ## mask variants segregating as an X-linked trait
  mask <- VariantFiltering:::.xLinkedMask(vObj=thisvr, bsgenome=Hsapiens, pedDf=pedDf)

  ## the first variant should segregate as X-linked while the second should not
  checkTrue(all(mask[1:3]) && all(!mask[4:6]))
}
