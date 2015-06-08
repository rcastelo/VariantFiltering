## additional VariantType methods for FiveSpliceSiteVariants and ThreeSpliceSiteVariants objects

## show methods

setMethod("show", "FiveSpliceSiteVariants",
          function(object) {
            cat("class:", class(object), "\n")
            cat("minIntronLength:", minIntronLength(object), "\n")
            cat("upstream:", upstream(object), "\n")
            cat("downstream:", downstream(object), "\n")
          })

setMethod("show", "ThreeSpliceSiteVariants",
          function(object) {
            cat("class:", class(object), "\n")
            cat("minIntronLength:", minIntronLength(object), "\n")
            cat("upstream:", upstream(object), "\n")
            cat("downstream:", downstream(object), "\n")
          })

## constructor methods

## 5' splice sites occur by default in introns that are at least 20nt long,
## have 3nt upstream and 4nt downstream from the first dinucleotide of the intron
FiveSpliceSiteVariants <- function(minIntronLength=20L,
                                   upstream=3, downstream=4) {
  minIntronLength <- VariantAnnotation:::.checkArgs(minIntronLength, "minIntronLength")
  upstream <- VariantAnnotation:::.checkArgs(upstream, "upstream")
  downstream <- VariantAnnotation:::.checkArgs(downstream, "downstream")

  new("FiveSpliceSiteVariants", minIntronLength=minIntronLength,
      upstream=upstream, downstream=downstream)
}

## 3' splice sites occur by default in introns that are at least 20nt long,
## have 18nt upstream and 3nt downstream from the last dinucleotide of the intron
ThreeSpliceSiteVariants <- function(minIntronLength=20L,
                                   upstream=18, downstream=3) {
  minIntronLength <- VariantAnnotation:::.checkArgs(minIntronLength, "minIntronLength")
  upstream <- VariantAnnotation:::.checkArgs(upstream, "upstream")
  downstream <- VariantAnnotation:::.checkArgs(downstream, "downstream")

  new("ThreeSpliceSiteVariants", minIntronLength=minIntronLength,
      upstream=upstream, downstream=downstream)
}

## getter and setters

setMethod("minIntronLength", "FiveSpliceSiteVariants",
          function(x) slot(x, "minIntronLength"))

setReplaceMethod("minIntronLength", "FiveSpliceSiteVariants",
                 function(x, value) {
                   slot(x, "minIntronLength") <- VariantAnnotation:::.checkArgs(value, "minIntronLength")
                   x
                 })

setMethod("upstream", "FiveSpliceSiteVariants",
          function(x) slot(x, "upstream"))

setReplaceMethod("upstream", "FiveSpliceSiteVariants",
                 function(x, value) {
                   slot(x, "upstream") <- VariantAnnotation:::.checkArgs(value, "upstream")
                   x
                 })

setMethod("downstream", "FiveSpliceSiteVariants",
          function(x) slot(x, "downstream"))

setReplaceMethod("downstream", "FiveSpliceSiteVariants",
                 function(x, value) {
                   slot(x, "downstream") <- VariantAnnotation:::.checkArgs(value, "downstream")
                   x
                 })

setMethod("minIntronLength", "ThreeSpliceSiteVariants",
          function(x) slot(x, "minIntronLength"))

setReplaceMethod("minIntronLength", "ThreeSpliceSiteVariants",
                 function(x, value) {
                   slot(x, "minIntronLength") <- VariantAnnotation:::.checkArgs(value, "minIntronLength")
                   x
                 })

setMethod("upstream", "ThreeSpliceSiteVariants",
          function(x) slot(x, "upstream"))

setMethod("downstream", "ThreeSpliceSiteVariants",
          function(x) slot(x, "downstream"))

setReplaceMethod("upstream", "ThreeSpliceSiteVariants",
                 function(x, value) {
                   slot(x, "upstream") <- VariantAnnotation:::.checkArgs(value, "upstream")
                   x
                 })

setReplaceMethod("downstream", "ThreeSpliceSiteVariants",
                 function(x, value) {
                   slot(x, "downstream") <- VariantAnnotation:::.checkArgs(value, "downstream")
                   x
                 })
