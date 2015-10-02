[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/VariantFiltering.svg)](http://bioconductor.org/packages/release/bioc/html/VariantFiltering.html "Bioconductor status")

[![Bioconductor Availability](http://bioconductor.org/shields/availability/release/VariantFiltering.svg)](http://bioconductor.org/packages/release/bioc/html/VariantFiltering.html#archives "Platform availability") 
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/VariantFiltering.svg)](http://bioconductor.org/packages/stats/bioc/VariantFiltering.html "Percentile downloads")
[![Bioconductor Commits](http://bioconductor.org/shields/commits/bioc/VariantFiltering.svg)](http://bioconductor.org/packages/release/bioc/html/VariantFiltering.html#svn_source "svn commits")
[![Support posts](http://bioconductor.org/shields/posts/VariantFiltering.svg)](https://support.bioconductor.org/t/VariantFiltering/ "Bioconductor support posts")

[![Bioconductor Release Build](http://bioconductor.org/shields/build/release/bioc/VariantFiltering.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/VariantFiltering/ "Bioconductor release build")
[![Bioconductor Devel Build](http://bioconductor.org/shields/build/devel/bioc/VariantFiltering.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/VariantFiltering/ "Bioconductor devel build")

# VariantFiltering: Filtering of coding and non-coding genetic variants

## Installation

This is the __development__ version of the R/Bioconductor package VariantFiltering. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [http://www.bioconductor.org/packages/release/bioc/html/VariantFiltering.html](http://www.bioconductor.org/packages/release/bioc/html/VariantFiltering.html) and follow the instructions there to install it.

If you were really looking for this development version, then to install it you
need first to install the development version of R that you can find at [http://cran.r-project.org](http://cran.r-project.org) and then type the following instructions from the R shell:

```r
source("http://bioconductor.org/biocLite.R")
library(BiocInstaller)
useDevel()
biocLite("VariantFiltering")
```

Alternatively, you can install it from GitHub using
the [devtools](https://github.com/hadley/devtools "devtools") package.

```r
install.packages("devtools")
library(devtools)
install_github("VariantFiltering", "rcastelo")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **VariantFiltering**
please use the [Bioconductor support site][http://support.bioconductor.org "Bioconductor support site"].
For bug reports and issues regarding this __development__ version of **VariantFiltering**
please use the GitHub issues link at the top-right of this page
([https://github.com/rcastelo/VariantFiltering/issues](https://github.com/rcastelo/VariantFiltering/issues)).
