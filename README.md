# VariantFiltering: Filtering of coding and non-coding genetic variants

[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/VariantFiltering.svg)](http://bioconductor.org/packages/release/bioc/html/VariantFiltering.html "How long has been VariantFiltering in a release of Bioconductor")
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/VariantFiltering.svg)](http://bioconductor.org/packages/stats/bioc/VariantFiltering.html "Percentile (top 5/20/50% or 'available') of downloads over the last 6 full months")
[![Bioconductor Commits](http://bioconductor.org/shields/commits/bioc/VariantFiltering.svg)](http://bioconductor.org/packages/devel/bioc/html/VariantFiltering.html#svn_source "Average SVN commits (to the devel branch) per month over the last 6 months")
[![Support posts](http://bioconductor.org/shields/posts/VariantFiltering.svg)](https://support.bioconductor.org/t/VariantFiltering/ "Bioconductor support site activity on VariantFiltering, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")

**Current build status**
- `release` [![Bioconductor Availability](http://bioconductor.org/shields/availability/release/VariantFiltering.svg)](http://bioconductor.org/packages/release/bioc/html/VariantFiltering.html#archives "Whether VariantFiltering release is available on all platforms") 
[![Bioconductor Release Build](http://bioconductor.org/shields/build/release/bioc/VariantFiltering.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/VariantFiltering/ "Bioconductor release build")
- `development` [![Bioconductor Availability](http://bioconductor.org/shields/availability/devel/VariantFiltering.svg)](http://bioconductor.org/packages/devel/bioc/html/VariantFiltering.html#archives "Whether VariantFiltering devel is available on all platforms") 
[![Bioconductor Devel Build](http://bioconductor.org/shields/build/devel/bioc/VariantFiltering.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/VariantFiltering/ "Bioconductor devel build")

## Installation

This is the __development__ version of the R/Bioconductor package VariantFiltering. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [http://bioconductor.org/packages/VariantFiltering](http://bioconductor.org/packages/VariantFiltering) and follow the instructions there to install it.

If you were really looking for this development version, then to install it you
need first to install the development version of R that you can find at [http://cran.r-project.org](http://cran.r-project.org) and then type the following instructions from the R shell:

```r
install.packages("BiocManager")
BiocManager::install("BiocInstaller", version="devel")
BiocInstaller::useDevel()
BiocManager::install("VariantFiltering")
```

Alternatively, you can install it from GitHub using the [devtools](https://github.com/hadley/devtools "devtools") package.

```r
install.packages("devtools")
library(devtools)
install_github("rcastelo/VariantFiltering")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **VariantFiltering**
please use the [Bioconductor support site](http://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **VariantFiltering**
please use the GitHub issues link at the top-right of this page
([https://github.com/rcastelo/VariantFiltering/issues](https://github.com/rcastelo/VariantFiltering/issues)).
