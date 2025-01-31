# VariantFiltering: Filtering of coding and non-coding genetic variants

[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/VariantFiltering.svg)](http://bioconductor.org/packages/release/bioc/html/VariantFiltering.html "How long has been VariantFiltering in a release of Bioconductor")
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/release/VariantFiltering.svg)](http://bioconductor.org/packages/stats/bioc/VariantFiltering/ "Ranking by the number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months.")
[![Support posts](http://bioconductor.org/shields/posts/VariantFiltering.svg)](https://support.bioconductor.org/t/VariantFiltering/ "Support site activity on VariantFiltering, last 6 months: answered posts/total posts.")
[![R-CMD-check-bioc](https://github.com/rcastelo/VariantFiltering/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/rcastelo/VariantFiltering/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/rcastelo/VariantFiltering/coverage.svg?branch=devel)](https://app.codecov.io/github/rcastelo/VariantFiltering?branch=devel)


**Current build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/VariantFiltering.svg)](https://bioconductor.org/packages/release/bioc/html/VariantFiltering.html#archives "Whether VariantFiltering release is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/VariantFiltering.svg)](https://bioconductor.org/packages/release/bioc/html/VariantFiltering.html#since "Number of recursive dependencies needed to install package") 
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/VariantFiltering.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/VariantFiltering "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/VariantFiltering.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/VariantFiltering/ "Bioconductor release build")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/VariantFiltering.svg)](https://bioconductor.org/packages/devel/bioc/html/VariantFiltering.html#archives "Whether VariantFiltering devel is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/VariantFiltering.svg)](https://bioconductor.org/packages/devel/bioc/html/VariantFiltering.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/devel/bioc/VariantFiltering.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/VariantFiltering "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/VariantFiltering.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/VariantFiltering/ "Bioconductor devel build")

The `VariantFiltering` package aids at filtering genetic variants using
different criteria such as inheritance model, amino acid change consequence,
minor allele frequency across human populations, splice site strength,
conservation, etc.

## Installation

You can install the `VariantFiltering` package from this GitHub repo using
the [remotes](https://cran.r-project.org/package=remotes) package, as follows:

```
library(remotes)

install_github("rcastelo/VariantFiltering")
```

Provided that you have installed first all its Bioconductor dependencies;
see the `DESCRIPTION` file. The vignette contains an example on how to use it.

## Questions, bug reports and issues

For bug reports and feature requests, please use the GitHub issues
[tab](https://github.com/rcastelo/VariantFiltering/issues)
at the top-left of this page.

## Contributing

Contributions to the software codebase of `VariantFiltering` are welcome as
long as contributors abide to the terms of the
[Bioconductor Contributor Code of Conduct](https://bioconductor.org/about/code-of-conduct).
If you want to contribute to the development of VariantFiltering please open an
[issue](https://github.com/functionalgenomics/VariantFiltering/issues) to start discussing
your suggestion or, in case of a bugfix or a straightforward feature, directly a
[pull request](https://github.com/functionalgenomics/VariantFiltering/pulls).

## Funding

This software project is supported in part by the Spanish
[Ministry of Science, Innovation and Universities](https://www.ciencia.gob.es).

