browseVariants <- function(vfResultsObj) {
  appDir <- system.file("shinyApp", package="VariantFiltering")
  if (appDir == "")
    stop("The VariantFiltering Shiny app cannot be found within the package.")

  argumentMissing <- TRUE
  if (!missing(vfResultsObj)) {
    if (class(vfResultsObj) != "VariantFilteringResults")
      stop("The input argument to 'browseVariants()' must be a 'VariantFilteringResults' object.")
    argumentMissing <- FALSE
  }

  runApp(appDir)
}
