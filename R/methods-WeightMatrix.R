readWm <- function(fname, locations=levels(.location()), strictLocations=FALSE) {

  if (any(!locations %in% levels(.location())))
    stop(sprintf("values in argument 'locations' must be one of: %s\n",
                 paste(levels(.location()), collapse=", ")))

  if (length(strictLocations) > 1 && length(strictLocations != length(locations)))
    stop("argument 'strictLocations' must be either one truth value (TRUE or FALSE) for all locations, or as many truth values as values given in the 'locations' argument\n")

  if (length(strictLocations) == 1)
    strictLocations <- rep(strictLocations, times=length(locations))

  if (is.null(names(strictLocations)))
    names(strictLocations) <- locations

  if (any(!names(strictLocations) %in% locations))
    stop("value names in argument 'strictLocations' must match one of the values in argument 'locations'\n")

  strictLocations <- strictLocations[locations]

  return(new("WeightMatrix", wm=.Call("scoss_read_wm", fname),
             locations=locations, strictLocations=strictLocations))
}

setMethod("wmLocations", signature(x="WeightMatrix"),
          function(x) {
            x@locations
          })

setReplaceMethod("wmLocations", signature(x="WeightMatrix"),
                 function(x, value) {
                   if (class(value) != "character")
                     stop("the assigned value must be a character vector\n")

                   if (any(!value %in% levels(.location())))
                     stop(sprintf("values must be one of: %s\n",
                                  paste(levels(.location()), collapse=", ")))

                   sl <- wmStrictLocations(x)
                   mask <- value %in% names(sl)
                   slxtra <- logical(0)
                   if (any(!mask)) {
                     slxtra <- rep(FALSE, sum(!mask))
                     names(slxtra) <- value[!mask]
                   }
                   sl <- c(sl, slxtra)
                   sl <- sl[value]

                   x@locations <- value
                   x@strictLocations <- sl

                   x
                 })

setMethod("wmStrictLocations", signature(x="WeightMatrix"),
          function(x) {
            x@strictLocations
          })

setReplaceMethod("wmStrictLocations", signature(x="WeightMatrix"),
                 function(x, value) {
                   if (class(value) != "logical")
                     stop("the assigned value must be a logical vector\n")

                   locations <- wmLocations(x)

                   if ((length(value) > 1 && length(value) > length(locations)))
                     stop("more values than locations in wmLocations(x)\n")

                   if ((length(value) > 1 && length(value) < length(locations)))
                     stop("when assigning > 1 fewer values than locations in wmLocations(x), values should have names\n")

                   if (length(value) == 1 && is.null(names(value)))
                     value <- rep(value, times=length(locations))

                   if (is.null(names(value)))
                     names(value) <- locations

                   if (any(!names(value) %in% locations))
                     stop("value names must match one of the values returned by wmLocations(x)\n")

                   x@strictLocations[locations[locations %in% names(value)]] <-
                     value[locations[locations %in% names(value)]]

                   x
                 })

setMethod("width", signature(x="WeightMatrix"),
          function(x) {
            .Call("scoss_width_wm", x@wm)
          })

setMethod("conservedPositions", signature(x="WeightMatrix"),
          function(x) {
            .Call("scoss_conserved_positions_wm", x@wm)
          })

setMethod("wmName", signature(x="WeightMatrix"),
          function(x) {
            .Call("scoss_name_wm", x@wm)
          })

setMethod("show", signature(object = "WeightMatrix"),
          function(object) {
            .Call("scoss_show_wm", object@wm)
            cat(sprintf("  locations: %s\n", paste(object@locations, collapse=", ")))
            cat(sprintf("  strict locations: %s\n", paste(object@strictLocations, collapse=", ")))
          })


setMethod("wmScore", signature(object="WeightMatrix", dnaseqs="character"),
          function(object, dnaseqs) {
            nsites <- nchar(dnaseqs) - width(object) + 1L

            if (any(nsites < 1))
              stop(sprintf("strings %s have fewer characters than the width of the input weight matrix (%d)",
                           paste(head(which(nsites < 1)), sep=", "), width(object)))

            return(.Call("scoss_wm_score", object@wm, dnaseqs, sum(nsites)))
          })

setMethod("wmScore", signature(object="WeightMatrix", dnaseqs="DNAStringSet"),
          function(object, dnaseqs) {
            nsites <- elementLengths(dnaseqs) - width(object) + 1L

            if (any(nsites < 1))
              stop(sprintf("DNA strings %s have fewer nucleotides than the width of the input weight matrix (%d)",
                           paste(head(which(nsites < 1)), sep=", "), width(object)))

            return(.Call("scoss_wm_score_DNAStringSet", object@wm, dnaseqs, sum(nsites)))
          })
