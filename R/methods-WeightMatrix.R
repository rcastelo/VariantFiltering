readWm <- function(fname, locations=levels(.location()), strictLocations=FALSE) {
  return(new("WeightMatrix", wm=.Call("scoss_read_wm", fname),
             locations=locations, strictLocations=strictLocations))
}

setMethod("wmLocations", signature(x="WeightMatrix"),
          function(x) {
            x@locations
          })

setMethod("wmStrictLocations", signature(x="WeightMatrix"),
          function(x) {
            x@strictLocations
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
            cat(sprintf("  strict locations: %s\n", object@strictLocations))
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
