CutoffsList <- function(...) {
  arg <- list(...)
  if (length(arg) == 1L && extends(class(arg[[1L]]), "list"))
    arg <- arg[[1L]]
  new("CutoffsList", listData=arg)
}

setGeneric("change<-", function(x, value, ...) standardGeneric("change<-"))

.checkValueCutoffArgs <- function(x, value, cutoff) {
  if (missing(cutoff))
    stop("please specify a 'cutoff' name as a second argument.")

  if (!is.character(cutoff))
    stop("the 'cutoff' argument must be a character string.")

  if (!cutoff %in% names(x))
    stop(sprintf("'%s' does not form part of the available cutoffs.", cutoff))
}

setReplaceMethod("change", signature(x="CutoffsList", value="integer"),
                 function(x, value, cutoff) {
                   change(x, cutoff) <- as.numeric(value)
                 })

setReplaceMethod("change", signature(x="CutoffsList", value="numeric"),
                 function(x, value, cutoff) {
                   .checkValueCutoffArgs(x, value, cutoff)

                   if (!is.numeric(x[[cutoff]]) && !is.integer(x[[cutoff]]))
                       stop("this cutoff does not take numeric values.")

                   x@listData[[cutoff]][1] <- value
                   x
                 })

setReplaceMethod("change", signature(x="CutoffsList", value="logical"),
                 function(x, value, cutoff) {
                   .checkValueCutoffArgs(x, value, cutoff)

                   if (!is.logical(x[[cutoff]]))
                     stop("this cutoff does not take a logical values.")

                   mtvalue <- integer(0)
                   if (is.null(names(value))) {
                     if (length(value) > 1)
                       stop("multiple values must have names.")
                     if (is.null(names(x[[cutoff]])))
                       mtvalue <- 1:length(value)
                     else {
                       value <- do.call("names<-", list(rep(value, length(x[[cutoff]])), names(x[[cutoff]])))
                       mtvalue <- match(names(value), names(x[[cutoff]]))
                     }
                   } else {
                     mask <- !names(value) %in% names(x[[cutoff]])
                     if (any(mask))
                       stop(sprintf("value names %s do not exist for cutoff '%s'",
                                    paste(names(value)[mask], collapse=", "), cutoff))
                     mtvalue <- match(names(value), names(x[[cutoff]]))
                   }

                   x@listData[[cutoff]][mtvalue] <- value
                   x
                 })
  
setReplaceMethod("change", signature(x="CutoffsList", value="character"),
                 function(x, value, cutoff) {
                   .checkValueCutoffArgs(x, value, cutoff)

                   if (!is.character(x[[cutoff]]) && !is.factor(x[[cutoff]]))
                     stop("this cutoff does not take character string values.")

                   if (length(value) > 1)
                     stop("a string value for a cutoff must be a singleton.")

                   curvalues <- x[[cutoff]]
                   if (is.factor(curvalues)) {
                     if (!value %in% levels(curvalues))
                       stop(sprintf("invalid value for cutoff '%s'", cutoff))
                   }

                   x@listData[[cutoff]][1] <- value
                   x
                 })

setReplaceMethod("cutoffs", signature(x="VariantFilteringResults", value="CutoffsList"),
                 function(x, value) {
                   x@cutoffs <- value
                   x
                 })

setReplaceMethod("sortings", signature(x="VariantFilteringResults", value="CutoffsList"),
                 function(x, value) {
                   x@sortings <- value
                   x
                 })
