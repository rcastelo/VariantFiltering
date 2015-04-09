.defaultFilters <- list(dbSNP=function(x) !is.na(allVariants(x, groupBy="nothing")$dbSNP),
                        OMIM=function(x) !is.na(allVariants(x, groupBy="nothing")$OMIM),
                        variantType=function(x) allVariants(x, groupBy="nothing")$TYPE %in% names(cutoffs(x)$variantType)[cutoffs(x)$variantType],
                        aaChangeType=function(x) {
                          mask <- rep(TRUE, length(x))
                          if (cutoffs(x)$aaChangeType %in% c("Conservative", "Radical")) {
                            aachangetype <- allVariants(x, groupBy="nothing")$AAchangeType
                            mask <- is.na(aachangetype) | achangetype %in% cutoffs(x)$aaChangeType
                          }
                          mask
                        },
                        SOterms=function(x) {
                          soterms <- cutoffs(x)$SOterms
                          if (is.null(soterms) || all(is.na(soterms))) ## no SO terms specified imply no restriction
                            soterms <- nodes(sog(x))[sapply(nodeData(sog(x), nodes(sog(x)), "varIdx"), length) > 0]
                          else {
                            soterms <- .findSOIDs(cutoffs(x)$SOterms, sog(x))

                            ## restrict to the given SO terms and all their more specific ancestors with 
                            ## at least one annotated variant
                            soterms <- unique(c(soterms, .ancestorsSO(param(x), soterms)))
                            soterms <- soterms[sapply(nodeData(sog(x), soterms, "varIdx"), length) > 0]
                          }

                          mask <- rep(TRUE, length(x))
                          if (length(soterms) > 0) {
                            varidxinsog <- unique(unlist(nodeData(sog(x), soterms, "varIdx"), use.names=FALSE))
                            mask <- rep(FALSE, length(x))
                            mask[varidxinsog] <- TRUE
                          }

                          mask
                        }
                        )

.defaultCutoffs <- list(dbSNP=NA,
                        OMIM=NA,
                        variantType=c(SNV=TRUE, Insertion=TRUE, Deletion=TRUE, MNV=TRUE, Delins=TRUE),
                        aaChangeType="Any",
                        SOterms=NA
                        )
