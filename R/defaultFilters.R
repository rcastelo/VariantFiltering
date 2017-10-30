.defaultFilters <- list(dbSNP=function(x) !is.na(VariantFiltering::allVariants(x, groupBy="nothing")$dbSNP),
                        OMIM=function(x) !is.na(VariantFiltering::allVariants(x, groupBy="nothing")$OMIM),
                        variantType=function(x) VariantFiltering::allVariants(x, groupBy="nothing")$TYPE %in% names(VariantFiltering::cutoffs(x)$variantType)[VariantFiltering::cutoffs(x)$variantType],
                        aaChangeType=function(x) {
                          mask <- rep(TRUE, length(x))
                          if (VariantFiltering::cutoffs(x)$aaChangeType %in% c("Conservative", "Radical")) {
                            aachangetype <- VariantFiltering::allVariants(x, groupBy="nothing")$AAchangeType
                            mask <- is.na(aachangetype) | achangetype %in% VariantFiltering::cutoffs(x)$aaChangeType
                          }
                          mask
                        },
                        SOterms=function(x) {
                          soterms <- VariantFiltering::cutoffs(x)$SOterms
                          if (is.null(soterms) || all(is.na(soterms))) ## no SO terms specified imply no restriction
                            soterms <- graph::nodes(VariantFiltering::sog(x))[sapply(graph::nodeData(VariantFiltering::sog(x), graph::nodes(VariantFiltering::sog(x)), "varIdx"), length) > 0]
                          else {
                            soterms <- VariantFiltering:::.findSOIDs(VariantFiltering::cutoffs(x)$SOterms, VariantFiltering::sog(x))

                            ## restrict to the given SO terms and all their more specific ancestors with 
                            ## at least one annotated variant
                            soterms <- unique(c(soterms, VariantFiltering:::.ancestorsSO(VariantFiltering::param(x), soterms)))
                            soterms <- soterms[sapply(graph::nodeData(VariantFiltering::sog(x), soterms, "varIdx"), length) > 0]
                          }

                          mask <- rep(TRUE, length(x))
                          if (length(soterms) > 0) {
                            varidxinsog <- unique(unlist(graph::nodeData(VariantFiltering::sog(x), soterms, "varIdx"), use.names=FALSE))
                            mask <- rep(FALSE, length(x))
                            mask[varidxinsog] <- TRUE
                          }

                          mask
                        }
                        )

.defaultFilterDescriptions <- DataFrame(Description=c("Presence in dbSNP",
                                                      "Presence in OMIM",
                                                      "Type of variant (SVN, Insertion, Deletion, MNV, Delins)",
                                                      "Type of amino acid change (conservative, radical)",
                                                      "Annotated to given Sequence Ontology terms"),
                                        row.names=c("dbSNP", "OMIM", "variantType", "aaChangetype", "SOterms"))

.defaultCutoffs <- list(variantType=c(SNV=TRUE, Insertion=TRUE, Deletion=TRUE, MNV=TRUE, Delins=TRUE),
                        aaChangeType=c("Any", "Conservative", "Radical"),
                        SOterms="Any"
                        )
