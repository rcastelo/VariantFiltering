.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  data_dir <- system.file("extdata", package=pkgname, lib.loc=libname)

  ## read human genes phylostrata
  assign("humanGenesPhylostrata", GenePhylostrataDb("human",
                                                    "http://genomebiology.com/content/supplementary/1471-2164-14-117-s1.xlsx",
                                                    "Nov 21, 2013",
                                                    "Neme R and Tautz D, BMC Genomics, 14:117, 2013",
                                                    pkgname, data_dir),
         envir=ns)

  ## read Sequence Ontology graph for the 'sequence_variant' ontology
  assign("sequence_variant.gSOXP", readRDS(file.path(data_dir, "sequence_variant.gSOXPmarch2015.rds")),
         envir=ns)

  namespaceExport(ns, "humanGenesPhylostrata")
  namespaceExport(ns, "sequence_variant.gSOXP")
}
