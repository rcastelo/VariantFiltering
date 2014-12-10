.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  data_dir <- system.file("extdata", package=pkgname, lib.loc=libname)

  ## read idlbn-weight matrices for human donor and acceptor sites
  ## wmDonorSites <- readWm(file.path(data_dir, "hsap.donors.hcmc10_15_1.ibn"))
  ## wmAcceptorSites <- readWm(file.path(data_dir, "hsap.acceptors.hcmc10_15_1.ibn"))

  ## assign("wmDonorSites", wmDonorSites, envir=ns)
  ## assign("wmAcceptorSites", wmAcceptorSites, envir=ns)

  assign("humanGenesPhylostrata", GenePhylostrataDb("human",
                                                    "http://genomebiology.com/content/supplementary/1471-2164-14-117-s1.xlsx",
                                                    "Nov 21, 2013",
                                                    "Neme R and Tautz D, BMC Genomics, 14:117, 2013",
                                                    pkgname, data_dir),
         envir=ns)

  ## namespaceExport(ns, "wmDonorSites")
  ## namespaceExport(ns, "wmAcceptorSites")
  namespaceExport(ns, "humanGenesPhylostrata")
}
