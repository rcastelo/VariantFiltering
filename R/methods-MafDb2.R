## constructor
MafDb2 <- function(provider, provider_version, download_url,
                   download_date, reference_genome,
                   rsIDgpSNVs, rsIDSNVs, rsIDidxSNVs,
                   data_pkgname, data_dirpath) {
  data_cache <- new.env(hash=TRUE, parent=emptyenv())
  data_pops <- list.files(path=data_dirpath, pattern="AF$")
  data_serialized_objnames <- sub(".rds", "",
                                  list.files(path=data_pops, pattern="*.rds"))

  assign(data_pkgname, list(), envir=data_cache)

  new("MafDb2", provider=provider,
                provider_version=provider_version,
                download_url=download_url,
                download_date=download_date,
                reference_genome=reference_genome,
                data_pkgname=data_pkgname,
                data_dirpath=data_dirpath,
                data_serialized_objnames=data_serialized_objnames,
                data_pops=data_pops,
                .data_cache=data_cache)
}

## accessors
setMethod("provider", "MafDb2", function(x) x@provider)

setMethod("providerVersion", "MafDb2", function(x) x@provider_version)

setMethod("referenceGenome", "MafDb2", function(x) x@reference_genome)

setMethod("organism", "MafDb2", function(object) organism(referenceGenome(object)))

setMethod("seqinfo", "MafDb2", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "MafDb2", function(x) seqnames(referenceGenome(x)))

setMethod("seqlengths", "MafDb2", function(x) seqlengths(referenceGenome(x)))

setMethod("seqlevelsStyle", "MafDb2", function(x) seqlevelsStyle(referenceGenome(x)))

setMethod("populations", "MafDb2", function(x) x@data_pops)


.decodeRAW2AF2 <- function(x) {
  x <- as.numeric(x)
  mask0s <- x == 0 | x == 254    ## decode both, raw byte values 0 and 254 as 0
  maskNAs <- x == 255            ## decode raw byte value 255 as NAs
  z <- x[!maskNAs & !mask0s]
  mask <- z <= 100 
  z[mask] <- z[mask] / 100
  z[!mask] <- ((z[!mask]-100)%%10)/10^(floor((z[!mask]-100)/10)+3)
  x[!maskNAs & !mask0s] <- z
  x[mask0s] <- NA_real_
  x[maskNAs] <- NA_real_
  x
}


## adapted from VariantTools::extractCoverageForPositions()
.extractRawFromRleList <- function(rlelst, pos) {
  if (any(!unique(seqnames(pos)) %in% names(rlelst)))
    stop("Some sequence names from input positions are missing from rlelst")
  if (any(width(pos) > 1L))
    stop("Some ranges are of width > 1")
  seqlevels(pos) <- names(rlelst)
  ord <- order(seqnames(pos))
  ans <- raw(length(pos))
  ans[ord] <- unlist(mapply(function(v, p) {
    runValue(v)[findRun(p, v)]
  }, rlelst, split(start(pos), seqnames(pos)), SIMPLIFY=FALSE), use.names=FALSE)
  ans
}

## adapted from BSgenome:::.normarg_ranges()
.str2gr <- function(ranges) {
  if (class(ranges) == "character") {
    spl1 <- strsplit(ranges, ":", fixed=TRUE)
    if (!all(sapply(spl1, length) == 2L))
      stop("Genomic ranges given in a character string should have the format CHR:START[-END]\n")
    ranges <- sapply(spl1, function(ranges) {
                             ss <- strsplit(ranges[2], "-", fixed=TRUE)[[1]]
                             if (length(ss) == 1L)
                               ss <- c(ss, ss)
                             c(ranges[1], ss)
                           })
    ranges <- t(ranges)
    ranges <- GRanges(seqnames=ranges[, 1],
                      IRanges(start=as.integer(ranges[, 2]),
                              end=as.integer(ranges[, 3])))
  } else if (class(ranges) == "VRanges") {
    ranges <- as(ranges, "GRanges")
    mcols(ranges) <- NULL
  } else if (class(ranges) != "GRanges" && class(ranges) != "GPos")
    stop("argument 'ranges' must be either a GRanges object, a GPos object or a character string with the format CHR:START[-END]\n")

  ranges
}

setMethod("mafByOverlaps", signature="MafDb2",
          function(x, ranges, pop="AF", caching=TRUE) {
            ranges <- .str2gr(ranges)

            if (seqlevelsStyle(ranges)[1] != seqlevelsStyle(x)[1])
              seqlevelsStyle(ranges) <- seqlevelsStyle(x)[1]

            if (class(pop) != "character")
              stop("argument 'pop' must be a character vector")

            snames <- unique(as.character(runValue(seqnames(ranges))))
            if (any(!snames %in% seqnames(x)))
              stop("Sequence names %s in 'ranges' not present in MafDb2 object.",
                   paste(snames[!snames %in% seqnames(x)], collapse=", "))

            if (any(!pop %in% populations(x)))
              stop(sprintf("population %s must be one of %s\n", pop, paste(populations(x), collapse=", ")))

            maflist <- get(x@data_pkgname, envir=x@.data_cache)
            missingMask <- !pop %in% names(maflist)
            for (popname in pop[missingMask])
              maflist[[popname]] <- RleList(compress=FALSE)
            anyMissing <- any(missingMask)

            ans <- as.data.frame(matrix(NA_real_, nrow=length(ranges), ncol=length(pop),
                                        dimnames=list(NULL, pop)))
            for (popname in pop) {
              missingMask <- !snames %in% names(maflist[[popname]])
              anyMissing <- anyMissing || any(missingMask)
              for (sname in snames[missingMask]) {
                fname <- file.path(x@data_dirpath, popname,
                                   sprintf("snvs.%s.%s.%s.rds", providerVersion(referenceGenome(x)), popname, sname))
                if (file.exists(fname))
                  maflist[[popname]][[sname]] <- readRDS(fname)
                else {
                  warning(sprintf("No ExAC data for chromosome %s.", sname))
                  maflist[[popname]][[sname]] <- Rle(raw(0))
                }
              }
              ans[[popname]] <- .decodeRAW2AF2(.extractRawFromRleList(maflist[[popname]], ranges))
            }

            if (anyMissing && caching)
              assign(x@data_pkgname, maflist, envir=x@.data_cache)
            rm(maflist)

            ans
          })

setMethod("mafById", signature="MafDb2",
          function(x, ids, pop="AF", caching=TRUE) {
            if (class(ids) != "character")
              stop("argument 'ids' must be a character string vector.")

            if (!exists("rsIDSNVs", envir=x@.data_cache)) {
              message("Loading first time annotations of rs identifiers to variants, produced by data provider.")
              rsIDSNVs <- readRDS(file.path(x@data_dirpath, "rsIDSNVs.rds"))
              assign("rsIDSNVs", rsIDSNVs, envir=x@.data_cache)
            }

            rsIDSNVs <- get("rsIDSNVs", envir=x@.data_cache)
            if (is.character(rsIDSNVs))      ## old inefficient storage and retrieval
              mt <- match(ids, rsIDSNVs)
            else if (is.integer(rsIDSNVs)) { ## more efficient storage and retrieval
              idsint <- rep(NA_integer_, length(ids))
              rsMask <- regexpr("^rs", ids) == 1
              idsint[rsMask] <- as.integer(sub(pattern="^rs", replacement="", x=ids[rsMask]))
              mt <- rep(NA_integer_, length(idsint))
              if (any(!is.na(idsint))) {
                idsintnoNAs <- idsint[!is.na(idsint)]
                ord <- order(idsintnoNAs)                        ## order ids to speed up
                mtfi <- findInterval(idsintnoNAs[ord], rsIDSNVs) ## call to findInterval()
                mtfi[ord] <- mtfi                                ## put matches into original order
                mt[!is.na(idsint)] <- mtfi                       ## integrate matches into result
              }
              mt[mt == 0] <- 1
              if (any(!is.na(mt))) {
                maskNAs <- idsint != rsIDSNVs[mt]
                mt[maskNAs] <- NA
              }
              if (any(!is.na(mt))) {
                if (!exists("rsIDidxSNVs", envir=x@.data_cache)) {
                  rsIDidxSNVs <- readRDS(file.path(x@data_dirpath, "rsIDidxSNVs.rds"))
                  assign("rsIDidxSNVs", rsIDidxSNVs, envir=x@.data_cache)
                }
                rsIDidxSNVs <- get("rsIDidxSNVs", envir=x@.data_cache)
                mt <- rsIDidxSNVs[mt]
              }
            } else
              stop("internal object 'rsIDSNVs' of unknown class.")

            if (any(!pop %in% populations(x)))
              stop(sprintf("population %s must be one of %s\n", pop, paste(populations(x), collapse=", ")))

            ans <- as.data.frame(matrix(NA_real_, nrow=length(ids), ncol=length(pop),
                                        dimnames=list(NULL, pop)))
            ans <- cbind(ID=ids, ans, stringsAsFactors=FALSE)

            if (any(!is.na(mt))) {
              if (!exists("rsIDgpSNVs", envir=x@.data_cache)) {
                rsIDgpSNVs <- readRDS(file.path(x@data_dirpath, "rsIDgpSNVs.rds"))
                assign("rsIDgpSNVs", rsIDgpSNVs, envir=x@.data_cache)
              }
              rsIDgpSNVs <- get("rsIDgpSNVs", envir=x@.data_cache)
              rng <- rsIDgpSNVs[mt[!is.na(mt)]]
              ans[!is.na(mt), pop] <- mafByOverlaps(x, rsIDgpSNVs[mt[!is.na(mt)]], pop, caching)
            }

            ans
          })

## show method
setMethod("show", "MafDb2",
          function(object) {
            cat(class(object), " object for ", organism(object), " (",
                provider(object), ")\n", sep="")
          })

## $ method
setMethod("$", signature(x="MafDb2"),
          function(x, name) {
            switch(name,
                   tag=gsub("MafDb.|.hs37d5", "", x@data_pkgname), ## by now just hardcoded, this should be part of the annotation package
                   stop("uknown MafDb2 slot.")
                   )
          })
