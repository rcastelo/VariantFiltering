\name{VariantFiltering-deprecated}
\alias{VariantFiltering-deprecated}
\alias{dbSNPpresent}
\alias{dbSNPpresent<-}
\alias{dbSNPpresent,VariantFilteringResults-method}
\alias{dbSNPpresent<-,VariantFilteringResults-method}
\alias{dbSNPpresent<-,VariantFilteringResults,ANY-method}
\alias{variantType}
\alias{variantType<-}
\alias{variantType,VariantFilteringResults-method}
\alias{variantType<-,VariantFilteringResults-method}
\alias{variantType<-,VariantFilteringResults,logical-method}
\alias{variantType<-,VariantFilteringResults,ANY,logical-method}
\alias{variantLocation}
\alias{variantLocation<-}
\alias{variantLocation,VariantFilteringResults-method}
\alias{variantLocation<-,VariantFilteringResults-method}
\alias{variantLocation<-,VariantFilteringResults,logical-method}
\alias{variantLocation<-,VariantFilteringResults,ANY,logical-method}
\alias{variantConsequence}
\alias{variantConsequence<-}
\alias{variantConsequence,VariantFilteringResults-method}
\alias{variantConsequence<-,VariantFilteringResults-method}
\alias{variantConsequence<-,VariantFilteringResults,logical-method}
\alias{variantConsequence<-,VariantFilteringResults,ANY,logical-method}
\alias{aaChangeType}
\alias{aaChangeType<-}
\alias{aaChangeType,VariantFilteringResults-method}
\alias{aaChangeType<-,VariantFilteringResults,character-method}
\alias{OMIMpresent}
\alias{OMIMpresent<-}
\alias{OMIMpresent,VariantFilteringResults-method}
\alias{OMIMpresent<-,VariantFilteringResults-method}
\alias{OMIMpresent<-,VariantFilteringResults,character-method}
\alias{naMAF}
\alias{naMAF<-}
\alias{naMAF,VariantFilteringResults-method}
\alias{naMAF<-,VariantFilteringResults,logical-method}
\alias{maxMAF}
\alias{maxMAF<-}
\alias{maxMAF,VariantFilteringResults-method}
\alias{maxMAF<-,VariantFilteringResults,numeric-method}
\alias{MAFpop,VariantFilteringResults-method}
\alias{MAFpop}
\alias{MAFpop<-}
\alias{MAFpop<-,VariantFilteringResults,logical-method}
\alias{MAFpop<-,VariantFilteringResults-method}
\alias{MAFpop<-,VariantFilteringResults,ANY,logical-method}
\alias{minPhastCons}
\alias{minPhastCons<-}
\alias{minPhastCons,VariantFilteringResults-method}
\alias{minPhastCons<-,VariantFilteringResults-method}
\alias{minPhastCons<-,VariantFilteringResults,ANY-method}
\alias{minPhylostratum}
\alias{minPhylostratum<-}
\alias{minPhylostratum,VariantFilteringResults-method}
\alias{minPhylostratum<-,VariantFilteringResults-method}
\alias{minPhylostratum<-,VariantFilteringResults,ANY-method}
\alias{minScore5ss}
\alias{minScore5ss<-}
\alias{minScore5ss,VariantFilteringResults-method}
\alias{minScore5ss<-,VariantFilteringResults,ANY-method}
\alias{minScore5ss<-,VariantFilteringResults-method}
\alias{minScore3ss}
\alias{minScore3ss<-}
\alias{minScore3ss,VariantFilteringResults-method}
\alias{minScore3ss<-,VariantFilteringResults,ANY-method}
\alias{minScore3ss<-,VariantFilteringResults-method}
\alias{minCUFC}
\alias{minCUFC<-}
\alias{minCUFC,VariantFilteringResults-method}
\alias{minCUFC<-,VariantFilteringResults,numeric-method}

\title{Deprecated functions in package \sQuote{VariantFiltering}}

\description{
  There are currently no deprecated functions in \sQuote{VariantFiltering}.
}
%%\description{
%%  These functions are provided for compatibility with older versions
%%  of \sQuote{VariantFiltering} only, and will be defunct at the next release.
%%}

\details{
  There are currently no deprecated functions.
%%  The following functions are deprecated and will be made defunct; use
%%  the replacement indicated below:
%%
%%  \itemize{
%%    \item{dbSNPpresent(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{OMIMpresent(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{variantType(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{aaChangeType(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{variantLocation(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{variantConsequence(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{naMAF(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{MAFpop(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{maxMAF(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{minPhastCons(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{minPhylostratum(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{minScore5ss(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{minScore3ss(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%    \item{minCUfc(): This functionality is now replaced by the combined used of functions \code{\link{filters}()}, \code{\link{active}()}, \code{\link{cutoffs}()}, code{\link{change}()}}
%%  }
}
\author{
  Robert Castelo <robert.castelo@upf.edu>
}
\keyword{internal}
\keyword{misc}
