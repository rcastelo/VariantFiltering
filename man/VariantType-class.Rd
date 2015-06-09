\name{VariantType-class}
\docType{class}

% Subclasses
\alias{FiveSpliceSiteVariants-class}
\alias{ThreeSpliceSiteVariants-class}

% Constructors for subclasses :
\alias{FiveSpliceSiteVariants}
\alias{ThreeSpliceSiteVariants}

% accessors :
\alias{minIntronLength}
\alias{minIntronLength,FiveSpliceSiteVariants-method}
\alias{minIntronLength,ThreeSpliceSiteVariants-method}
\alias{minIntronLength<-}
\alias{minIntronLength<-,FiveSpliceSiteVariants-method}
\alias{minIntronLength<-,ThreeSpliceSiteVariants-method}
\alias{upstream}
\alias{upstream,FiveSpliceSiteVariants-method}
\alias{upstream,ThreeSpliceSiteVariants-method}
\alias{upstream<-}
\alias{upstream<-,FiveSpliceSiteVariants-method}
\alias{upstream<-,ThreeSpliceSiteVariants-method}
\alias{downstream}
\alias{downstream,FiveSpliceSiteVariants-method}
\alias{downstream,ThreeSpliceSiteVariants-method}
\alias{downstream<-}
\alias{downstream<-,FiveSpliceSiteVariants-method}
\alias{downstream<-,ThreeSpliceSiteVariants-method}

% show method:
\alias{show,FiveSpliceSiteVariants-method}
\alias{show,ThreeSpliceSiteVariants-method}

\title{FiveSpliceSiteVariants and ThreeSpliceSiteVariants subclasses of VariantType}

\description{
  \code{VariantType} subclasses specify the type of variant to be located with
  the function \code{\link[VariantAnnotation]{locateVariants}()} from the
  \code{VariantAnnotation} package. Here we describe two such subclasses defined
  within the \code{VariantFiltering} package: \code{FiveSpliceSiteVariants} and
  \code{ThreeSpliceSiteVariants}. For further details on the
  \code{\link[VariantAnnotation]{VariantType-class}} please consult the corresponding
  manual page in the \code{VariantAnnotation} package.
}

\usage{
    FiveSpliceSiteVariants(minIntronLength=20L, upstream = 3L, downstream = 4L)
    ThreeSpliceSiteVariants(minIntronLength=20L, upstream = 18L, downstream = 3L)
}

\details{
  The subclasses \code{FiveSpliceSiteVariants} and \code{ThreeSpliceSiteVariants}
  can be used as in the \code{region} argument to the
  \code{\link[VariantAnnotation]{locateVariants}()} function from the \code{VariantAnnotation}
  package or as values in the \code{regionAnnotations} argument to the
  \code{\link{VariantFilteringParam}()} function that constructs a
  \code{\link{VariantFilteringParam-class}} object. These two subclasses allow
  the package \code{VariantFiltering} and the function \code{\link[VariantAnnotation]{locateVariants}()}
  in the \code{VariantAnnotation} package, to annotate variants located within
  a splice site region defined by the user. The boundaries of this region as well
  as the minimum intron length can be defined with the following arguments to the
  constructor functions:
}

\arguments{
  \item{minIntronLength}{
    Minimum intron length in nucleotides, by default \code{minIntronLength=20L}.
    Splice sites annotated within introns smaller than this length will not be
    used to annotate variants.
  }
  \item{upstream, downstream}{
    Single \code{integer} values representing the number of nucleotides
    upstream and downstream of the 5' and 3' dinucleotide ends of introns.
    By default, the constructor function \code{FiveSpliceSiteVariants} has
    values \code{upstream=3L} and \code{downstream=4L}, corresponding to
    3 nucleotides upstream and 4 nucleotides downstream of the first dinucleotide
    of the intron, while \code{ThreeSpliceSiteVariants} has values
    \code{upstream=18L} and \code{downstream=3L}, corresponding to 18 nucleotides
    upstream and 4 nucleotides downstream of the last dinucleotide of the intron.
  }
}

\section{Accessors}{
  In the following code, \code{x} is a \code{PromoterVariants} or a
  \code{AllVariants} object.
  \describe{
    \item{}{
      \code{minIntronLength(x)}, \code{minIntronLength(x) <- value}:
      Gets or sets the minimum number of nucleotides long of introns
      to be considered in the annotation of splice sites.
    }
    \item{}{
      \code{upstream(x)}, \code{upstream(x) <- value}:
      Gets or sets the number of nucleotides upstream of the 5' and
      3' end of the intron, excluding the first and last dinucleotide,
      respectively.
    }
    \item{}{
      \code{downstream(x)}, \code{downstream(x) <- value}:
      Gets or sets the number of nucleotides downstream of the 5' and
      3' end of the intron, excluding the first and last dinucleotide,
      respectively.
    }
  }
}

\author{Robert Castelo <robert.castelo@upf.edu>}

\seealso{
  \itemize{
    \item The \link[VariantAnnotation]{VariantType-class} man page in the 
          \code{VariantAnnotation} package.
  }
}

\examples{
  FiveSpliceSiteVariants()
  ThreeSpliceSiteVariants()
}
