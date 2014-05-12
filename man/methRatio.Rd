\docType{methods}
\name{methRatio}
\alias{methRatio}
\title{convert ampliconView objects to a table of methylation ratio}
\usage{
methRatio(x,per.region=FALSE,asGRanges=TRUE,dup.resolve=TRUE,simple.tags=TRUE,coverage.th=0)
}
\arguments{
  \item{x}{an AmpliconViews object}

  \item{per.region}{logical, if FALSE base-pair methylation
  ratio (def:FALSE)}

  \item{asGRanges}{logical, if TRUE object returned is a
  GRanges object (def:TRUE)}

  \item{dup.resolve}{logical, if TRUE , and if
  per.region=FALSE, the duplicated CpGs will be resolved by
  taking the one with highest coverage. Since amplicon
  designs can overlap, there might be CpGs that are covered
  by different amplicon designs at the same time, this will
  result in duplicated CpGs in the resulting object.
  (def:TRUE)}

  \item{simple.tags}{logical. If TRUE, tags associated with
  amplicons will be output as a single column. If FALSE, a
  tag column will be output for each experiment.}

  \item{coverage.th}{default to 0. If this is set to a
  value larger than zero, the bases/regions that have have
  coverage below this value will have NA methylation and NA
  coverage.}
}
\value{
data frame or GRanges object with locations of CpGs and
methylation ratio and coverage
}
\description{
convert ampliconView objects to a table of methylation
ratio
}
\examples{
# methRatio(a.list,per.region=FALSE,asGRanges=TRUE,dup.resolve=TRUE)
}

