\docType{methods}
\name{getMethMat}
\alias{getMethMat}
\alias{getMethMat,AmpliconViews,character,character}
\alias{getMethMat,AmpliconViews,character,character-method}
\title{get methylation matrix}
\usage{
getMethMat(x,sampleName,ampliconName)

\S4method{getMethMat}{AmpliconViews,character,character}(x, sampleName,
  ampliconName)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}} object}

  \item{sampleName}{name of the sample}

  \item{ampliconName}{name of the amplicon}
}
\description{
If \code{AmpliconViews} object has the methylation call
matrix this function returns it.
}

