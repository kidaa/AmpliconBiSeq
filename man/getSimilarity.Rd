\docType{methods}
\name{getSimilarity}
\alias{getSimilarity}
\alias{getSimilarity,AmpliconViews,character,character}
\alias{getSimilarity,AmpliconViews,character,character-method}
\title{get similarity for base methylation profiles in an amplicon}
\usage{
getSimilarity(x,sampleName,ampliconName,method="basic")

\S4method{getSimilarity}{AmpliconViews,character,character}(x, sampleName,
  ampliconName, method)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}} object}

  \item{sampleName}{name of the sample}

  \item{ampliconName}{name of the amplicon}
}
\description{
The function returns different metrics for methylation
profile similarity for bases in an amplicon
}

