\docType{methods}
\name{getAvMeth}
\alias{getAvMeth}
\alias{getAvMeth,AmpliconViews,character,character}
\alias{getAvMeth,AmpliconViews,character,character-method}
\title{get average methylation for bases in an amplicon}
\usage{
getAvMeth(x,sampleName,ampliconName)

\S4method{getAvMeth}{AmpliconViews,character,character}(x, sampleName,
  ampliconName)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}} object}

  \item{sampleName}{name of the sample}

  \item{ampliconName}{name of the amplicon}
}
\description{
The function returns average methylation for each base for
a given amplicon in a given sample.
}

