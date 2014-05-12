\docType{methods}
\name{getMetaMethylation}
\alias{getMetaMethylation}
\alias{getMetaMethylation,AmpliconViews,character,character-method}
\alias{getMetaMethylation,AmpliconViews,character,missing}
\title{get meta-methylation info from AmpliconView object}
\usage{
getMetaMethylation(x,sampleName,ampliconName)

\S4method{getMetaMethylation}{AmpliconViews,character,character}(x, sampleName,
  ampliconName)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}} object}

  \item{sampleName}{name of the sample}

  \item{ampliconName}{name of the amplicon}
}
\value{
a list containing meta-methylation information: patterns,
explained percentage of data SVD(PCA) results
}
\description{
get meta-methylation info from AmpliconView object
}

