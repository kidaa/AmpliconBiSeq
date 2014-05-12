\docType{methods}
\name{getAmplicon}
\alias{getAmplicon}
\alias{getAmplicon,AmpliconViews,character,character}
\alias{getAmplicon,AmpliconViews,character,character-method}
\title{subset AmpliconViews object}
\usage{
getAmplicon(x,sampleNames=NULL,ampliconNames=NULL)

\S4method{getAmplicon}{AmpliconViews,character,character}(x, sampleNames,
  ampliconNames)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}} object}

  \item{sampleNames}{character vector of sample names}

  \item{ampliconNames}{character vector of amplicon names}
}
\value{
\code{\link{AmpliconViews}} object
}
\description{
The function returns data stored in any of the
\code{\link{AmpliconViews}}
}
\examples{
data(ampViewEx)
          getAmplicon(ampViewEx,"mock4","chr18_69674375_69674775")
}

