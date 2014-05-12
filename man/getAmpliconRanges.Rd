\docType{methods}
\name{getAmpliconRanges}
\alias{getAmpliconRanges}
\alias{getAmpliconRanges,AmpliconViews}
\alias{getAmpliconRanges,AmpliconViews-method}
\title{get locations of amplicons from AmpliconViews object}
\usage{
getAmpliconRanges(x)

\S4method{getAmpliconRanges}{AmpliconViews}(x)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}}}
}
\value{
character vector for sample Names
}
\description{
The function returns amplicon locations of stored in any of
the \code{\link{AmpliconViews}}
}
\examples{
data(ampViewEx);getAmpliconRanges(ampViewEx)
}

