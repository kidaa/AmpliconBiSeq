\docType{methods}
\name{getSampleNames}
\alias{getSampleNames}
\alias{getSampleNames,AmpliconViews}
\alias{getSampleNames,AmpliconViews-method}
\title{get sample Names from AmpliconViews object}
\usage{
getSampleNames(x)

\S4method{getSampleNames}{AmpliconViews}(x)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}}}
}
\value{
character vector for sample Names
}
\description{
The function returns the sample Names stored in any of the
\code{\link{AmpliconViews}}
}
\examples{
data(ampViewEx);getSampleNames(ampViewEx)
}

