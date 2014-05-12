\docType{methods}
\name{ampData}
\alias{ampData}
\alias{ampData,AmpliconViews}
\alias{ampData,AmpliconViews-method}
\title{get data from AmpliconViews object}
\usage{
ampData(x)

\S4method{ampData}{AmpliconViews}(x)
}
\arguments{
  \item{x}{an \code{\link{AmpliconViews}}}
}
\value{
list of amplicon data
}
\description{
The function returns data stored in any of the
\code{\link{AmpliconViews}}
}
\examples{
data(ampViewEx);ampData(ampViewEx)
}

