\docType{methods}
\name{primers2ranges}
\alias{primers2ranges}
\title{convert primer list to genomic intervals}
\usage{
primers2ranges(primers, as.data.frame = FALSE)
}
\arguments{
  \item{primers}{a list of primers returned by
  \code{designPrimers} function}

  \item{as.data.frame}{logical indicating if a data frame
  should be returned instead of \code{GRanges} object.}
}
\description{
function returns a GRanges or data frame object from a list
of primers designed by \code{designPrimers} function after
calculating genomic location of the amplicon targeted by
the primers.
}
\examples{
data(bisPrimers)
 gr.pr=primers2ranges(bisPrimers)
}
\seealso{
\code{\link{filterPrimers}}, \code{\link{designPrimers}}
}

