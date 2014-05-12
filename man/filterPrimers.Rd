\docType{methods}
\name{filterPrimers}
\alias{filterPrimers}
\title{filter primers based on thresholds}
\usage{
filterPrimers(primers, minConPrimer = 1, minCGonAmp = 2)
}
\arguments{
  \item{primers}{a list of primers returned by
  \code{designPrimers} or a data frame returned by
  \code{filterPrimers} function with as.data.frame=TRUE
  argument.}

  \item{minConPrimer}{minimum number of non-CpG Cs on the
  primers}

  \item{minCGonAmp}{minimum number of CpGs covered on the
  targeted region}
}
\description{
\code{filterPrimers} filters the primers based on covered
CGs on amplicons and number of non-CpG Cs on the primers
}
\examples{
data(bisPrimers)
          filt.primers=filterPrimers(bisPrimers,minConPrimer=1,minCGonAmp=1)
}
\author{
Based on Arnaud Krebs' function, modified by Altuna Akalin
}
\seealso{
\code{\link{primers2ranges}}, \code{\link{designPrimers}}
}

