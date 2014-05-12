\docType{methods}
\name{designPrimers}
\alias{designPrimers}
\title{Design primers for amplicons}
\usage{
designPrimers(target, genome, tag = NULL, primer3 = NULL, settings = NULL,
  thermoParam = NULL, sizeRange = c(151, 500), Tm = c(55, 57, 58),
  ncores = 4)
}
\arguments{
  \item{target}{\code{DNAStringSet} or \code{GRanges}
  object. If it is a GRanges object you need to provide
  BSgenome package with it. If it is a \code{DNAStringSet}
  object then names argument should contain
  chromomsome,start and end, by following "chr_start_end"
  naming convetion.}

  \item{genome}{default is NULL, it needs to be set only if
  target is set to be a \code{GRanges} object}

  \item{tag}{extra information about amplicons, this will
  be concatanated to the names of BioStringSet object which
  is either created internally or provided by the user by
  supplying 'target' argument with a \code{BioStringSet}
  object}

  \item{primer3}{the filesystem path to primer3 executable,
  default value NULL will use the version installed with
  the package}

  \item{settings}{text file for p3 settings. Default value
  NULL will use the default settings file installed with
  the package.}

  \item{termoParam}{Location for thermodynamic parameters
  for primer3, this should designate a directory, not a
  single file. Default value NULL, will use the directory
  that is installed with the package.}

  \item{sizeRange}{a two element vector of integers for
  size range of the amplicons}

  \item{Tm}{a three element vector having minimum,optimum
  and maximum melting tempratures for primer desing (in
  that order).}

  \item{ncores}{number of cores to run on. default=1}
}
\description{
Designs primers for Amplicon Bisulfite sequencing
experiments
}
\author{
Based on Arnaud Krebs' function, modified by Altuna Akalin
}
\seealso{
\code{\link{primers2ranges}}, \code{\link{filterPrimers}}
}

