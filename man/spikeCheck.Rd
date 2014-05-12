\docType{methods}
\name{spikeCheck}
\alias{spikeCheck}
\title{check the quality of the spike-ins}
\usage{
spikeCheck(proj, auxName = "all", sampleName = NULL, coverage = 0,
  targets = NULL, ...)
}
\arguments{
  \item{proj}{\code{\link[QuasR]{qProject}} object from
  QuasR preferably produced by ampBiSeqAlign}

  \item{auxName}{string for which spike-ins should be
  checked, if equals "all" everything is checked}

  \item{sampleName}{string for which sample the spike-ins
  should be checked}

  \item{coverage}{minumum coverage before calculation of %
  methylation}

  \item{targets}{named list of GRanges objects identifies
  targeted regions in given auxiliary genomes}

  \item{...}{arguments to be passed to \code{\link{hist}}
  function}
}
\value{
a list of percentages of methylated CpGs in all reads.
}
\description{
This function provides methylation statistics for
spike-ins. It plots a histogram or set of histograms from
spike-in experiments which are helpful to deduce conversion
efficiency of the experiment.
}
\examples{
# spikeCheck(proj,auxName="all",sampleName="ES.2")
          # spikeCheck(proj,auxName="T7",sampleName="ES.2")

          #la.gr=readRDS("/work2/gschub/Juliane/HTS/QuasRBiSeq/lamdaGR.rds")
          # par(mfrow=c(2,2))
          #spikeCheck(proj,auxName="all",sampleName="ES.2",targets=list(lambda=la.gr))
}

