\docType{class}
\name{AmpliconViews-class}
\alias{AmpliconViews-class}
\title{An S4 class for amplicon bisulfite sequencing}
\description{
This class is designed to hold methylat ion statistics and
locations for targeted amplicons from a amplicon bisulfite
sequencing experiment .
\code{\link[AmpliconBiSeq]{AmpliconViews}} function returns
an object of \code{AmpliconViews} class.
}
\section{Slots}{
  \describe{ \item{\code{sampleNames}}{Names of samples in
  a vector} \item{\code{amplicons}}{GRanges object for
  locations of the amplicons} \item{\code{data}}{list
  structure holding the information on amplicons}

  }
}

\section{Details}{
  \code{AmpliconViews} class is desingned to contain
  amplicon information such as methylation, coverage,
  meta-patterns and similarities between CpG profiles.
}

\section{Constructor}{
  see \code{\link{AmpliconViews}}
}

\section{Subsetting}{
  an AmpliconViews object containing multiple amplicons and
  samples can be subsetted using
  \code{\link[AmpliconBiSeq]{getAmplicon}} function.
}

\section{Accessors}{
  The following functions provides access to data slots of
  methylDiff:
  \code{\link[AmpliconBiSeq]{getSamples}},\code{\link[AmpliconBiSeq]{getAmpliconRanges}},
  \code{\link[AmpliconBiSeq]{getAmplicon}},
  \code{\link[AmpliconBiSeq]{getMetaMethylation}}
}
\examples{
library(GenomicRanges)
}

