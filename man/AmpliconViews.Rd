\docType{methods}
\name{AmpliconViews}
\alias{AmpliconViews}
\title{Makes an AmpliconViews object}
\usage{
AmpliconViews(proj, range, tag = NULL, sampleNames, conv = NULL,
  exp.var = 80, noise.tol = "auto", example.size = 100,
  call.matrix = FALSE, verbose = FALSE)
}
\arguments{
  \item{proj}{qProject object from QuasR}

  \item{range}{GRanges object for the amplicon locations}

  \item{tag}{a character vector, containing a additional
  tags that describe the amplicons. The vector length
  should be same as \code{range} object.}

  \item{sampleNames}{sample name as character}

  \item{conv}{minimum conversion efficiency a numeric value
  between 0 and 100, reads below this conversion efficiency
  will be discarded}

  \item{exp.var}{percentage of reads explained by
  meta-methylation patterns. This helps to infer top
  meta-methylation patterns. If set to 90, meta-patterns
  are ranked by their cluster sizes, and patterns that
  attributes to 90 percent of the data are returned.}

  \item{noise.tol}{percentage of variation to be considered
  as noise. This usually corresponds to PCA components that
  explain small percentage of the data to be removed. If
  equals to 'auto', the noise is automatically removed by
  removing components that explain low amount of
  variation.}

  \item{example.size}{size of the example matrix returned.
  example matrix contains example reads from the
  experiment.}

  \item{call.matrix}{\code{logical}, default FALSE, if TRUE
  methylation call matrix that contains methylation calls
  per fragment are returned. This might result in a large
  AmpliconViews object, should be set to FALSE normally.}
}
\value{
returns an \link{AmpliconViews-class} object
}
\description{
Makes an AmpliconViews object from qProject object from
QuasR.
}
\examples{
# a=GRanges(seqnames=,range=IRanges() )
 # av=ampliconView(proj,range=a[2,],samp="amSeq",conv=80,exp.var=80,example.size=100)
}

