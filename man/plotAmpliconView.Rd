\docType{methods}
\name{plotAmpliconView}
\alias{plotAmpliconView}
\title{plot ampliconView object with Gviz}
\usage{
plotAmpliconView(obj, example.reads = FALSE, h.panel = NULL,
  sim.heat = TRUE, h.heat = 0.5, sim.col = colorRampPalette(c("blue",
  "yellow", "red"))(50), sim.na.threshold = 10, sim.measure = "similarity",
  newpage = TRUE, meta.exp = 80, ...)
}
\arguments{
  \item{obj}{\code{\link{AmpliconViews}} object with one
  amplicon.  Use \code{\link{getAmplicon}} function to get
  one amplicon from one sample if necessary.}

  \item{example.reads}{if TRUE (default: FALSE), example
  reads sampled from actual reads will be displayed. Each
  \code{AmpliconViews} object will have a set of example
  reads that correspond to meta-methylation profiles.}

  \item{h.panel}{height of the panels, this should be a
  numeric vector with length equaling to the number of
  panels plotted}

  \item{sim.heat}{logical, if TRUE a similarity heatmap is
  drawn under the tracks, showing similarity of base
  methylation profiles}

  \item{h.heat}{relative height of the similarity heatmap,
  should be a numeric value between 0 and 0.5. Smaller the
  value, smaller the heatmap on the plot.}

  \item{sim.col}{colors for similarity heatmap default:
  colorRampPalette(c("black", "yellow","purple"))(50)}

  \item{sim.na.threshold}{number of pair-wise observation
  needed to calculate similarity scores from CpGs with low
  number of overlap may be unreliable}

  \item{sim.measure}{similarity measure, should be one of
  "similarity","tanay","msimilarity",}

  \item{meta.exp}{% of data set explained by meta-profiles,
  the meta-profiles explaining data set will be plotted}

  \item{newpage}{if TRUE the plot will start a fresh
  graphical device.}

  \item{...}{Other arguments to Gviz::plotTracks}
}
\description{
plot ampliconView object with Gviz
}
\examples{
data(ampViewEx) # load example data
      myAmp=getAmplicon(ampViewEx,"mock4","chr18_69674375_69674775")
      plotAmpliconView(myAmp)
}

