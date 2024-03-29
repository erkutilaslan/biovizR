% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maplot_dge.R
\name{maplot_dge}
\alias{maplot_dge}
\title{MA-plot visualization of DGE analysis.}
\usage{
maplot_dge(
  dge_data,
  FDR = 0.01,
  FC = 0.5,
  TOP = 10,
  type = "deseq2",
  header = ""
)
}
\arguments{
\item{dge_data}{Path to the input file.}

\item{FDR}{Default 0.01. Adjust FDR value.}

\item{FC}{Default 0.5. Adjust log2FC threshold.}

\item{TOP}{Default 10. Adjust top number of DE genes to visualize.}

\item{type}{Default deseq2. Specify input datatype.}

\item{header}{Default is empty. Set a title for the MA-plot.}
}
\value{
An ma-plot of RNA-Seq DGE data.
}
\description{
This function visualizes DGE results as an MA-plot.
}
