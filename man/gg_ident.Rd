% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gg_ident.R
\name{gg_ident}
\alias{gg_ident}
\title{Line of Identity Plots}
\usage{
gg_ident(
  x,
  y,
  data,
  group = NULL,
  panel = NULL,
  point_geom = c("point", "jitter"),
  line_color = "black",
  line_type = "solid",
  line_alpha = 1,
  line_size = 1,
  point_color = "black",
  point_alpha = 1,
  point_size = 1
)
}
\arguments{
\item{x, y}{Columns to be plotted on x or y axis (Both should be numeric).}

\item{data}{A data frame.}

\item{group}{A grouping factor for plotting different colors/dodge on the same axis.}

\item{panel}{A factor by which to facet the plot.}

\item{point_geom}{The type of geom to represent the data. Choices include "point" (default) or "jitter".}

\item{line_size, line_color, line_alpha, line_type}{Change line of identity parameters.}

\item{point_size, point_alpha, point_color}{Change parameters for point data.}
}
\value{
A ggplot object
}
\description{
Function for line-of-identity visualizations.
Something that should save time for plotting data for paired samples.
}
\details{
Creates simple "line of identity" plots for summarizing paired differences.
The ggsides pacakge is a good complementary addition to this function.
}
