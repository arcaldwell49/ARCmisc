#' Summary plots
#'
#' @description Function for simple factorial experiment plots.
#' Something that should save time for plotting data.
#'
#' @param x,y Columns to be plotted on x or y axis (x should be a factor).
#' @param data A data frame.
#' @param group A grouping factor for plotting different colors/dodge on the same axis.
#' @param panel A factor by which to facet the plot.
#' @param trace TRUE/FALSE to connect summaries with line.
#' @param err_width Width of error bars, or separation of groups.
#' @param sum_stat Type of summary statistic. Options include mean (defualts) or median. Error bars will be the SD or IQR (25th and 75th percentile) for mean and median, respectively.
#' @param show_points,show_summary,show_slab Logical indicator of what to show on plot.
#' @param sum_size,sum_alpha Change parameters for summary statistics.
#' @param point_size,point_alpha Change parameters for point data.
#' @param sum_color Color for summary statistic (only used with groups == NULL).
#' @param point_color Color for data points (only used with groups == NULL).
#'
#' @return A ggplot object.
#'
#' @details
#' Creats simple plots for summarizing factorial designs. See vignettes for examples.
#'
#' @importFrom ggdist geom_dots stat_slab geom_slabinterval
#' @importFrom stats quantile
#' @import ggplot2
#' @export
gg_sum = function(data,
                  x,
                  y,
                  group = NULL,
                  panel = NULL,
                  trace = TRUE,
                  sum_stat = c("mean","median"),
                  show_points = TRUE,
                  show_summary = TRUE,
                  show_slab = FALSE,
                  err_width = .2,
                  sum_size = 1.25,
                  sum_alpha = 1,
                  point_size = 1,
                  point_alpha = 1,
                  sum_color = "darkgrey",
                  point_color = "black"){

  # get data frame -----
  group = deparse(substitute(group))
  panel = deparse(substitute(panel))

  if(group == "NULL"){
    group = NULL
  }
  if(panel == "NULL"){
    panel = NULL
  }

  if(!is.null(group) && !is.null(panel)){
    df = data[c(deparse(substitute(y)),
                deparse(substitute(x)),
                panel,
                group
    )]
    colnames(df) = c("y","x","panel","group")
  }

  if(is.null(group) && !is.null(panel)){
    df = data[c(deparse(substitute(y)),
                deparse(substitute(x)),
                panel
    )]
    colnames(df) = c("y","x","panel")
  }

  if(!is.null(group) && is.null(panel)){
    df = data[c(deparse(substitute(y)),
                deparse(substitute(x)),
                group
    )]
    colnames(df) = c("y","x","group")
  }

  if(is.null(group) && is.null(panel)){
    df = data[c(deparse(substitute(y)),
                deparse(substitute(x))
    )]
    colnames(df) = c("y","x")
  }

  # Settings -----
  sum_stat = match.arg(sum_stat)
  tpanel = ifelse(is.null(panel),
                  FALSE,
                  TRUE)
  tgroup= ifelse(is.null(group),
                 FALSE,
                 TRUE)
  if(!tgroup){
    p1 = plt_ngrp_pt(data = df,
                     tpanel = tpanel,
                     trace = trace,
                     show_points = show_points,
                     sum_size = sum_size,
                     sum_alpha = sum_alpha,
                     sum_stat = sum_stat,
                     point_size = point_size,
                     point_alpha = point_alpha,
                     err_width = err_width,
                     show_summary = show_summary,
                     show_slab = show_slab,
                     sum_color = sum_color,
                     point_color = point_color)
  } else {
    p1 = plt_grp_pt(data = df,
                    tpanel = tpanel,
                    trace = trace,
                    show_points = show_points,
                    sum_size = sum_size,
                    sum_alpha = sum_alpha,
                    point_size = point_size,
                    point_alpha = point_alpha,
                    err_width = err_width,
                    show_summary = show_summary,
                    show_slab = show_slab)
  }


  return(p1)
}

