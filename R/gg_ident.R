#' Line of Identity Plots
#'
#' @description Function for line-of-identity visualizations.
#' Something that should save time for plotting data for paired samples.
#'
#' @param x,y Columns to be plotted on x or y axis (x should be a factor).
#' @param data A data frame.
#' @param group A grouping factor for plotting different colors/dodge on the same axis.
#' @param panel A factor by which to facet the plot.
#' @param point_geom The type of geom to represent the data. Choices include "point" (default) or "jitter".
#' @param point_size,point_alpha,point_color Change parameters for point data.
#' @param line_size,line_color,line_alpha,line_type Change line of identity parameters.
#'
#' @return A ggplot object
#'
#' @details
#' Creates simple "line of identity" plots for summarizing paired differences.
#'
#' @import ggplot2
#' @export
gg_ident = function(x, y,
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
                    point_size = 1){
  point_geom = match.arg(point_geom)

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
  tpanel = ifelse(is.null(panel),
                  FALSE,
                  TRUE)
  tgroup= ifelse(is.null(group),
                 FALSE,
                 TRUE)

  scalemin = min(c(min(df$x, na.rm = TRUE),min(df$y, na.rm = TRUE)))
  scalemax = max(c(max(df$x, na.rm = TRUE),max(df$y, na.rm = TRUE)))

  x_lab = deparse(substitute(x))
  y_lab = deparse(substitute(y))

  p1 = ggplot(df,
              aes(x=x, y=y)) +
    geom_abline(intercept = 0, slope = 1,
                linetype = line_type,
                color = line_color,
                alpha = line_alpha,
                size = line_size) +
    scale_y_continuous(limits = c(scalemin,scalemax))+
    scale_x_continuous(limits = c(scalemin,scalemax)) +
    labs(x = x_lab,
         y = y_lab)

  if(!tgroup){



    if(point_geom == "point"){
      p1 = p1 + geom_point(color = point_color,
                           alpha = point_alpha,
                           size = point_size)
    }

    if(point_geom == "jitter"){
      p1 = p1 + geom_jitter(color = point_color,
                            alpha = point_alpha,
                            size = point_size)
    }
  }

  if(tgroup){
    if(point_geom == "point"){
      p1 = p1 + geom_point(aes(color = group),
                           alpha = point_alpha,
                           size = point_size)
    }

    if(point_geom == "jitter"){
      p1 = p1 + geom_jitter(aes(color = group),
                            alpha = point_alpha,
                            size = point_size)
    }
  }

  if(tpanel){
    p1 = p1 +
      facet_wrap(~panel)

  }

  return(p1)
}

