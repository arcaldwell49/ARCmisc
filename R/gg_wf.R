#' Waterfall Plots
#'
#' @description Function for waterfall differences visualizations.
#' Something to make effects from paired (pre-post) intervention studies.
#'
#' @param y Column to be plotted on  y axis (should be numeric).
#' @param data A data frame.
#' @param group A grouping factor for plotting different colors/dodge on the same axis.
#' @param panel A factor by which to facet the plot.
#' @param point_geom The type of geom to represent the data. Choices include "bar" (default), "point", or "line".
#' @param point_size,point_alpha,point_color Change parameters for point data. Size does not work for when geom is "bar".
#' @param line_size,line_color,line_alpha,line_type Change horizontal line parameters.
#' @param hline The position of the horizontal line on the y-axis (default = 0).
#' @param line_size,line_color,line_alpha,line_type Change line parameters.
#'
#' @return A ggplot object
#'
#' @details
#' Creates a simple waterfall plot.
#'
#' @references
#' Gillespie T. W. (2012). Understanding waterfall plots. Journal of the advanced practitioner in oncology, 3(2), 106–111. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4093310/
#'
#' Shao et al. (2014) Use and Misuse of Waterfall Plots, JNCI: Journal of the National Cancer Institute, 106(12), 331. https://doi.org/10.1093/jnci/dju331
#' @import ggplot2
#' @importFrom dplyr arrange
#' @export


gg_wf = function(y,
                 data,
                 group = NULL,
                 panel = NULL,
                 point_geom = c("bar","point","line"),
                 point_color = "black",
                 point_alpha = 1,
                 point_size = 1,
                 hline = 0,
                 line_color = "black",
                 line_type = "solid",
                 line_alpha = 1,
                 line_size = 1){
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
                panel,
                group
    )]
    colnames(df) = c("y","panel","group")
  }

  if(is.null(group) && !is.null(panel)){
    df = data[c(deparse(substitute(y)),
                panel
    )]
    colnames(df) = c("y","panel")
  }

  if(!is.null(group) && is.null(panel)){
    df = data[c(deparse(substitute(y)),
                group
    )]
    colnames(df) = c("y","group")
  }

  if(is.null(group) && is.null(panel)){
    df = data[c(deparse(substitute(y)))]
    colnames(df) = c("y")
  }

  # Settings -----
  tpanel = ifelse(is.null(panel),
                  FALSE,
                  TRUE)
  tgroup = ifelse(is.null(group),
                  FALSE,
                  TRUE)

  df2 <- arrange(df, y)

  df2$x = 1:nrow(df2)

  if(!tgroup){
    p1 = ggplot(df2,
                aes(x=x,y=y))
  } else {
    p1 = ggplot(df2,
                aes(x=x,y=y,
                    group=group,
                    color = group,
                    fill = group))
  }

  if(point_geom == "point"){
    p1 = p1 +
      if(tgroup){
        geom_point(size = point_size,
                   alpha = point_alpha)
      } else{
        geom_point(size = point_size,
                   alpha = point_alpha,
                   color = point_color)
      }

  }

  if(point_geom == "line"){
    p1 = p1 +
      if(tgroup){
        geom_line(size = point_size,
                  alpha = point_alpha)
      } else{
        geom_line(size = point_size,
                  alpha = point_alpha,
                  color = point_color)
      }
  }

  if(point_geom == "bar"){
    p1 = p1 +
      if(tgroup){
        geom_bar(stat = "identity",
                 alpha = point_alpha)
      } else{
        geom_bar(stat = "identity",
                 alpha = point_alpha,
                 color = point_color)
      }

  }

  if(tpanel){
    p1 = p1 +
      facet_wrap(~panel)
  }

  p1 = p1 +
    geom_hline(yintercept = hline,
               color = line_color,
               linetype = line_type,
               alpha = line_alpha,
               size = line_size) +
    labs(x = "") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  return(p1)
}
