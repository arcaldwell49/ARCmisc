# Functions to help other functions
#' @keywords internal
# Data ------
.get_data_2_samples <- function(x, y = NULL, data = NULL,
                                paired = FALSE, allow_ordered = FALSE,
                                verbose = TRUE, ...) {
  if (inherits(x, "formula")) {
    # Validate:
    if (length(x) != 3L) {
      stop("Formula must have one of the following forms:",
           "\n\ty ~ group,\n\ty ~ 1,\n\tPair(x,y) ~ 1",
           call. = FALSE
      )
    }

    # Pull columns
    mf <- .resolve_formula(x, data, ...)

    if (ncol(mf) > 2L) {
      stop("Formula must have only one term on the RHS.", call. = FALSE)
    }

    x <- mf[[1]]
    y <- NULL
    if (ncol(mf) == 2L) {
      y <- mf[[2]]
      if (!is.factor(y)) y <- factor(y)
    }
  } else {
    # Test if they are they are column names
    x <- .resolve_char(x, data)
    y <- .resolve_char(y, data)
  }


  # If x is ordered and allowed to be...
  if (allow_ordered && is.ordered(x)) {
    if (is.ordered(y)) {
      if (!isTRUE(all.equal(levels(y),levels(x)))) {
        stop("x and y are ordered, but do not have the same levels.", call. = FALSE)
      }
      y <- as.numeric(y)
    }

    x <- as.numeric(x)
  }

  # x should be a numeric vector or a Pair:
  if (!is.numeric(x)) {
    stop("Cannot compute effect size for a non-numeric vector.", call. = FALSE)
  } else if (inherits(x, "Pair")) {
    x <- x[, 1] - x[, 2]
    y <- NULL
  }


  # y should be NULL, numeric, or a factor:
  if (!is.null(y)) {
    if (!is.numeric(y)) {
      if (insight::n_unique(y) != 2) {
        stop("Grouping variable y must have exactly 2 levels.", call. = FALSE)
      }

      if (length(x) != length(y)) {
        stop("Grouping variable must be the same length.", call. = FALSE)
      }

      data <- Filter(length, split(x, y))
      x <- data[[1]]
      y <- data[[2]]
    }

    if (verbose && insight::n_unique(y) == 2) {
      warning("'y' is numeric but has only 2 unique values.",
              "\nIf this is a grouping variable, convert it to a factor.",
              call. = FALSE
      )
    }
  }

  if (verbose && (anyNA(x) || anyNA(y))) {
    warning("Missing values detected. NAs dropped.", call. = FALSE)
  }

  if (paired && !is.null(y)) {
    o <- stats::complete.cases(x, y)
    x <- x[o]
    y <- y[o]
  } else {
    x <- stats::na.omit(x)
    y <- stats::na.omit(y)
  }


  list(x = x, y = y)
}


# Formula ----

#' @keywords internal
#' @importFrom stats model.frame na.pass
.resolve_formula <- function(formula, data, subset, na.action, ...) {
  cl <- match.call(expand.dots = FALSE)
  cl[[1]] <- quote(stats::model.frame)
  if ("subset" %in% names(cl)) {
    cl$subset <- substitute(subset)
  }
  cl$... <- NULL
  cl$na.action <- stats::na.pass
  eval.parent(cl)
}

#' @keywords internal
.resolve_char <- function(nm, data) {
  if (is.character(nm) && length(nm) == 1L) {
    if (is.null(data)) {
      stop("Please provide data argument.", call. = FALSE)
    } else if (!nm %in% names(data)) {
      stop("Column ", nm, " missing from data.", call. = FALSE)
    }

    return(data[[nm]])
  }
  nm
}

# conversions -----

odds_to_pr <- function(x, log = FALSE) {
  if (log) {
    stats::plogis(x)
  } else {
    stats::plogis(log(x))
  }
}

pr_to_odds <- function(x, log = FALSE) {
  if (log) {
    stats::qlogis(x)
  } else {
    exp(stats::qlogis(x))
  }
}

rb_to_odds <- function(x) {
  pr_to_odds(rb_to_cstat(x))
}

rb_to_cstat <- function(x) {
  (x + 1) / 2
}

cstat_to_rb <- function(x){
  2*x-1
}

z_to_rho <- function(x){
  tanh(x)
}

rho_to_z <- function(x){
  atanh(x)
}

p_from_z <- function(x, alternative = "two.sided", se = 1){

  p = switch(alternative,
         "two.sided" = 2*pnorm(-abs(x), sd = se),
         "greater" = pnorm(x, sd = se, lower.tail = FALSE),
         "less" = pnorm(x, sd = se, lower.tail = TRUE))

  return(p)
}

p_from_t <- function(x, alternative = "two.sided", df){

  p = switch(alternative,
             "two.sided" = 2*pt(abs(x), df = df, lower.tail = FALSE),
             "greater" = pt(x, df = df, lower.tail = FALSE),
             "less" = pt(x, df = df, lower.tail = TRUE))

  return(p)
}

p_from_odds = function(x, alternative = "two.sided", se = 1){
  p = switch(alternative,
             "two.sided" = pnorm(abs(log(x)),sd=se,lower.tail=FALSE)*2,
             "greater" = pnorm(log(x),sd=se,lower.tail=FALSE),
             "less" = pnorm(log(x),sd=se,lower.tail=TRUE))
}


# Matrix functions ------

.rd_mat = function(p) {

  nc = nrow(p) # Number of categories
  ng = ncol(p) # Number of groups
  minc = nc - 1

  Rd = matrix(0, nrow = nc, ncol = ng)

  for (c in 1:nc) {
    pc = c+1
    if(pc <= nc){
      for (i in pc:nc) {
        Rd[c, 2] = Rd[c, 2] + p[i, 1]
      }
    }

    if(c-1 >=1){
      for (j in 1:(c-1)) {
        Rd[c, 1] = Rd[c, 1] + p[j, 2]
      }
    }

  }

  return(Rd)
}

.rs_mat = function(p) {

  nc = nrow(p) # Number of categories
  ng = ncol(p) # Number of groups
  minc = nc - 1

  Rs = matrix(0, nrow = nc, ncol = ng)

  for (c in 1:nc) {
    pc = c+1
    if(pc <= nc){
      for (i in pc:nc) {

        Rs[c, 1] = Rs[c, 1] + p[i, 2]
      }
    }

    if(c-1 >=1){
      for (j in 1:(c-1)) {
        Rs[c, 2] = Rs[c, 2] + p[j, 1]
      }
    }

  }

  return(Rs)
}

zsimp = function(z){
  sr = ranktransform(z,verbose = FALSE,
                                 sign = TRUE)
  sr = ifelse(is.na(sr),0,sr)
  z <- sum(sr)
  return(z)
}
zse_simp = function(z){
  sr = datawizard::ranktransform(z,verbose = FALSE,
                                 sign = TRUE)
  sr = ifelse(is.na(sr),0,sr)
  return(sqrt(sum(sr ^ 2)))
}

# plot summaries ----

## need to fix summary functions

### This... does not work...
sum_sdl = function(x){
  mu = mean(x, na.rm = TRUE)
  lower = mu - sd(x, na.rm = TRUE)
  upper = mu + sd(x, na.rm = TRUE)
  vec = data.frame(y =  mu,
                   ymin = lower,
                   ymax = upper)

  return(vec)
}

# Median + IQR summary
#stat_summary(
#  mapping = aes(x = cut, y = depth),
#  fun.min = function(z) { quantile(z,0.25) },
#  fun.max = function(z) { quantile(z,0.75) },
#  fun = median)

plt_ngrp_pt = function(data,
                       tpanel = FALSE,
                       trace = FALSE,
                       sum_stat = "mean",
                       err_stat = "mean_sdl",
                       show_points = TRUE,
                       sum_size = 1.25,
                       sum_alpha = 1,
                       point_size = 1,
                       point_alpha = 1,
                       err_width = .2,
                       show_summary = TRUE,
                       show_slab = FALSE,
                       sum_color = "black",
                       point_color = "darkgrey"){
  if(show_points){
    x_nudge = err_width
  } else {
    x_nudge = 0
  }
  g1 = ggplot(data,
              aes(x=x,y=y,
                  group = 1))

  if(show_points){
    g1 = g1 + geom_dots(dotsize = 1,
                        binwidth = .5*point_size,
                        side = "topleft",
                        fill = point_color,
                        scale = .25,
                        alpha = point_alpha,
                        color = point_color,
                        layout = "weave")
    #geom_point(position=position_jitter(height = 0,
    #                          width = .2),
    #           alpha = point_alpha,
    #           size = point_size)
  }

  if(show_slab){
    g1 = g1 + stat_slab(scale = err_width*2,
                        #binwidth = .5*point_size,
                        side = "topleft",
                        #fill = "black",
                        #scale = err_width,
                        alpha = point_alpha,
                        color = "black",
                        fill = sum_color,
                        #layout = "weave",
                        position=position_dodge(err_width*2)
    )
  }

  if(trace){
    g1 = g1 + stat_summary(
      fun = sum_stat,
      geom = "line",
      #shape = "square",
      #size = 1.25 * 2,
      alpha = .75*sum_alpha,
      size = sum_size*.75,
      color = sum_color,
      position = position_nudge(x = x_nudge, y = 0)
    )
  }

  if(show_summary){

    if(sum_stat == "mean"){
      g1 = g1 +
        stat_summary(fun.data = err_stat,
                     fun.args = list(mult = 1),
                     alpha = sum_alpha,
                     width = err_width,
                     geom = "errorbar",
                     colour = sum_color, size = sum_size,
                     position=position_nudge(x = x_nudge, y = 0)) +
        stat_summary(fun = sum_stat, geom = "point",
                     shape = "square",
                     size = sum_size*2,
                     alpha = sum_alpha,
                     color = sum_color,
                     position=position_nudge(x = x_nudge, y = 0))
    } else {

      g1 = g1 +
        stat_summary(
          fun.min = function(z) {
            quantile(z, 0.25)
          },
          fun.max = function(z) {
            quantile(z, 0.75)
          },
          geom = "errorbar",
          alpha = sum_alpha,
          width = err_width,
          colour = sum_color,
          size = sum_size,
          position = position_nudge(x = x_nudge,
                                    y = 0)
        ) +
        stat_summary(
          fun = "median",
          geom = "point",
          shape = "square",
          size = sum_size * 2,
          alpha = sum_alpha,
          color = sum_color,
          position = position_nudge(x = x_nudge, y = 0)
        )
    }

  }

  if(tpanel){
    g1 = g1 + facet_wrap(~panel)
  }

  return(g1)
}

plt_grp_pt = function(data,
                      tpanel = FALSE,
                      trace = FALSE,
                      sum_stat = "mean",
                      err_stat = "mean_sdl",
                      show_points = TRUE,
                      sum_size = 1.25,
                      sum_alpha = 1,
                      point_size = 1,
                      point_alpha = 1,
                      err_width = .2,
                      show_summary = TRUE,
                      show_slab = show_slab){
  if(show_points){
    x_nudge = err_width
  } else {
    x_nudge = 0
  }
  g1 = ggplot(data,
              aes(x=x,y=y,
                  color = group,
                  fill = group))

  if(show_points){
    g1 = g1 + geom_dots(dotsize = 1,
                        binwidth = .5*point_size,
                        side = "topleft",
                        #fill = "black",
                        scale = err_width,
                        alpha = point_alpha,
                        color = "black",
                        layout = "weave",
                        position=position_dodge(err_width*2)
    )
    #geom_point(position=position_jitter(height = 0,
    #                          width = .2),
    #           alpha = point_alpha,
    #           size = point_size)
  }

  if(show_slab){
    g1 = g1 + stat_slab(scale = err_width*2,
                        #binwidth = .5*point_size,
                        side = "topleft",
                        #fill = "black",
                        #scale = err_width,
                        alpha = point_alpha,
                        color = "black",
                        #layout = "weave",
                        position=position_dodge(err_width*2)
    )
  }

  if(trace){
    g1 = g1 + stat_summary(
      aes(group = group),
      fun = sum_stat,
      geom = "line",
      #shape = "square",
      #size = 1.25 * 2,
      alpha = .75*sum_alpha,
      size = sum_size*.75,
      #color = "darkgray",
      position=position_dodge(err_width*2)
    )
  }

  if(show_summary){
  g1 = g1 +  stat_summary(fun.data = err_stat,
                 fun.args = list(mult = 1),
                 alpha = sum_alpha,
                 #width = err_width,
                 geom = "pointrange",
                 shape = "square",
                 #colour = "darkgray",
                 size = .75*sum_size,
                 position=position_dodge(err_width*2))
  }

  if(tpanel){
    g1 = g1 + facet_wrap(~panel)
  }

  return(g1)
}

# Global variables -----

utils::globalVariables(c("x","y", "group"))
