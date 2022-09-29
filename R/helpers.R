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
#' @keywords internal
.get_data_nested_groups <- function(x, groups = NULL, blocks = NULL, data = NULL,
                                    wide = TRUE, allow_ordered = FALSE,
                                    verbose = TRUE, ...) {
  if (inherits(x, "formula")) {
    if (length(x) != 3L ||
        x[[3L]][[1L]] != as.name("|")) {
      stop("Formula must have the 'x ~ groups | blocks'.", call. = FALSE)
    }

    x[[3L]][[1L]] <- as.name("+")

    x <- .resolve_formula(x, data, ...)

    if (ncol(x) != 3L) {
      stop("Formula must have only two term on the RHS.", call. = FALSE)
    }
  } else if (inherits(x, "data.frame")) {
    x <- as.matrix(x)
  } else if (!inherits(x, c("table", "matrix", "array"))) {
    x <- .resolve_char(x, data)
    groups <- .resolve_char(groups, data)
    blocks <- .resolve_char(blocks, data)

    if (length(x) != length(groups) || length(x) != length(blocks)) {
      stop("x, groups and blocks must be of the same length.", call. = FALSE)
    }

    x <- data.frame(x, groups, blocks)
  }


  if (inherits(x, c("matrix", "array"))) {
    x <- as.table(x)
  }

  if (inherits(x, c("table"))) {
    x <- as.data.frame(x)[, c(3, 2, 1)]
  }

  colnames(x) <- c("x", "groups", "blocks")

  if (allow_ordered && is.ordered(x$x)) {
    x$x <- as.numeric(x$x)
  }
  if (!is.numeric(x$x)) {
    stop("Cannot compute effect size for a non-numeric vector.", call. = FALSE)
  }
  if (!is.factor(x$groups)) x$groups <- factor(x$groups)
  if (!is.factor(x$blocks)) x$blocks <- factor(x$blocks)


  if (verbose && anyNA(x)) {
    warning("Missing values detected. NAs dropped.", call. = FALSE)
  }
  x <- stats::na.omit(x)

  # By this point, the data is in long format
  if (wide) {
    x <- datawizard::data_to_wide(x,
                                  values_from = "x",
                                  id_cols = "blocks",
                                  names_from = "groups"
    )
    x <- x[, -1]
  }
  x
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

odds_to_probs <- function(odds, log = FALSE) {
  if (log) {
    stats::plogis(odds)
  } else {
    stats::plogis(log(odds))
  }
}

probs_to_odds <- function(probs, log = FALSE) {
  if (log) {
    stats::qlogis(probs)
  } else {
    exp(stats::qlogis(probs))
  }
}

rb_to_odds <- function(rb) {
  probs_to_odds(rb_to_cstat(rb))
}

rb_to_cstat <- function(rb) {
  (rb + 1) / 2
}

cstat_to_rb <- function(cstat){
  2*cstat-1
}

z_to_rho <- function(z){
  tanh(z)
}

rho_to_z <- function(rho){
  atanh(rho)
}

p_from_z <- function(z, alternative = "two.sided", se = 1){

  p = switch(alternative,
         "two.sided" = 2*pnorm(-abs(z), sd = se),
         "greater" = pnorm(z, sd = se, lower.tail = FALSE),
         "less" = pnorm(z, sd = se, lower.tail = TRUE))

  return(p)
}

p_from_odds = function(odds, alternative = "two.sided", se = 1){
  p = switch(alternative,
             "two.sided" = pnorm(abs(log(odds)),sd=se,lower.tail=FALSE)*2,
             "greater" = pnorm(log(odds),sd=se,lower.tail=FALSE),
             "less" = pnorm(log(odds),sd=se,lower.tail=TRUE))
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
  sr = datawizard::ranktransform(z,verbose = FALSE,
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

