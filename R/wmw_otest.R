#' Rank-based test using WMWodds for one, two, and paired samples
#'
#'
#' @description Performs Wilcoxon-Mann-Whitney odds (WMWodds) hypothesis test,
#' otherwise know as Agresti's Generalized Odds Ratios (GenOR), for two-group, paired samples, and one sample data.
#'
#'
#' @param x,y A numeric or ordered vector, or a character name of one in `data`.
#'   Any missing values (`NA`s) are dropped from the resulting vector. `x` can
#'   also be a formula (see [`stats::wilcox.test()`]), in which case `y` is
#'   ignored.
#' @param mu a number indicating the value around which (a-)symmetry (for
#'   one-sample or paired samples) or shift (for independent samples) is to be
#'   estimated. See [stats::wilcox.test].
#' @param ci Confidence level of the interval.
#' @param ci_method Method by which confidence intervals are calculated. When "normal", the Fisher r-to-z transformation is utilized. For two-sample cases, this can be set to "gamma" to use the method outlined by Agresti using the Goodman-Kruskal gamma approximation.
#' @param alternative a character string specifying the alternative hypothesis;
#'   Controls the type of CI returned: `"two.sided"` (default, two-sided CI),
#'   `"greater"` or `"less"` (one-sided CI). Partial matching is allowed (e.g.,
#'   `"g"`, `"l"`, `"two"`...).
#' @param data An optional data frame containing the variables.
#' @param paired If `TRUE`, the values of `x` and `y` are considered as paired.
#'   This produces an effect size that is equivalent to the one-sample effect
#'   size on `x - y`.
#' @param verbose Toggle warning messages on or off.
#' @param ... Arguments passed to or from other methods. When `x` is a formula,
#'   these can be `subset` and `na.action`.
#'
#' @return An object of the class `htest`.
#' \describe{
#'     \item{statistic}{The calculated odds.}
#'     \item{estimate}{The unnamed calculated odds.}
#'     \item{p.value}{The p-value for the test.}
#'     \item{null.value}{The location parameter mu.}
#'     \item{alternative}{A character string describing the alternative hypothesis.}
#'     \item{method}{The name of the method to report.}
#'     \item{data.name}{The confidence interval for the odds.}
#'     \item{conf.int}{The location parameter mu.}
#'
#'}
#'
#' @details
#' Agresti's generalized odds ratios (GenOR) calculates the odds that,
#' if a pair of observations are randomly selected from
#' two groups, the outcome in one group is higher than the other.
#' This implementation also allows for paired samples comparisons wherein
#' the odds of a random observation being greater than zero/mu is calculated.
#'
#'
#' @examples
#' # Use the sleep data
#'
#' wmw_otest(extra ~ group, data = sleep, paired = TRUE)
#'
#' @references
#' Agresti, A. (1980). Generalized odds ratios for ordinal data.
#' \emph{Biometrics}, 59-67.
#'
#' O'Brien, R. G., & Castelloe, J. (2006, March).
#' Exploiting the link between the Wilcoxon-Mann-Whitney test and a simple odds statistic.
#' In \emph{Thirty-first Annual SAS Users Group International Conference}.
#'
#' Churilov, L., Arnup, S., Johns, H., Leung, T., Roberts,
#' S., Campbell, B. C., Davis, S. M. & Donnan, G. A. (2014).
#' An improved method for simple, assumption-free ordinal analysis of the
#' modified Rankin Scale using generalized odds ratios.
#' \emph{International Journal of Stroke}, 9(8), 999-1005.
#'
#' Howard, G., Waller, J. L., Voeks, J. H., Howard, V. J., Jauch, E. C.,
#' Lees, K. R., ... & Hess, D. C. (2012). A simple, assumption-free,
#' and clinically interpretable approach for analysis of modified Rankin
#' outcomes. \emph{Stroke}, 43(3), 664-669.
#'
#' @export
wmw_otest <- function(x,
                      y = NULL,
                      data = NULL,
                      mu = 0,
                      ci = 0.95,
                      ci_method = "normal",
                      alternative = "two.sided",
                      paired = FALSE,
                      verbose = TRUE,
                      ...) {


  alternative <- match.arg(alternative)
  if(!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")
  if (!is.numeric(ci)) {
    ci <- NULL
    ci_method <- NULL
  }
  if(is.numeric(ci)){
    if(!((length(ci) == 1L)
         && is.finite(ci)
         && (ci > 0)
         && (ci < 1)))
      stop("'ci' must be a single number between 0 and 1")
    ci.level <- if (alternative == "two.sided") ci else 2 * ci - 1
    alpha <- 1 - ci.level
  } else{
    stop("'ci' must be a single number between 0 and 1")
  }

  #if (is.null(ci)) {
  #  alternative <- NULL
  #  #interval = c(NA,NA)

  #} else {
  #  ci_method <- list(method = "normal")
  #}
  if(is.null(data) & is.null(y)){
    DNAME <- deparse(substitute(x))
  } else if(!is.null(y)) {
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  } else {

    DNAME <- gsub("~", "by", Reduce(paste, deparse(x)))
  }
  ## Prep data
  out <- .get_data_2_samples(x, y, data, verbose, ...)
  x <- out$x
  y <- out$y

  if(is.ordered(x)){
    x = as.numeric(x)
  }
  if(!is.numeric(x)) {
    stop("'x' must be numeric or ordered factor")
  }
  if(!is.null(y)) {
    if(is.ordered(y)){
      y = as.numeric(y)
    }
    if(!is.numeric(y)) stop("'y' must be numeric or ordered factor")

    if(paired) {
      if(length(x) != length(y))
        stop("'x' and 'y' must have the same length")
      OK <- complete.cases(x, y)
      x <- x[OK] - y[OK]
      y <- NULL

    }
    else {
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
    }
  } else {
    x <- x[is.finite(x)]
  }

  if(length(x) < 1L)
    stop("not enough (finite) 'x' observations")

  if(is.null(y)) { ##------------------ 1-sample/paired case -------------------
    z <- x - mu

    n_x = as.double(length(x))
    n_a <- sum(z > 0) + 0.5*sum(z == 0)

    cstat = n_a / n_x
    #if(cstat == 0 || cstat == 1){
    #  stop("Odds ratio cannot be estimated. No overlap with zero/null.")
    #}
    odds = probs_to_odds(cstat)
    rho = cstat_to_rb(cstat)
    zstat = rho_to_z(rho)



    if (ci_method == "normal") {
      rho = cstat_to_rb(cstat)
      zstat = rho_to_z(rho)
      # Stolen from effectsize
      maxw <- (n_x ^ 2 + n_x) / 2
      SE <- sqrt((2 * n_x ^ 3 + 3 * n_x ^ 2 + n_x) / 6) / maxw
      interval <- z_to_rho(zstat + c(-1, 1) * qnorm(1 - alpha / 2) * SE) |>
        rb_to_cstat() |> probs_to_odds()
      p_value = p_from_z(zstat, alternative, SE)

      STATISTIC = zstat
      names(STATISTIC) = "z"
      PARA = SE
      names(PARA) = "SD for z"


    } else {
      stop("ci_method must be set to \"normal\" for one/paired sample(s)")
    }
    METHOD = "One-Sample Wilcoxon-Mann-Whitney Odds"



  } else { ##------------------------ 2-sample case -------------------------
    if(length(y) < 1L)
      stop("not enough 'y' observations")
    #r <- rank(c(x - mu, y))
    #n_x <- as.double(length(x))
    #n_y <- as.double(length(y))
    x2 = x - mu

    # Get Mann-Whitney U
    #Ustat <-  sum(r[seq_along(x)]) - n_x * (n_x + 1) / 2
    # Calc c-index
    #cstat = Ustat / (n_x * n_y)
    #if(cstat == 0 || cstat == 1){
    #  stop("Odds ratio cannot be estimated. No overlap between groups so concordance is 0% or 100%!")
    #}
    response = c(y,x2)
    group = c(rep(1,length(y)),rep(2,length(x)))
    crosstab = as.matrix(table(response, group))
    N = sum(crosstab)
    p = crosstab/N
    Rt = p[, 2:1]
    Rs = .rs_mat(p)
    Rd = .rd_mat(p)
    Rs = Rs + (0.5) * Rt
    Rd = Rd + 0.5 * Rt
    Pc = sum(p * Rs)
    Pd = sum(p * Rd)

    odds = (Pc/Pd)
    cstat = odds_to_probs(odds)
    #if(odds == 0 || odds == Inf){
    #  stop("Odds ratio cannot be estimated. No overlap between groups.")
    #}

    ### ----- ci: normal approx ----
    ### bootstrap removed
    if(ci_method == "gamma"){
      n1 <- length(x)
      n2 <- length(y)
      #odds2 = (Pc/Pd)
      # Not to self: checks calcs from above match matrix calcs
      #if(round(odds2,7) != round(odds,7)){
      #  stop("Matrix broken; likely bug in code.")
      #}
      SEodds = 2/Pd * (sum(p * (odds * Rd - Rs)^2)/N)^0.5
      SElnodds = SEodds/odds

      interval = exp(qnorm(c(alpha/2, 1 - alpha/2), mean = log(odds),
                           sd = SElnodds))
      p_value = p_from_odds(odds, alternative, SElnodds)
      SE = 2/Pd*(sum(p*(odds*Rd-Rs)^2)/sum(n1,n2))^0.5

      STATISTIC = log(odds)
      names(STATISTIC) = "log odds"
      PARA = SElnodds
      names(PARA) = "SE of log odds"

    } else if(ci_method == "normal"){
      n1 <- length(x)
      n2 <- length(y)
      SE <- sqrt((n1 + n2 + 1) / (3 * n1 * n2))
      rho = cstat_to_rb(cstat)
      zstat = rho_to_z(rho)
      interval <- z_to_rho(zstat + c(-1, 1) * qnorm(1 - alpha / 2) * SE) |>
        rb_to_cstat() |> probs_to_odds()
      p_value = p_from_z(zstat, alternative, SE)
      STATISTIC = zstat
      names(STATISTIC) = "z"
      PARA = SE
      names(PARA) = "SD for z"


    }

    #else if(ci_method == "percent"){
    #
    #if (insight::check_if_installed("boot", "for estimating CIs", stop = FALSE)) {
    #  data2 = data.frame(
    #    response = c(x,y),
    #    group = c(rep("x", length(x)), rep("y", length(y)))
    #  )
    #  interval <-  .wmw_odds_ci(
    #    data = data2,
    #    ci = ci.level,
    #    alternative = alternative,
    #    iterations = iterations,
    #    mu = mu,
    #    sample = 2
    #  )
    #
    #  ci_method <- list(method = "percentile bootstrap", iterations = iterations)
    #}
    METHOD = "Two-Sample Wilcoxon-Mann-Whitney Odds"

  }
  #out = data.frame(odds = odds)
  #if(!is.null(ci_method)){
  #  out$CI = ci.level
  #  out$CI_low = if (alternative == "less") 0 else interval[1]
  #  out$CI_high = if (alternative == "greater") Inf else interval[2]
  #}
  interval = switch(alternative,
                    "two.sided" = interval,
                    "less" = c(0,interval[2]),
                    "greater" = c(interval[1],+Inf))

  #class(out) <- c("effectsize_difference", "effectsize_table", "see_effectsize_table", class(out))
  #attr(out, "paired") <- paired
  #attr(out, "mu") <- mu
  #attr(out, "ci") <- ci
  #attr(out, "ci_method") <- ci_method
  #attr(out, "approximate") <- FALSE
  #attr(out, "alternative") <- alternative
  #return(out)
  attr(cint, "conf.level") <- ci

  names(odds) = "odds"
  names(mu) <- if(paired || !is.null(y)) "location shift" else "location"
  RVAL <- list(statistic = STATISTIC,
               parameter = PARA,
               p.value = p_value,
               null.value = mu,
               alternative = alternative,
               method = METHOD,
               estimate = odds,
               conf.int = interval,
               data.name = DNAME)
  class(RVAL) <- "htest"
  RVAL
}

