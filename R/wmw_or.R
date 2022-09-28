#' Generalized Odds Ratio for Ordinal Data
#'
#' @aliases print.Genodds
#'
#' @description Performs Agresti's Generalized Odds Ratios (GenOR) for two-group, paired samples, and one sample data..
#'
#' @usage
#' genodds(response, group, strata=NULL,
#'         alpha=0.05,ties="split",
#'         nnt=FALSE,verbose=FALSE,upper=TRUE, suppress_interpretation=FALSE,
#'         assume_no_effect=FALSE,
#'         permutation_test=FALSE, nPermutations=5000)
#'

#'
#' @return An object of the class `htest`.
#' \describe{
#'     \item{pooled_lnodds}{The pooled log(odds).}
#'     \item{pooled_lnconf.int}{(1-\code{alpha})\% Confidence intervals for pooled log(odds).}
#'     \item{pooled_SElnodds}{Standard error of pooled log(odds).}
#'     \item{pooled_SElnnull}{Standard error of pooled log(odds) under the null hypothesis.}
#'     \item{pooled_p}{The p-value of the test of pooled log(odds) = 1.}
#'     \item{pooled_rel_statistic}{Statistic of test that strata odds are equal.}
#'     \item{pooled_rel_p}{p-value for test that strata odds are equal.}
#'     \item{relative_lnodds}{A matrix giving the log of the ratio of odds between strata (generalised relative risk ratio).}
#'     \item{relative_selnodds}{A matrix containing the standard error of the log(relative risk ratio).}
#'     \item{results}{A list containing a summary of each strata measure.}
#'     \item{param.record}{A list containing parameters used in the test.}
#'}
#'
#' @details
#' Agresti's generalized odds ratios (GenOR) calculates the odds that,
#' if a pair of observations are randomly selected from
#' two groups, the outcome in one group is higher than the other.
#' This implementation determines the direction of this comparison
#' using factor levels. Odds are given with reference to
#' observations corresponding to the higher \code{group} level
#' resulting in a higher value in \code{response}.
#' The opposite direction can be calculated by either calculating 1/genodds,
#' or by specifying \code{response=1-response} in function input.
#'
#' If \code{nnt=TRUE}, the Number Needed to Treat (NNT) is printed.
#' NNT is a health economics measure and is related to generalised
#' odds ratios through the formula NNT=1+2/(GenOR-1).
#' It measures the expected number of patients required for a
#' treatment to have impacted a patient's outcome.
#' In this implementation, a positive NNT occurs when GenOR>1
#' and corresponds to the number needed to treat in the higher
#' \code{group} level to observe a higher \code{response} value,
#' while a negative NNT occurs when GenOR<1 and corresponds
#' to the number needed to treat in the higher \code{group}
#' level to observe a lower \code{response} value.
#' If the confidence interval for GenOR straddles 1,
#' the confidence interval for NNT is given as the union of disjoint
#' intervals.
#'

#'
#' @examples
#' # Use the alteplase dataset provided by package and calculate genodds
#' df <- alteplase
#' x <- genodds(df$mRS,df$treat,df$time)
#' x
#' print(x,nnt=TRUE)
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
  }

  if (is.null(ci)) {
    alternative <- NULL
    #interval = c(NA,NA)

  } else {
    ci_method <- list(method = "normal")
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
  CORRECTION <- 0
  if(is.null(y)) { ##------------------ 1-sample/paired case -------------------
    z <- x - mu

    n_x = as.double(length(x))
    n_a <- sum(z > 0) + 0.5*sum(z == 0)

    cstat = n_a / n_x
    #if(cstat == 0 || cstat == 1){
    #  stop("Odds ratio cannot be estimated. No overlap with zero/null.")
    #}
    odds = probs_to_odds(cstat)


    if(is.numeric(ci)){
      # Wilson interval for binomial prob
      # Then converted to odds
      zstar <- qnorm(1-alpha/2)
      zstar2 = zstar^2
      p1 <- cstat + 0.5 * zstar2/n_x
      p2 <- zstar * sqrt((cstat * (1 - cstat) + 0.25 * zstar2/n_x)/n_x)
      p3 <- 1 + zstar2/n_x
      lcl <- (p1 - p2)/p3
      ucl <- (p1 + p2)/p3
      interval <- probs_to_odds(c(lcl,ucl))
    }




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
    Rs = Rs + (1 - 0.5) * Rt
    Rd = Rd + 0.5 * Rt
    Pc = sum(p * Rs)
    Pd = sum(p * Rd)

    odds = (Pc/Pd)

    #if(odds == 0 || odds == Inf){
    #  stop("Odds ratio cannot be estimated. No overlap between groups.")
    #}

    ### ----- ci: normal approx ----
    ### bootstrap removed
    if(is.numeric(ci)){

      #odds2 = (Pc/Pd)
      # Not to self: checks calcs from above match matrix calcs
      #if(round(odds2,7) != round(odds,7)){
      #  stop("Matrix broken; likely bug in code.")
      #}
      SEodds = 2/Pd * (sum(p * (odds * Rd - Rs)^2)/N)^0.5
      SElnodds = SEodds/odds
      interval = exp(qnorm(c(alpha/2, 1 - alpha/2), mean = log(odds),
                           sd = SElnodds))

    } #else if(ci_method == "percent"){
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


  }
  out = data.frame(odds = odds)
  if(!is.null(ci_method)){
    out$CI = ci.level
    out$CI_low = if (alternative == "less") 0 else interval[1]
    out$CI_high = if (alternative == "greater") Inf else interval[2]
  }

  class(out) <- c("effectsize_difference", "effectsize_table", "see_effectsize_table", class(out))
  attr(out, "paired") <- paired
  attr(out, "mu") <- mu
  attr(out, "ci") <- ci
  attr(out, "ci_method") <- ci_method
  attr(out, "approximate") <- FALSE
  attr(out, "alternative") <- alternative
  return(out)
}

