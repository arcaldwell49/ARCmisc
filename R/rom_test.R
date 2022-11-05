#' Ratio of Means Test
#'
#' @description Performs a test on a ratio of means.
#'
#'
#' @param x,y A numeric or ordered vector, or a character name of one in `data`.
#'   Any missing values (`NA`s) are dropped from the resulting vector. `x` can
#'   also be a formula (see [`stats::wilcox.test()`]), in which case `y` is
#'   ignored.
#' @param data an optional matrix or data frame (or similar: see model.frame) containing the variable. By default the variables are taken from environment.
#' @param null a number indicating the value of the null hypothesis. See [stats::wilcox.test].
#' @param ci Confidence level of the interval.
#' @param ci_method Method for calculating confidence levels.
#' @param vtype Method by which variance is calculated for two independent samples. Ignored in cases of paired samples. Default is "LS".
#' @param bias_c Bias correction for small samples. Default is TRUE.
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
#'     \item{statistic}{The test statistic (z or t).}
#'     \item{estimate}{The ratio of means or "response rato" (RR).}
#'     \item{p.value}{The p-value for the test.}
#'     \item{null.value}{The specified hypothesized value of the RR.}
#'     \item{alternative}{A character string describing the alternative hypothesis.}
#'     \item{method}{The name of the method to report.}
#'     \item{data.name}{The title of the data.}
#'     \item{conf.int}{The confience interval for the RR.}
#'     \item{stderr}{The standard error of the RR.}
#'
#'}
#'
#' @details
#' Ratio of means... To be added
#' vtype details
#'
#' @examples
#' # Use the sleep data
#'
#' rom_otest(extra ~ group, data = sleep, paired = TRUE)
#'
#' @references
#' Lajeunesse, M. J. (2011). On the meta‐analysis of response ratios for studies with correlated and multi‐group designs. Ecology, 92(11), 2049-2055 .https://doi.org/10.1890/11-0423.1
#'
#' Lajeunesse, M. J. (2015). Bias and correction for the log response ratio in ecological meta‐analysis. Ecology, 96(8), 2056-2063. https://doi.org/10.1890/14-2402.1
#'
#' Kieser, M. and Hauschke, D. (1999) Approximate Sample Sizes for Testing Hypotheses about the Ratio and Difference of Two Means. Journal of Biopharmaceutical Studies, 9(4)  641-650.
#'
#' Hauschke, D. et al  (1999). Sample Size Determination for Proving Equivalence Based on the Ratio of Two Means for Normally Distributed Data. Statistics in Medicine, 18,  93-105.
#' @importFrom datawizard ranktransform
#' @importFrom stats complete.cases pnorm qnorm sd
#' @export

rom_test <- function(x,
                     y = NULL,
                     data = NULL,
                     null = 1,
                     ci = 0.95,
                     ci_method = c("t","normal"),
                     vtype = c("LS", "AV", "HO", "AVHO"),
                     alternative = c("two.sided", "less", "greater"),
                     bias_c = TRUE,
                     paired = FALSE,
                     verbose = TRUE,
                     ...) {
  ci_method = match.arg(ci_method)
  alternative <- match.arg(alternative)
  vtype = match.arg(vtype)
  if(paired && vtype != "HO"){
    message("vtype has no effect on paired results.")
  }
  if(!missing(null) && ((length(null) > 1L) || !is.finite(null)) || null <= 0)
    stop("'null' must be a single, positive number")
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
    stop("Only one sample provided. y or data must be provided.")
  } else if(!is.null(y)) {
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  } else {

    DNAME <- gsub("~", "by", Reduce(paste, deparse(x)))
  }
  ## Prep data
  out <- .get_data_2_samples(x=x, y=y, data=data,
                             verbose = TRUE,
                             paired = FALSE,
                             allow_ordered = FALSE,
                             ...)
  x <- out$x
  y <- out$y

  if(is.null(y)){
    stop("Missing y. Only one sample provided.")
  }

  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }

  if (!is.numeric(y)){
    stop("'y' must be numeric")
  }


  if(any(x <=0) || any(y <= 0)){
    stop("Data cannot have values less than or equal to 1.")
  }
  if(length(x) < 1L)
    stop("not enough (finite) 'x' observations")

  if(paired) { ##------------------ paired case -------------------
    m1i = mean(x)
    sd1i = sd(x)
    n1i = length(x)
    m2i = mean(y)
    sd2i = sd(y)
    n2i = length(y)
    ri = cor(x,y)

    if(n1i != n2i){
      warning("Paired samples of varying lengths. Results likely bogus.")
    }

    df1 = min(c(n1i,n2i))

    log_val = logrom_calc(
      paired = TRUE,
      bias_c = bias_c,
      vtype = vtype,
      m1i = m1i,
      sd1i = m2i,
      n1i = n1i,
      m2i = m2i,
      sd2i = sd2i,
      n2i = n2i,
      ri = ri
    )

    if(bias_c){
      METHOD = "Paired Ratio of Means (Bias Corrected)"
    } else {
      METHOD = "Paired Ratio of Means"
    }

  } else { ##------------------------ 2-sample case -------------------------
    if(length(y) < 1L)
      stop("not enough 'y' observations")
    m1i = mean(x)
    sd1i = sd(x)
    n1i = length(x)
    m2i = mean(y)
    sd2i = sd(y)
    n2i = length(y)
   # ri = cor(x,y)
    df1 = n1i+n2i-2

    log_val = logrom_calc(
      paired = FALSE,
      bias_c = bias_c,
      vtype = vtype,
      m1i = m1i,
      sd1i = m2i,
      n1i = n1i,
      m2i = m2i,
      sd2i = sd2i,
      n2i = n2i,
      ri = NULL
    )


    if(bias_c){
      METHOD = "Two-Sample Ratio of Means (Bias Corrected)"
    } else {
      METHOD = "Two-Sample Ratio of Means"
    }

  }

  if (ci_method == "normal") {

    SE <- sqrt(log_val$var_rom)
    interval <- exp(log_val$log_rom + c(-1, 1) * qnorm(1 - alpha / 2) * SE)
    p_value = p_from_z(log_val$log_rom-null, alternative, SE)
    stderr = SE
    names(stderr) = "SE[log(rom)]"
    STATISTIC = log_val$log_rom
    names(STATISTIC) = "log(rom)"
    PARA = SE
    names(PARA) = "SE[log(rom)]"

  }

  if (ci_method == "t"){
    SE <- sqrt(log_val$var_rom)
    interval <- exp(log_val$log_rom + c(-1, 1) * qt(1 - alpha / 2,df1) * SE)
    p_value = p_from_t((log_val$log_rom-null)/SE, alternative, df = df1)

    stderr = SE
    names(stderr) = "SE[log(rom)]"
    STATISTIC = log_val$log_rom
    names(STATISTIC) = "t"
    PARA = df1
    names(PARA) = "df"
  }

  interval = switch(alternative,
                    "two.sided" = interval,
                    "less" = c(0,interval[2]),
                    "greater" = c(interval[1],+Inf))

  attr(interval, "conf.level") <- ci

  rom = exp(log_val$log_rom)
  names(rom) = "RR"
  names(null) <- "Ratio of Means"
  RVAL <- list(statistic = STATISTIC,
               parameter = PARA,
               p.value = p_value,
               null.value = null,
               alternative = alternative,
               method = METHOD,
               estimate = rom,
               conf.int = interval,
               stderr = stderr,
               data.name = DNAME)
  class(RVAL) <- "htest"
  RVAL
}


