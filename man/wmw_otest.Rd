% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wmw_otest.R
\name{wmw_otest}
\alias{wmw_otest}
\title{Rank-based test using WMWodds for one, two, and paired samples}
\usage{
wmw_otest(
  x,
  y = NULL,
  data = NULL,
  mu = 0,
  ci = 0.95,
  ci_method = c("normal", "gamma"),
  alternative = c("two.sided", "less", "greater"),
  paired = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x, y}{A numeric or ordered vector, or a character name of one in `data`.
Any missing values (`NA`s) are dropped from the resulting vector. `x` can
also be a formula (see [`stats::wilcox.test()`]), in which case `y` is
ignored.}

\item{data}{An optional data frame containing the variables.}

\item{mu}{a number indicating the value around which (a-)symmetry (for
one-sample or paired samples) or shift (for independent samples) is to be
estimated. See [stats::wilcox.test].}

\item{ci}{Confidence level of the interval.}

\item{ci_method}{Method by which confidence intervals are calculated. When "normal", the Fisher r-to-z transformation is utilized. For two-sample cases, this can be set to "gamma" to use the method outlined by Agresti using the Goodman-Kruskal gamma approximation.}

\item{alternative}{a character string specifying the alternative hypothesis;
Controls the type of CI returned: `"two.sided"` (default, two-sided CI),
`"greater"` or `"less"` (one-sided CI). Partial matching is allowed (e.g.,
`"g"`, `"l"`, `"two"`...).}

\item{paired}{If `TRUE`, the values of `x` and `y` are considered as paired.
This produces an effect size that is equivalent to the one-sample effect
size on `x - y`.}

\item{verbose}{Toggle warning messages on or off.}

\item{...}{Arguments passed to or from other methods. When `x` is a formula,
these can be `subset` and `na.action`.}
}
\value{
An object of the class `htest`.
\describe{
    \item{statistic}{The calculated odds.}
    \item{estimate}{The unnamed calculated odds.}
    \item{p.value}{The p-value for the test.}
    \item{null.value}{The location parameter mu.}
    \item{alternative}{A character string describing the alternative hypothesis.}
    \item{method}{The name of the method to report.}
    \item{data.name}{The data name.}
    \item{conf.int}{The confidence interval fro the odds.}

}
}
\description{
Performs Wilcoxon-Mann-Whitney odds (WMWodds) hypothesis test,
otherwise know as Agresti's Generalized Odds Ratios (GenOR), for two-group, paired samples, and one sample data.
}
\details{
Agresti's generalized odds ratios (GenOR) calculates the odds that,
if a pair of observations are randomly selected from
two groups, the outcome in one group is higher than the other.
This implementation also allows for paired samples comparisons wherein
the odds of a random observation being greater than zero/mu is calculated.
}
\examples{
# Use the sleep data

wmw_otest(extra ~ group, data = sleep, paired = TRUE)

}
\references{
Agresti, A. (1980). Generalized odds ratios for ordinal data.
\emph{Biometrics}, 59-67.

O'Brien, R. G., & Castelloe, J. (2006, March).
Exploiting the link between the Wilcoxon-Mann-Whitney test and a simple odds statistic.
In \emph{Thirty-first Annual SAS Users Group International Conference}.

Churilov, L., Arnup, S., Johns, H., Leung, T., Roberts,
S., Campbell, B. C., Davis, S. M. & Donnan, G. A. (2014).
An improved method for simple, assumption-free ordinal analysis of the
modified Rankin Scale using generalized odds ratios.
\emph{International Journal of Stroke}, 9(8), 999-1005.

Howard, G., Waller, J. L., Voeks, J. H., Howard, V. J., Jauch, E. C.,
Lees, K. R., ... & Hess, D. C. (2012). A simple, assumption-free,
and clinically interpretable approach for analysis of modified Rankin
outcomes. \emph{Stroke}, 43(3), 664-669.
}
