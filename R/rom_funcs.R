
logrom_calc = function(paired = FALSE,
                       bias_c = TRUE,
                       vtype = "HO",
                       m1i,
                       sd1i,
                       n1i,
                       m2i,
                       sd2i,
                       n2i,
                       ri = NULL) {
  if (!paired) {
    yi <- log(m1i / m2i)


    ### sample size weighted average of the coefficient of variation in group 1
    mn1wcvi <- .wmean(sd1i / m1i,
                      n1i,
                      na.rm = TRUE)
    ### sample size weighted average of the coefficient of variation in group 2
    mn2wcvi <- .wmean(sd2i / m2i,
                      n2i,
                      na.rm = TRUE)
    ### sample size weighted average of the two CV values

    mnwcvi  <-
      (sum(n1i * (sd1i / m1i)) + sum(n2i * (sd2i / m2i))) / sum((n1i +
                                                                   n2i))

    ### large sample approximation to the sampling variance (does not assume homoscedasticity)
    if (vtype == "LS") {
      vi <-
        sd1i ^ 2 / (n1i * m1i ^ 2) + sd2i ^ 2 / (n2i * m2i ^
                                                               2)
    }
    ### estimator assuming homoscedasticity
    if (vtype == "HO") {
      mi   <- n1i+n2i - 2
      sdpi <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/mi)
      vi <-
        sdpi ^ 2 / (n1i * m1i ^ 2) + sdpi ^ 2 / (n2i * m2i ^
                                                               2)
    }
    ### estimator using the weighted averages of the CV values
    if (vtype == "AV") {
      vi <- mn1wcvi ^ 2 / n1i + mn2wcvi ^ 2 / n2i
    }
    ### estimator using the weighted average of two weighted averages of the CV values
    if (vtype == "AVHO"){
      vi <- mnwcvi ^ 2 * (1 / n1i + 1 / n2i)
    }

  }

  if (paired) {
    yi <- log(m1i / m2i)
    vi <-
      sd1i ^ 2 / (n1i * m1i ^ 2) + sd2i ^ 2 / (n1i * m2i ^ 2) - 2 * ri * sd1i *
      sd2i / (m1i * m2i * n1i)

  }

  if(bias_c){
    J = 0.5 * (sd1i^2 / (n1i * m1i^2) - sd2i^2 / (n2i * m2i^2))
    yi = yi + J

    Jvar = 0.5 * (sd1i^4 / (n1i^2 * m1i^4) - sd2i^4 / (n2i^2 * m2i^4))
    vi = vi + Jvar
  }


  rval = list(
    log_rom = yi,
    var_rom = vi
  )
  return(rval)
}

.wmean = function (x, w, na.rm = FALSE) {
    if (na.rm) {
      i <- !(is.na(x) | is.na(w))
      x <- x
      w <- w
    }
    sum(x * w)/sum(w)
  }

