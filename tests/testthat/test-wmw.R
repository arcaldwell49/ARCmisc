test_that("basic run", {
  data(sleep)
  data(mtcars)
  expect_error(wmw_otest())
  mt2 = mtcars
  mt2$am = as.factor(mt2$am)
  test1 = wmw_otest(extra ~ group,
                    data = sleep)
  test2 = wmw_otest(extra ~ group,
                    paired = TRUE,
                    data = sleep)
  test3 = wmw_otest(extra ~ group,
                    paired = TRUE,
                    mu = 1,
                    data = sleep)
  test4 = wmw_otest(sleep$extra)

  test5 = wmw_otest(mpg ~ 1,
                    mu = 20,
                    data = mtcars)

  test6 = wmw_otest(mpg ~ am,
                    paired = FALSE,
                    data = mt2)
  expect_error(wmw_otest(mpg ~ am,
                    paired = TRUE,
                    data = mt2))
  })
