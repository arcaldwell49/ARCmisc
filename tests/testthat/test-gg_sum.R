test_that("basic run", {
  mt2 = mtcars
  mt2$cyl = as.factor(mtcars$cyl)
  mt2$am = as.factor(mtcars$am)
  test1 = gg_sum(data = mt2,
                  y = mpg,
                  x = cyl,
                  point_alpha = .7) #+ ggprism::theme_prism()

  test2 = gg_sum(data = mt2,
                  y = mpg,
                  x = cyl,
                  show_points = FALSE,
                  trace = FALSE) #+ ggprism::theme_prism()

  test3 = gg_sum(data = mt2,
                  y = mpg,
                  x = cyl,
                  group = am,
                  show_points = TRUE,
                  trace = TRUE,
                  point_alpha = .5) #+ ggprism::theme_prism()

  test4 =  gg_sum(data = mt2,
                           y = mpg,
                           x = cyl,
                           group = am,
                           show_points = TRUE,
                           trace = FALSE,
                           point_alpha = .5,
                  show_summary = FALSE) #+ ggprism::theme_prism()

  test5 =  gg_sum(data = mt2,
                  y = mpg,
                  x = cyl,
                  group = am,
                  show_points = FALSE,
                  trace = FALSE,
                  point_alpha = .5,
                  show_slab = TRUE)

  test5 =  gg_sum(data = mt2,
                  y = mpg,
                  x = cyl,
                  show_points = FALSE,
                  trace = FALSE,
                  point_alpha = .5,
                  show_slab = TRUE)

  test6 = gg_sum(data = mt2,
                 y = mpg,
                 x = cyl,
                 point_alpha = .7,
                 point_color = "skyblue")


})
