# --------------------------------------------------- #
# Author: Marius D. Pascariu
# Last update: Wed Nov 14 15:40:37 2018
# --------------------------------------------------- #
remove(list = ls())


x  <- 0:110
y  <- 1985:2014
dx <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
M  <- fit_MEM(data = dx, n = 6)
P1 <- predict(M, h = 16, x.h = 0:110, jumpchoice = 'fit')
P2 <- predict(M, h = 16, x.h = 0:130, jumpchoice = 'actual')


test_that("Test for positive densities.", {
  expect_true(all(is.na(M$fitted.values[,1])))
  expect_true(all(M$fitted.values[,-1] >= 0))
  expect_true(all(P1$predicted.values >= 0))
  expect_true(all(P2$predicted.values >= 0))
})


test_that("Test that fitMaxEntMortality plots are produced.", {
  res <- resid(M)
  expect_error(plot(res, plotType = 'scatterxxx'))
  expect_false(is.null(plot(res)))
  expect_false(is.null(plot(res, plotType = 'scatter')))
  expect_false(is.null(plot(res, plotType = 'colourmap')))
  expect_false(is.null(plot(res, plotType = 'signplot')))
  expect_false(is.null(plot(M, plotType = 'observed')))
  expect_false(is.null(plot(M, plotType = 'fitted')))
  expect_false(is.null(plot(P1)))
  expect_false(is.null(plot(P1, plotType = "upper")))
  expect_false(is.null(plot(P1, plotType = "lower")))
})

