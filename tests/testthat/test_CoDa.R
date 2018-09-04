remove(list = ls())
library(MortalityForecast)

# Test model fitting
x = 0:110
y = 1960:2016
D <- MortalityForecast.data$dx[paste(x), paste(y)]
M1 <- fitOeppen(D, x = x, y = y)
M2 <- fitOeppen(D)
vsn <- 1e-200

testCodaFit <- function(M){
  test_that("Test model fitting",{
    expect_s3_class(M, "fitOeppen")
    expect_output(print(M))
    expect_output(print(summary(M)))
    expect_warning(print(M), regexp = NA) # Expect no warning
    expect_warning(summary(M), regexp = NA)
    expect_true(all(fitted(M) >= 0))
    expect_true(all(round(colSums(fitted(M)), vsn) == 1)) # we apply a very small rounding to avoid a false negative
    expect_false(all(is.infinite(fitted(M))))
    expect_false(all(is.na(coef(M))))
    expect_false(all(is.na(resid(M))))
    expect_equal(nrow(M$input$data), length(M$x))
    expect_equal(ncol(M$input$data), length(M$y))
    expect_identical(dim(fitted(M)), dim(resid(M)), dim(M$input$data))
  })
}

for (i in 1:2) testCodaFit(get(paste0("M", i)))


# Test model prediction
P1 <- predict(M1, h = 20)
P2 <- predict(M2, h = 10)

testCodaPred <- function(P){
  test_that("Test model prediction", {
    expect_s3_class(P, "predict.fitOeppen")
    expect_output(print(P))
    expect_true(all(P$predicted.values >= 0))
    expect_true(all(P$conf.intervals$L80 >= 0))
    expect_true(all(P$conf.intervals$L95 >= 0))
    expect_true(all(P$conf.intervals$U80 >= 0))
    expect_true(all(P$conf.intervals$U95 >= 0))
    expect_true(all(round(colSums(P$predicted.values), vsn) == 1))
    expect_true(all(round(colSums(P$conf.intervals$L80), vsn) == 1))
    expect_true(all(round(colSums(P$conf.intervals$L95), vsn) == 1))
    expect_true(all(round(colSums(P$conf.intervals$U80), vsn) == 1))
    expect_true(all(round(colSums(P$conf.intervals$U95), vsn) == 1))
    expect_equal(length(P$y), ncol(P$predicted.values))
    expect_equal(nrow(P$kt), ncol(P$predicted.values))
  })
}

for (i in 1:2) testCodaPred(get(paste0("P", i)))

# ----------------------------------------------
# Test plots
test_that("Test that plots are produced",{
  expect_false(is.null(plot(M1)))
  expect_false(is.null(plot(M1, plotType = "observed")))
  expect_false(is.null(plot(M1, plotType = "fitted")))
  expect_false(is.null(plot(resid(M1))))
  expect_false(is.null(plot(resid(M1), plotType = "scatter")))
  expect_false(is.null(plot(resid(M1), plotType = "colourmap")))
  expect_false(is.null(plot(resid(M1), plotType = "signplot")))
})

# ----------------------------------------------

# Validate input tests

expect_error(fitOeppen(D, x = c(NA,1:109), y = 1960:2016))
expect_error(fitOeppen(D, x = 0:110, y = c(NA,1961:2016)))
expect_error(fitOeppen(D, y = 1960:20160))
expect_error(fitOeppen(D, x = 0:1000))

dNA <- D
dNA[1,1] <- NA
expect_error(fitOeppen(dNA))

dNeg <- dNA
dNeg[1,1] <- -1
expect_error(fitOeppen(dNeg))
