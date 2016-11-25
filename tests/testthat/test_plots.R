library(RunATAC)
library(testthat)
context("Plots produce expected output")

df <- data.frame(x=rnorm(n = 10000, mean = 0, sd = 100),
                 y=rnorm(n = 10000, mean=100, sd = 5))

test_that("Object returned is class ggplot", {
        expect_equal(class(plot_v(df)), c("gg", "ggplot"))
})