library(testthat)
source("./R/functions.R")

test_results <- auto_test("./R", "./tests", reporter = "summary")
