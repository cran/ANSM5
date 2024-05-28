test_that("Example 11.3", {
  expect_equal(spearman.beta(ch11$reportedtime, ch11$parentlimit, H0 = 1)$stat,
               0.56666667)
  expect_equal(spearman.beta(ch11$reportedtime, ch11$parentlimit,
                             H0 = 1)$pval.exact, 0.59484127)
  expect_equal(spearman.beta(ch11$reportedtime, ch11$parentlimit, H0 = 1,
                             do.CI = TRUE)$CI.sample.lower, 0.5)
  expect_equal(spearman.beta(ch11$reportedtime, ch11$parentlimit, H0 = 1,
                             do.CI = TRUE)$CI.sample.upper, 1.925)
})
