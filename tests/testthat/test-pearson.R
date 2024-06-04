test_that("Section 10.1.2", {
  expect_equal(pearson(ch10$q1, ch10$q2, alternative = "greater",
                       do.asymp = TRUE, do.exact = FALSE)$stat, 0.37293673)
  expect_equal(pearson(ch10$q1, ch10$q2, alternative = "greater",
                       do.asymp = TRUE, do.exact = FALSE)$pval.asymp,
               0.116244594)
})

test_that("Example 11.2", {
  H0 <- 1
  tmp <- pearson(ch11$parentlimit, ch11$reportedtime - H0 * ch11$parentlimit,
                 alternative = "two.sided")
  expect_equal(tmp$stat, 0.139175550)
  expect_equal(tmp$pval.exact, 0.79246032)
})
