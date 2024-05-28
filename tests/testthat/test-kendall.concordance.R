test_that("Exercise 10.11", {
  expect_equal(kendall.concordance(ch10$marks, ch10$script, ch10$examiner,
                                   do.exact = FALSE,
                                   do.asymp = TRUE)$pval.asymp, 0.00049976417)
  expect_equal(kendall.concordance(ch10$marks, ch10$examiner, ch10$script,
                                   do.exact = FALSE,
                                   do.asymp = TRUE)$pval.asymp, 0.00247788356)
})
