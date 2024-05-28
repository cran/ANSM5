test_that("Example 13.2", {
  tmp <- oddsratio.2x2diff(ch13$physical.activity, ch13$tv.viewing,
                           ch13$gender, do.exact = FALSE,
                           do.asymp = TRUE)
  expect_equal(tmp$pval.asymp.stat, 0.7918946)
  expect_equal(tmp$pval.asymp, 0.42842214)
  expect_equal(tmp$CI.asymp.lower, -0.201675165)
  expect_equal(tmp$CI.asymp.upper, 0.475127192)
  tmp <- oddsratio.2x2diff(ch13$physical.activity, ch13$tv.viewing,
                            ch13$gender)
  expect_equal(tmp$pval.exact.stat, 0.136726013)
  expect_equal(tmp$pval.exact, 0.44031545)
  expect_equal(tmp$CI.exact.lower, -0.34815705)
  expect_equal(tmp$CI.exact.upper, 0.34589494)
})
