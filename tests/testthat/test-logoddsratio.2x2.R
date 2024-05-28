test_that("Exercise 13.2", {
  expect_equal(logoddsratio.2x2(ch13$physical.activity[ch13$gender == "Boy"],
                                  ch13$tv.viewing[ch13$gender == "Boy"],
                                  do.exact = FALSE,
                                  do.asymp = TRUE)$pval.asymp.stat, 0.23135076)
  expect_equal(logoddsratio.2x2(ch13$physical.activity[ch13$gender == "Boy"],
                                  ch13$tv.viewing[ch13$gender == "Boy"],
                                  do.exact = FALSE,
                                  do.asymp = TRUE)$pval.asymp, 0.040533743)
  expect_equal(logoddsratio.2x2(ch13$physical.activity[ch13$gender == "Girl"],
                                  ch13$tv.viewing[ch13$gender == "Girl"],
                                  do.exact = FALSE,
                                  do.asymp = TRUE)$pval.asymp.stat, 0.36807677)
  expect_equal(logoddsratio.2x2(ch13$physical.activity[ch13$gender == "Girl"],
                                  ch13$tv.viewing[ch13$gender == "Girl"],
                                  do.exact = FALSE,
                                  do.asymp = TRUE)$pval.asymp, 0.0048226686)
})
