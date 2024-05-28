test_that("Example 13.3", {
  expect_equal(breslow.day(ch13$machine, ch13$output.status,
                           ch13$material.source)$pval.asymp.stat, 0.089150303)
  expect_equal(breslow.day(ch13$machine, ch13$output.status,
                           ch13$material.source)$pval.asymp, 0.95640372)
  expect_equal(breslow.day(ch13$machine, ch13$output.status,
                           ch13$material.source)$CI.asymp.lower, 2.02585526)
  expect_equal(breslow.day(ch13$machine, ch13$output.status,
                           ch13$material.source)$CI.asymp.upper, 10.913474)
})

test_that("Example 13.4", {
  expect_equal(breslow.day(ch13$drug, ch13$side.effects,
                           ch13$age.group)$pval.asymp.stat, 6.4969442)
  expect_equal(breslow.day(ch13$drug, ch13$side.effects,
                           ch13$age.group)$pval.asymp, 0.089783097)
  expect_equal(breslow.day(ch13$drug, ch13$side.effects,
                           ch13$age.group)$CI.asymp.lower, 0.40369219)
  expect_equal(breslow.day(ch13$drug, ch13$side.effects,
                           ch13$age.group)$CI.asymp.upper,2.64532585)
})

test_that("Exercise 13.7", {
  expect_equal(breslow.day(ch13$medicine, ch13$response,
                           ch13$location)$pval.asymp, 0.52001332)
})
