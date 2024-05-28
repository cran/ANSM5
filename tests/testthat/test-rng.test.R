test_that("Example 4.17", {
  expect_equal(rng.test(ch4$dates.as.degrees)$pval.stat, 1.01229097)
  expect_equal(rng.test(ch4$dates.as.degrees)$pval, 0.000022540754)
})

test_that("Exercise 4.13", {
  expect_equal(rng.test(ch4$accident.bearings)$pval, 0.0074382381)
})
