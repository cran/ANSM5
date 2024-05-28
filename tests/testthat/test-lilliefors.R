test_that("Example 4.4", {
  expect_equal(lilliefors(ch4$ages, seed = 1)$pval.mc, 0.0057)
})
