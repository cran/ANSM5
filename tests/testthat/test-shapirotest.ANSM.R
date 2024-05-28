test_that("Example 4.4", {
  expect_equal(shapirotest.ANSM(ch4$ages)$pval.asymp, 0.00117134809)
})

test_that("Exercise 4.1", {
  expect_equal(shapirotest.ANSM(ch4$McGamma)$pval.asymp, 0.00039674968)
})

test_that("Exercise 4.3", {
  expect_equal(shapirotest.ANSM(ch4$visiting.supporters)$pval.asymp,
               0.0011657189)
})

test_that("Example 5.1", {
  expect_equal(shapirotest.ANSM(ch5$LVF - ch5$RVF)$pval.asymp, 0.7289256)
})

test_that("Example 5.3", {
  expect_equal(shapirotest.ANSM(ch5$bp)$pval.asymp.stat, 0.90126828)
  expect_equal(shapirotest.ANSM(ch5$bp)$pval.asymp, 0.022900098)
  expect_equal(shapirotest.ANSM(ch5$bp.incorrect)$pval.asymp.stat,
               0.79465777)
  expect_equal(shapirotest.ANSM(ch5$bp.incorrect)$pval.asymp,
               0.000241231707)
})
