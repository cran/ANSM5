test_that("Example 6.11", {
  tmp <- siegel.tukey(ch6$typeA, ch6$typeB)
  expect_equal(tmp$pval.exact.stat,
               paste0("\n","106 (rank sum from ch6$typeA), ",
                      "104 (rank sum from ch6$typeB)", "\n",
                      "70 (Mann-Whitney U from ch6$typeA), ",
                      "26 (Mann-Whitney U from ch6$typeB)"))
  expect_equal(tmp$pval.exact, 0.097880448)
})

test_that("Exercise 6.11", {
  expect_equal(siegel.tukey(ch6$typeA, ch6$typeB,
                            mean.shift = TRUE)$pval.exact, 0.20828769)
})

test_that("Exercise 6.14", {
  expect_equal(siegel.tukey(ch6$time.withoutLD,
                            ch6$time.withLD)$pval.exact, 0.86651127)
})

test_that("Exercise 6.16", {
  expect_equal(siegel.tukey(ch6$travel, ch6$politics)$pval.exact,
               0.040821408)
})
