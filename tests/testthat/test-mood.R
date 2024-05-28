test_that("Example 6.12", {
  tmp <- mood(ch6$typeA, ch6$typeB)
  expect_equal(tmp$pval.exact.stat, 164)
  expect_equal(tmp$pval.exact, 0.12933238)
  expect_equal(mood(ch6$typeA, ch6$typeB, do.exact = FALSE,
                    do.asymp = TRUE)$pval.asymp.stat, 164)
  expect_equal(mood(ch6$typeA, ch6$typeB, do.exact = FALSE,
                    do.asymp = TRUE)$pval.asymp, 0.12562273)
})
