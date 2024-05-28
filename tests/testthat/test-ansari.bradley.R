test_that("Example 6.12", {
  tmp <- ansari.bradley(ch6$typeA, ch6$typeB)
  expect_equal(tmp$pval.exact.stat, 30)
  expect_equal(tmp$pval.exact,
               0.143732635)
  expect_equal(ansari.bradley(ch6$typeA, ch6$typeB, do.exact = FALSE,
                              do.asymp = TRUE)$pval.asymp.stat, 30)
  expect_equal(ansari.bradley(ch6$typeA, ch6$typeB, do.exact = FALSE,
                              do.asymp = TRUE)$pval.asymp, 0.121413624)
})

test_that("Exercise 6.14", {
  expect_equal(ansari.bradley(ch6$time.withoutLD,
                              ch6$time.withLD)$pval.exact, 1)
})

test_that("Exercise 6.16", {
  expect_equal(ansari.bradley(ch6$travel, ch6$politics)$pval.exact,
               0.032810842)
})
