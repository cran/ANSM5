test_that("Example 6.9", {
  expect_equal(control.median(ch6$sampleI, ch6$sampleII,
                              alternative = "greater")$pval.exact, 0.0077527593)
})

test_that("Exercise 9.8", {
  expect_equal(control.median(ch9$bulbA, ch9$bulbB,
                              alternative = "greater")$pval.exact, 0.014207735)
})
