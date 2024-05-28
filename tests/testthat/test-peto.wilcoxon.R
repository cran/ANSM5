test_that("Example 9.4", {
  tmp <- peto.wilcoxon(ch9$sampleI.survtime, ch9$sampleII.survtime,
                       ch9$sampleI.censor, ch9$sampleII.censor,
                       alternative = "less")
  expect_equal(tmp$pval.exact.stat, 66.25)
  expect_equal(tmp$pval.exact, 0.048174048)
})
