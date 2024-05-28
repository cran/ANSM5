test_that("Example 6.16", {
  expect_equal(cramer.von.mises(ch6$salivaF, ch6$salivaM)$pval.stat,
               0.33418367)
  expect_equal(cramer.von.mises(ch6$salivaF, ch6$salivaM)$pval,
               "p-value >= 0.10")
  expect_equal(cramer.von.mises(ch6$salivaF, ch6$salivaM,
                                alternative = "greater")$pval,
               "p-value >= 0.05")
})
