test_that("Example 9.1", {
  tmp <- gehan.wilcoxon(ch9$symp.survtime, ch9$asymp.survtime,
                        ch9$symp.censor, ch9$asymp.censor,
                        alternative = "less", do.exact = FALSE,
                        do.asymp = TRUE)
  expect_equal(tmp$pval.asymp.stat, 217.5)
  expect_equal(tmp$pval.asymp, 0.0050981517)
})
