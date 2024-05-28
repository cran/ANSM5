test_that("Example 9.6", {
  tmp <- logrank(ch9$samplesAB.survtime, ch9$samplesAB.censor, ch9$samplesAB,
                 score.censored = FALSE)
  expect_equal(tmp$pval.exact.stat, 1.10595238)
  expect_equal(tmp$pval.exact, 0.33809524)
})

test_that("Example 9.7", {
  tmp <- logrank(ch9$samplesXYZ.survtime, ch9$samplesXYZ.censor, ch9$samplesXYZ)
  expect_equal(tmp$pval.exact.stat, 1.3702381)
  expect_equal(tmp$pval.exact, 0.29761905)
})
