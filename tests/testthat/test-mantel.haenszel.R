test_that("Example 13.4", {
  expect_equal(mantel.haenszel(ch13$drug, ch13$side.effects, ch13$age.group)$pval.asymp.stat,
               0.0051217907)
  expect_equal(mantel.haenszel(ch13$drug, ch13$side.effects, ch13$age.group)$pval.asymp,
               0.94294675)
})
