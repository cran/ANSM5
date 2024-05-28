test_that("Example 6.14", {
  tmp <- moses.extreme.reactions(ch6$groupI, ch6$groupII)
  expect_equal(tmp$pval.exact.stat, 12)
  expect_equal(tmp$pval.exact, 0.0148844963)
  tmp <- moses.extreme.reactions(ch6$groupI.trimmed, ch6$groupII)
  expect_equal(tmp$pval.exact, 0.036404772)
  expect_equal(tmp$pval.exact.stat, 10)
  expect_equal(
    moses.extreme.reactions(ch6$groupI.amended, ch6$groupII)$pval.exact,
    0.03989045)
})
