test_that("Section 13.2.5", {
  tmp <- zelen(ch13$drug, ch13$side.effects, ch13$age.group)
  expect_equal(tmp$pval.exact.stat, 0.0082126222)
  expect_equal(tmp$pval.exact, 0.068931048)
})

test_that("Example 13.3", {
  tmp <- zelen(ch13$machine, ch13$output.status, ch13$material.source)
  expect_equal(tmp$pval.exact.stat, 0.14161837)
  expect_equal(tmp$pval.exact, 1)
})
