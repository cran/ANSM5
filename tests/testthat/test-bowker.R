test_that("Example 12.12", {
  expect_equal(bowker(ch12$side.effect.new,
                      ch12$side.effect.old)$pval.asymp.stat, 9.3333333)
  expect_equal(bowker(ch12$side.effect.new,
                      ch12$side.effect.old)$pval.asymp, 0.02517214)
})

test_that("Exercise 12.12", {
  expect_equal(bowker(ch12$first.response,
                      ch12$second.response)$pval.asymp.stat, 23)
  expect_equal(bowker(ch12$first.response,
                      ch12$second.response)$pval.asymp, 0.000001620014)
})
