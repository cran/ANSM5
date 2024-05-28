test_that("Example 14.2", {
  bs.sample <- bs(ch14$example14.2, nsims.bs = 40, seed = 1)
  expect_equal(range(bs.sample$vector), c(0.805, 3.5))
  expect_equal(sd(bs.sample$vector), 0.46151393)
})

test_that("Example 14.5", {
  expect_equal(bs(ch14$example14.2, nsims.bs = 2000, CI.width = 0.95,
                  seed = 1)$CI, c(1.005, 2.74))
  expect_equal(bs(ch14$example14.2, nsims.bs = 2000, CI.width = 0.99,
                  seed = 1)$CI, c(0.69, 3.5))
  expect_equal(sd(bs(ch14$example14.2, nsims.bs = 2000, seed = 1)$vector),
               0.4144401)
})
