test_that("Example 4.6", {
  expect_equal(binom(3, 20)$CI.exact.lower, 0.03207)
  expect_equal(binom(3, 20)$CI.exact.upper, 0.37893)
  expect_equal(binom(3, 20)$actualCIwidth.exact, 0.95000353)
})

test_that("Example 4.8", {
  expect_equal(pbinom(4, 20, 0.1), 0.9568255)
  expect_equal(pbinom(4, 20, 0.1, lower.tail = FALSE), 0.043174495)
  expect_equal(pbinom(4, 20, 0.2, lower.tail = FALSE), 0.37035174)
  expect_equal(pnorm(0.1118, lower.tail = FALSE), 0.455491)
})

test_that("Example 4.10", {
  expect_equal(binom(1, 54, 0.07, alternative="less")$pval.exact, 0.100604448)
  expect_equal(binom(1, 54)$CI.exact.lower, 0.00046)
  expect_equal(binom(1, 54)$CI.exact.upper, 0.09892)
})

test_that("Example 4.11", {
  expect_equal(binom(15, 20, 0.5)$pval.exact, 0.041389465)
  expect_equal(binom(15, 20)$CI.exact.lower, 0.50895)
  expect_equal(binom(15, 20)$CI.exact.upper, 0.91343)
  expect_equal(binom(14, 20)$CI.exact.lower, 0.45721)
  expect_equal(binom(14, 20)$CI.exact.upper, 0.88107)
})

test_that("Example 4.12", {
  expect_equal(binom(13, 32, 0.25)$pval.exact, 0.062912758)
  expect_equal(binom(13, 32)$CI.exact.lower, 0.23698)
  expect_equal(binom(13, 32)$CI.exact.upper, 0.59356)
})

test_that("Exercise 4.6", {
  expect_equal(binom(26, 72)$CI.exact.lower, 0.25116)
  expect_equal(binom(26, 72)$CI.exact.upper, 0.48288)
})

test_that("Exercise 4.7", {
  expect_equal(binom(6, 18, 0.5)$pval.exact, 0.23788452)
  expect_equal(binom(75, 225, 0.5)$pval.exact, 6.433102e-07)
  expect_equal(binom(6, 18)$CI.exact.lower, 0.13342)
  expect_equal(binom(6, 18)$CI.exact.upper, 0.59008)
  expect_equal(binom(75, 225)$CI.exact.lower, 0.27208)
  expect_equal(binom(75, 225)$CI.exact.upper, 0.39905)
})

test_that("Example 5.4", {
  expect_equal(binom(14, 23, 0.5)$pval.exact, 0.4048729)
})

test_that("Exercise 5.7", {
  expect_equal(binom(32, 57, 0.5)$pval.exact, 0.42704286)
  expect_equal(binom(32, 70, 0.75, alternative = "less")$pval.exact,
               0.000000158118703)
  expect_equal(binom(32, 70, 0.75)$CI.exact.lower, 0.33744)
  expect_equal(binom(32, 70, 0.75)$CI.exact.upper, 0.58058)
})

test_that("Exercise 5.8", {
  expect_equal(binom(24, 40, 0.5)$pval.exact, 0.26818725)
})

test_that("Section 13.4", {
  expect_equal(binom.test(32, 100, 0.4)$p.value, 0.125293086)
  expect_equal(binom.test(48, 100, 0.4)$p.value, 0.10363706)
})
