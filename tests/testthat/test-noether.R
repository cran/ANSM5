test_that("Example 5.8", {
  #  expect_equal(VGAM::plaplace(0, 1, sqrt(2), lower.tail = FALSE), 0.75346565)
  expect_equal(noether(p1 = 0.7534, alpha = 0.05, power = 0.9), 33.342312)
  #
  expect_equal(pnorm(-0.5, lower.tail = FALSE), 0.69146246)
  expect_equal(noether(p1 = 0.6915, alpha = 0.05, power = 0.9), 58.380978)
  #
  expect_equal(pexp(0.3862, rate = 1/2, lower.tail = FALSE), 0.82439953)
  expect_equal(noether(p1 = 0.8244, alpha = 0.05, power = 0.9), 20.34451)
})

test_that("Example 5.9", {
  expect_equal(pnorm(-1 / sqrt(2), lower.tail = FALSE), 0.76024994)
  expect_equal(noether(p1 = 0.7602, alpha = 0.05, power = 0.8), 22.829351)
})

test_that("Exercise 5.16", {
  expect_equal(noether(p1 = 0.8, alpha = 0.05, power = 0.9), 23.788465)
  expect_equal(noether(p1 = 0.8, alpha = 0.025, power = 0.9), 29.187286)
})
