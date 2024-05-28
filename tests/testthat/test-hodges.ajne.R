test_that("Example 4.16", {
  expect_equal(hodges.ajne(ch4$times.as.degrees)$pval.stat, 1)
  expect_equal(hodges.ajne(ch4$times.as.degrees)$pval, 0.05859375)
})

test_that("Exercise 4.13", {
  expect_equal(hodges.ajne(ch4$accident.bearings)$pval, 0.109375)
})

test_that("Exercise 4.14", {
  expect_equal(hodges.ajne(ch4$board.angles)$pval, 0.9375)
})

test_that("Test overly-uniform data", {
  uniform_data <- c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330)
  expect_equal(hodges.ajne(uniform_data)$pval.stat, 6)
  expect_equal(hodges.ajne(uniform_data)$pval, NA)
})
