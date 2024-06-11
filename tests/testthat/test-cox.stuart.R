test_that("Example 4.13", {
  expect_equal(cox.stuart(ch4$precipitation)$pval.exact, 0.68750)
})

test_that("Exercise 4.8", {
  expect_equal(cox.stuart(ch4$days.waiting)$pval.exact, 0.0703125)
})

test_that("Exercise 4.9", {
  expect_equal(cox.stuart(ch4$rainfall.by.latitude)$pval.exact, 0.125)
})

test_that("Exercise 4.10", {
  expect_equal(cox.stuart(ch4$points)$pval.exact, 0.2890625)
})

test_that("Exercise 4.11", {
  expect_equal(cox.stuart(ch4$rainfall.DRC)$pval.exact, 0.2890625)
})

test_that("Exercise 4.12", {
  expect_equal(cox.stuart(ch4$piped.water.DRC)$pval.exact, 1)
})

test_that("Exercise 10.5", {
  expect_equal(cox.stuart(app1$McDelta[order(ch10$death.year)], alternative = "less")$pval.exact,
               0.1875)
})
