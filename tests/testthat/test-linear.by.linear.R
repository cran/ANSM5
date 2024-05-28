test_that("Example 13.8", {
  expect_equal(linear.by.linear(
    ch13$dose, ch13$dose.side.effect, do.mc = FALSE,
    do.asymp = TRUE)$pval.asymp.stat, 2.3888324)
  expect_equal(linear.by.linear(
    ch13$dose, ch13$dose.side.effect, do.mc = FALSE,
    do.asymp = TRUE)$pval.asymp, 0.0084510055)
})

test_that("Example 13.10", {
  expect_equal(linear.by.linear(ch13$frequency, ch13$person, do.mc = FALSE,
                                do.asymp = TRUE)$pval.asymp * 2, 0.0051610501)
})

test_that("Exercise 13.4", {
  expect_equal(linear.by.linear(ch13$SBP, ch13$cholesterol, do.mc = FALSE,
                                do.asymp = TRUE)$pval.asymp * 2, 0.58217184)
  expect_equal(linear.by.linear(ch13$SBP, ch13$cholesterol,
                   u = c(1, 2, 3, 6, 10), v = c(1, 2, 6), do.mc = FALSE,
                   do.asymp = TRUE)$pval.asymp * 2, 0.381711194)
})

test_that("Exercise 13.5", {
  expect_equal(linear.by.linear(ch13$schooling, ch13$abortion.attitude,
                                do.mc = FALSE, do.asymp = TRUE)$pval.asymp * 2,
               1.339675e-21)
})

test_that("Exercise 13.9", {
  expect_equal(linear.by.linear(ch13$PPI.ages, ch13$PPI.people,
                                do.mc = FALSE, do.asymp = TRUE)$pval.asymp * 2,
               8.1210037e-08)
})
