test_that("Example 11.4", {
  tmp <- theil.kendall(ch11$reportedtime, ch11$parentlimit, H0 = 1,
                       do.asymp = TRUE)
  expect_equal(tmp$stat, 0.56666667)
  expect_equal(tmp$pval.exact, 0.081746032)
  expect_equal(tmp$pval.asymp, 0.176474291)
})

test_that("Example 11.6", {
  expect_equal(theil.kendall(ch11$reportedtime, ch11$parentlimit,
                             do.alpha = TRUE)$stat.note,
               paste0("Estimate of alpha using median of d_i: 2.33333\n",
                      "Hodges-Lehmann estimator of alpha: 2.40000\n\n",
                      "Neither exact, asymptotic nor Monte Carlo test requested"))
})

test_that("Exercise 11.3", {
  expect_equal(theil.kendall(ch11$rotten, ch11$days.stored)$stat, 2.375)
})

test_that("Exercise 11.9", {
  expect_equal(theil.kendall(ch11$weight.gain.A, ch11$food.weight.A,
                             do.abbreviated = TRUE)$stat, 0.17105263)
  expect_equal(theil.kendall(ch11$weight.gain.B, ch11$food.weight.B,
                             do.abbreviated = TRUE)$stat, 0.35365854)
})

test_that("Exercise 11.10", {
  expect_equal(theil.kendall(ch11$SW.England, ch11$N.Scotland)$stat, 1.024)
  expect_equal(theil.kendall(ch11$N.Scotland, ch11$SW.England)$stat, 0.97196262)
  expect_equal(median(
    ch11$N.Scotland - theil.kendall(ch11$N.Scotland, ch11$SW.England)$stat *
      ch11$SW.England), 20.490654)
  expect_equal(median(
    ch11$N.Scotland - theil.kendall(ch11$N.Scotland, ch11$SW.England)$stat *
      ch11$SW.England) + theil.kendall(ch11$N.Scotland, ch11$SW.England)$stat *
      122, 139.070093)
})
