test_that("Example 7.8", {
  tmp <- cochran.q(ch7$outcome, ch7$climb, ch7$member)
  expect_equal(tmp$pval.exact.stat, 1.2)
  expect_equal(tmp$pval.exact, 0.85185185)
  expect_equal(cochran.q(ch7$outcome, ch7$climb, ch7$member,
                         do.exact = FALSE, do.asymp = TRUE)$pval.asymp.stat,
               1.2)
  expect_equal(cochran.q(ch7$outcome, ch7$climb, ch7$member,
                         do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.54881164)
})

test_that("Exercise 7.14", {
  expect_equal(cochran.q(ch7$soc.media.use, ch7$participant,
                         ch7$day, do.exact = FALSE,
                         do.asymp = TRUE)$pval.asymp, 0.69740213)
})
