test_that("Example 4.15", {
  expect_equal(runs.ncat(ch4$births, alternative = "less")$pval.asymp.stat,
               -1.15433076)
  expect_equal(runs.ncat(ch4$births, alternative = "less",
                         cont.corr = FALSE)$pval.asymp, 0.061889005)
  expect_equal(runs.ncat(ch4$births, alternative = "less")$pval.asymp,
               0.124182295)
})

test_that("Exercise 4.17", {
  expect_equal(runs.ncat(ch4$football.results, alternative = "less",
                         cont.corr = FALSE)$pval.asymp, 0.017170837)
  expect_equal(runs.ncat(ch4$football.results,
                         alternative = "less")$pval.asymp, 0.03603851)
})

test_that("Example 7.10", {
  expect_equal(runs.ncat(ch7$team.member[order(ch7$procedure.time)],
                         alternative = "less",
                         cont.corr = FALSE)$pval.asymp.stat, -4.2602817)
  expect_equal(runs.ncat(ch7$team.member[order(ch7$procedure.time)], alternative = "less",
                         cont.corr = FALSE)$pval.asymp, 0.0000102084714)
  expect_equal(runs.ncat(ch7$team.member[order(ch7$procedure.time)], alternative = "less",
                         cont.corr = TRUE)$pval.asymp.stat, -4.0666325)
  expect_equal(runs.ncat(ch7$team.member[order(ch7$procedure.time)], alternative = "less",
                         cont.corr = TRUE)$pval.asymp, 0.000023848684)
})

test_that("Exercise 7.16", {
  expect_equal(runs.ncat(ch7$regions[order(ch7$affordability)],
                         alternative = "less")$pval.asymp.stat, -1.41710858)
  expect_equal(runs.ncat(ch7$regions[order(ch7$affordability)],
                         alternative = "less")$pval.asymp, 0.078225593)
})
