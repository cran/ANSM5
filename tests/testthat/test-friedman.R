test_that("Example 7.6", {
  tmp <- friedman(ch7$pulse, ch7$time.period, ch7$student)
  expect_equal(tmp$pval.exact.stat, 10.5714286)
  expect_equal(tmp$pval.exact, 0.00272205075)
  expect_equal(friedman(ch7$pulse, ch7$time.period, ch7$student,
                        do.exact = FALSE, do.asymp = TRUE)$pval.asymp.stat,
               10.5714286)
  expect_equal(friedman(ch7$pulse, ch7$time.period, ch7$student,
                        do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.0050634142)
})

test_that("Example 7.7", {
  expect_equal(friedman(ch7$nodes, ch7$treatment, ch7$block,
                        use.Iman.Davenport = TRUE, do.exact = FALSE,
                        do.asymp = TRUE)$pval.asymp.stat, 6.4035088)
  expect_equal(friedman(ch7$nodes, ch7$treatment, ch7$block,
                        use.Iman.Davenport = TRUE, do.exact = FALSE,
                        do.asymp = TRUE)$pval.asymp, 0.0022524635)
})

test_that("Exercise 7.9", {
  expect_equal(friedman(ch7$prem.contractions, ch7$drug,
                        ch7$patient, do.exact = FALSE,
                        do.asymp = TRUE)$pval.asymp, 0.0179300609)
  expect_equal(friedman(ch7$prem.contractions, ch7$drug,
                        ch7$patient, use.Iman.Davenport = TRUE,
                        do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.0112272605)
})

test_that("Exercise 7.11", {
  #(i)
  expect_equal(friedman(ch7$births, ch7$weekday, ch7$week,
                        do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.26163381)
  expect_equal(friedman(ch7$births, ch7$weekday, ch7$week,
                        use.Iman.Davenport = TRUE, do.exact = FALSE,
                        do.asymp = TRUE)$pval.asymp, 0.26281532)
  #(ii)
  expect_equal(friedman(ch7$births, ch7$week, ch7$weekday,
                        do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.019992448)
  expect_equal(friedman(ch7$births, ch7$week, ch7$weekday,
                        use.Iman.Davenport = TRUE, do.exact = FALSE,
                        do.asymp = TRUE)$pval.asymp, 0.0086079546)
})

test_that("Exercise 7.12", {
  expect_equal(friedman(ch7$names.recalled, ch7$group,
                        ch7$medical.student, do.exact = FALSE,
                        do.asymp = TRUE)$pval.asymp, 0.0043202395)
  expect_equal(friedman(ch7$names.recalled, ch7$group,
                        ch7$medical.student, use.Iman.Davenport = TRUE,
                        do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.00084503)
})
