test_that("Example 4.14", {
  tmp <- runs.2cat(ch4$tosses1)
  expect_equal(tmp$pval.exact.stat, 3)
  expect_equal(tmp$pval.exact, 0.000216501764)
  expect_equal(runs.2cat(ch4$tosses1, do.exact = FALSE,
                         do.asymp = TRUE)$pval.asymp.stat, -3.4460122)
  expect_equal(runs.2cat(ch4$tosses1, do.exact = FALSE,
                         do.asymp = TRUE)$pval.asymp, 0.0005689247)
  tmp <- runs.2cat(ch4$tosses2)
  expect_equal(tmp$pval.exact.stat, 20)
  expect_equal(tmp$pval.exact, 0.0000216501764)
  expect_equal(runs.2cat(ch4$tosses2, do.exact = FALSE,
                         do.asymp = TRUE)$pval.asymp.stat, 3.9054805)
  expect_equal(runs.2cat(ch4$tosses2, do.exact = FALSE,
                         do.asymp = TRUE)$pval.asymp, 0.00009403835)
  tmp <- runs.2cat(ch4$tosses3)
  expect_equal(tmp$pval.exact.stat, 13)
  expect_equal(tmp$pval.exact, 0.36184806)
  expect_equal(runs.2cat(ch4$tosses3, do.exact = FALSE,
                         do.asymp = TRUE)$pval.asymp.stat, 0.74299208)
  expect_equal(runs.2cat(ch4$tosses3, do.exact = FALSE,
                         do.asymp = TRUE)$pval.asymp, 0.45748648)
})

test_that("Exercise 4.10", {
  expect_equal(runs.2cat(
    sign(ch4$points - median(ch4$points)))$pval.exact, 0.063403263)
})

test_that("Exercise 4.11", {
  expect_equal(runs.2cat(
    sign(ch4$rainfall.DRC - median(ch4$rainfall.DRC)))$pval.exact,
    0.063403263)
})

test_that("Exercise 4.12", {
  expect_equal(runs.2cat(
    sign(ch4$piped.water.DRC - median(ch4$piped.water.DRC)))$pval.exact,
    0.42890443)
})

test_that("Example 5.10", {
  expect_equal(runs.2cat(sign(ch5$yr0910 - ch5$yr1314),
                         alternative = "less")$pval.exact, 0.65151515)
})

test_that("Example 6.17", {
  code <- c(rep(0, length(ch6$groupA)), rep(1, length(ch6$groupB)))
  ranked_code <- code[order(c(ch6$groupA, ch6$groupB))]
  tmp <- runs.2cat(ranked_code, alternative = "less")
  expect_equal(tmp$pval.exact.stat, 12)
  expect_equal(tmp$pval.exact, 0.68004287)
})

test_that("Example 6.18", {
  tmp <- runs.2cat(ch6$sex, alternative = "less")
  expect_equal(tmp$pval.exact.stat, 5)
  expect_equal(tmp$pval.exact, 0.196969697)
})

test_that("Exercise 6.17", {
  expect_equal(runs.2cat(ch6$twins, alternative = "greater")$pval.exact,
               0.66408669)
})
