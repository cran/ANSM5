test_that("Example 6.13", {
  tmp <- conover(ch6$typeA, ch6$typeB)
  expect_equal(tmp$pval.exact.stat, 720)
  expect_equal(tmp$pval.exact, 0.126728586)
  expect_equal(conover(ch6$typeA, ch6$typeB, do.exact = FALSE,
                       do.asymp = TRUE)$pval.asymp.stat, 720)
  expect_equal(conover(ch6$typeA, ch6$typeB, do.exact = FALSE,
                       do.asymp = TRUE)$pval.asymp, 0.110647884)
})

test_that("Exercise 6.12", {
  expect_equal(conover(ch6$typeA, ch6$typeB,
                       abs.ranks = TRUE)$pval.exact, 0.13493689)
  expect_equal(conover(ch6$typeA, ch6$typeB, abs.ranks = TRUE, do.exact = FALSE,
                       do.asymp = TRUE)$pval.asymp, 0.10703901)
})

test_that("Exercise 6.14", {
  expect_equal(conover(ch6$time.withoutLD,
                       ch6$time.withLD)$pval.exact, 0.55602176)
})

test_that("Exercise 6.16", {
  expect_equal(conover(ch6$travel, ch6$politics)$pval.exact,
               0.0104786854)
})

test_that("Example 7.9", {
  expect_equal(conover(ch7$dementia.age, ch7$features, do.exact = FALSE,
                       do.asymp = TRUE)$pval.asymp.stat, 1.85924)
  expect_equal(conover(ch7$dementia.age, ch7$features, do.exact = FALSE,
                       do.asymp = TRUE)$pval.asymp, 0.39470367)
})

test_that("Exercise 7.15", {
  expect_equal(conover(ch7$braking.distance.2, ch7$initial.speed,
                       do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.0455638507)
})
