test_that("Example 7.1", {
  tmp <- kruskal.wallis(ch7$affordability, ch7$regions)
  expect_equal(tmp$pval.exact.stat, 6.7454545)
  expect_equal(tmp$pval.exact, 0.01)
  expect_equal(
    kruskal.wallis(ch7$affordability, ch7$regions, do.exact = FALSE,
                   do.asymp = TRUE)$pval.asymp.stat, 6.7454545)
  expect_equal(
    kruskal.wallis(ch7$affordability, ch7$regions, do.exact = FALSE,
                   do.asymp = TRUE)$pval.asymp, 0.034295975)
})

test_that("Example 7.2", {
  tmp <- kruskal.wallis(ch7$age, ch7$positions)
  expect_equal(tmp$pval.exact.stat, 2.734127)
  expect_equal(tmp$pval.exact, 0.270820846)
  expect_equal(kruskal.wallis(ch7$age, ch7$positions, do.exact = FALSE,
                              do.asymp = TRUE)$pval.asymp.stat, 2.734127)
  expect_equal(kruskal.wallis(ch7$age, ch7$positions, do.exact = FALSE,
                              do.asymp = TRUE)$pval.asymp, 0.25485424)
})

test_that("Exercise 7.2", {
  expect_equal(kruskal.wallis(ch7$sentences, ch7$author,
                              do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.010327934)
})

test_that("Exercise 7.4", {
  expect_equal(kruskal.wallis(ch7$head.width, ch7$species,
                              do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.10882766)
})

test_that("Exercise 7.16", {
  expect_equal(kruskal.wallis(ch7$affordability,
                              ch7$regions)$pval.exact, 0.01)
})

test_that("Example 12.6", {
  tmp <- kruskal.wallis(ch12$time.to.failure, ch12$cause, do.asymp = TRUE)
  expect_equal(tmp$pval.asymp.stat, 10.444765)
  expect_equal(tmp$pval.asymp, 0.033566401)
})
