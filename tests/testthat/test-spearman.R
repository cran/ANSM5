test_that("Example 10.2", {
  expect_equal(spearman(ch10$q1, ch10$q2, alternative = "greater",
                        do.asymp = TRUE, do.exact = FALSE)$stat, 0.43356643)
  expect_equal(spearman(ch10$q1, ch10$q2, alternative = "greater",
                        do.asymp = TRUE, do.exact = FALSE)$pval.asymp,
               0.079552885)
})

test_that("Exercise 10.1", {
  expect_equal(spearman(ch10$ERA, ch10$ESMS, do.exact = FALSE)$stat, 0.9030303)
})

test_that("Example 10.8", {
  tmp <- spearman(ch10$British, ch10$American)
  expect_equal(tmp$stat, 0.82142857)
  expect_equal(tmp$pval.exact, 0.034126984)
})

test_that("Example 10.9", {
  tmp <- spearman(ch10$Canadian, ch10$Australian)
  expect_equal(tmp$stat, 0)
  expect_equal(tmp$pval.exact, 1)
})
