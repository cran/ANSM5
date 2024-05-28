test_that("Example 10.8", {
  expect_equal(kendall.tau(ch10$death.year, ch10$age.at.death,
                           alternative = "greater", do.asymp = TRUE,
                           do.exact = FALSE)$stat, 0.38964325)
  expect_equal(kendall.tau(ch10$death.year, ch10$age.at.death,
                           alternative = "greater", do.asymp = TRUE,
                           do.exact = FALSE)$pval.asymp, 0.03309629)
})

test_that("Example 10.8", {
  tmp <- kendall.tau(ch10$British, ch10$American)
  expect_equal(tmp$stat, 0.61904762)
  expect_equal(tmp$pval.exact, 0.06904762)
})

test_that("Example 10.9", {
  tmp <- kendall.tau(ch10$Canadian, ch10$Australian)
  expect_equal(tmp$stat, 0.047619048)
  expect_equal(tmp$pval.exact, 1)
})
