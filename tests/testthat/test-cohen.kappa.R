test_that("Example 10.11", {
  expect_equal(cohen.kappa(ch10$dentistA, ch10$dentistB, do.exact = FALSE,
                           do.asymp = TRUE, alternative = "greater")$stat,
               0.41747573)
  expect_equal(cohen.kappa(ch10$dentistA, ch10$dentistB, do.exact = FALSE,
                           do.asymp = TRUE,
                           alternative = "greater")$pval.asymp.stat,
               4.5303334)
  expect_equal(cohen.kappa(ch10$dentistA, ch10$dentistB, do.exact = FALSE,
                           do.asymp = TRUE,
                           alternative = "greater")$pval.asymp, 0.0000029445342)
})

test_that("Example 10.12", {
  tmp <- cohen.kappa(ch10$questionnaire, ch10$demonstration, ch10$items)
  expect_equal(tmp$stat, c(0.3229219144, 0.2273352191, 0.0187919463))
  expect_equal(tmp$pval.asymp.stat, 7.1914942)
  expect_equal(tmp$pval.asymp, 0.027440175)
})

test_that("Exercise 10.12", {
  expect_equal(cohen.kappa(ch10$observerA, ch10$observerB)$stat, 0.39839572)
  expect_equal(cohen.kappa(ch10$observerA, ch10$observerB, do.exact = FALSE,
                           do.asymp = TRUE)$pval.asymp, 2.528931e-10)
})

test_that("Exercise 10.13", {
  tmp <- cohen.kappa(ch10$questionnaire, ch10$demonstration,
                     ch10$gender)
  expect_equal(tmp$stat, c(0.184024697, 0.245444989))
  expect_equal(tmp$pval.asymp.stat, 0.85621234)
  expect_equal(tmp$pval.asymp, 0.35480082)
})
