test_that("Example 12.1", {
  expect_equal(chisqtest.ANSM(ch12$feedback.freq, ch12$PPI.person,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 9.9061084)
  expect_equal(chisqtest.ANSM(ch12$feedback.freq, ch12$PPI.person,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 0.0070618077)
})

test_that("Example 12.3", {
  expect_equal(chisqtest.ANSM(ch12$infection.site, ch12$district, do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 14.959596)
  expect_equal(chisqtest.ANSM(ch12$infection.site, ch12$district, do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 0.133541873)
})

test_that("Example 12.4", {
  expect_equal(chisqtest.ANSM(ch12$drugYZ, ch12$side.effect, do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 11.160998)
  expect_equal(chisqtest.ANSM(ch12$drugYZ, ch12$side.effect, do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 0.0008353497)
  expect_equal(chisqtest.ANSM(ch12$drugYZ, ch12$side.effect, cont.corr = FALSE,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 14.9466695)
})

test_that("Example 12.5", {
  expect_equal(chisqtest.ANSM(ch12$drugAB, ch12$side.effect.level,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 11.7558974)
  expect_equal(chisqtest.ANSM(ch12$drugAB, ch12$side.effect.level,
                               do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.0192618608)
})

test_that("Example 12.10", {
  expect_equal(chisqtest.ANSM(table(ch12$last.digits), do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 11.6324786)
  expect_equal(chisqtest.ANSM(table(ch12$last.digits), do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 0.23483005)
})

test_that("Example 12.11", {
  lambda <- sum(0:6 * table(ch12$accidents)) / length(ch12$accidents)
  poisson.probs <- c(dpois(0:2, lambda), 1 - sum(dpois(0:2, lambda)))
  expect_equal(chisqtest.ANSM(table(ch12$accidents.reduced), p = poisson.probs,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 197.73701)
  expect_equal(chisqtest.ANSM(table(ch12$accidents.reduced), p = poisson.probs,
                               do.exact = FALSE, do.asymp = TRUE)$pval.asymp, 0)
})

test_that("Exercise 12.1", {
  expect_equal(chisqtest.ANSM(ch12$bronchitis, ch12$otitis.media,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 33.263889)
  expect_equal(chisqtest.ANSM(ch12$bronchitis, ch12$otitis.media,
                               do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.00000105472323)
})

test_that("Exercise 12.2", {
  expect_equal(chisqtest.ANSM(ch12$welsh.language, ch12$opportunities,
                               cont.corr = FALSE, do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 5.1437532)
  expect_equal(chisqtest.ANSM(ch12$welsh.language, ch12$opportunities,
                               cont.corr = FALSE, do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 0.0233301574)
  expect_equal(chisqtest.ANSM(ch12$welsh.language, ch12$opportunities,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 3.7356353)
  expect_equal(chisqtest.ANSM(ch12$welsh.language, ch12$opportunities,
                               do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.053263409)
})

test_that("Exercise 12.6", {
  expect_equal(chisqtest.ANSM(ch12$ethnic.group, ch12$diabetes.status,
                               do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.07932863)
})


test_that("Exercise 12.7", {
  expect_equal(chisqtest.ANSM(table(ch12$horse.wins), do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 16)
  expect_equal(chisqtest.ANSM(table(ch12$horse.wins), do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 0.0251163607)
})

test_that("Exercise 12.8", {
  expect_equal(chisqtest.ANSM(table(ch12$F1.wins), p = c(0.05, 0.15, 0.2, 0.6),
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 408.19048)
  expect_equal(chisqtest.ANSM(table(ch12$F1.wins), p = c(0.05, 0.15, 0.2, 0.6),
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 3.723882e-88)
})

test_that("Exercise 12.9", {
  expect_equal(chisqtest.ANSM(table(ch12$strokes), do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 24.2712504)
  expect_equal(chisqtest.ANSM(table(ch12$strokes), do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp, 0.0004655284)
})

test_that("Exercise 12.10", {
  poisson.probs <- c(dpois(0:3, 0.61), 1 - sum(dpois(0:3, 0.61)))
  expect_equal(chisqtest.ANSM(table(ch12$recurrent.visits), p = poisson.probs,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 15370.6084)
  expect_equal(chisqtest.ANSM(table(ch12$recurrent.visits), p = poisson.probs,
                               do.exact = FALSE, do.asymp = TRUE)$pval.asymp, 0)
})

test_that("Exercise 12.11", {
  binomial.probs <- c(dbinom(0:2, 5, 0.188), sum(dbinom(3:5, 5, 0.188)))
  expect_equal(chisqtest.ANSM(table(ch12$holes), p = binomial.probs,
                               do.exact = FALSE,
                               do.asymp = TRUE)$pval.asymp.stat, 5.0729944)
  expect_equal(chisqtest.ANSM(table(ch12$holes), p = binomial.probs,
                               do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
               0.1665295)
})

test_that("Exercise 13.3", {
  expect_equal(
    chisqtest.ANSM(ch13$breakfast.eaten[ch13$boys.girls == "Boys"],
                    ch13$VEL[ch13$boys.girls == "Boys"], do.exact = FALSE,
                    do.asymp = TRUE)$pval.asymp, 0.70749547)
  expect_equal(
    chisqtest.ANSM(ch13$breakfast.eaten[ch13$boys.girls == "Girls"],
                    ch13$VEL[ch13$boys.girls == "Girls"], do.exact = FALSE,
                    do.asymp = TRUE)$pval.asymp, 0.035047454)
})
