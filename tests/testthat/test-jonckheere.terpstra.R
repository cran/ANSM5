test_that("Example 7.3", {
  tmp <- jonckheere.terpstra(ch7$dementia.age, ch7$features,
                             alternative = "greater")
  expect_equal(tmp$pval.exact.stat, 9)
  expect_equal(tmp$pval.exact, 0.0046620047)
  expect_equal(
    jonckheere.terpstra(ch7$dementia.age, ch7$features,
                        alternative = "greater", do.exact = FALSE,
                        do.asymp = TRUE,
                        do.asymp.ties.adjust = FALSE)$pval.asymp.stat,
    -2.52810291)
  expect_equal(
    jonckheere.terpstra(ch7$dementia.age, ch7$features,
                        alternative = "greater", do.exact = FALSE,
                        do.asymp = TRUE,
                        do.asymp.ties.adjust = FALSE)$pval.asymp, 0.0057340365)
  expect_equal(
    jonckheere.terpstra(ch7$dementia.age, ch7$features,
                        alternative = "greater", do.exact = FALSE,
                        do.asymp = TRUE)$pval.asymp, 0.0054947789)
})

test_that("Exercise 7.5", {
  tmp <- jonckheere.terpstra(ch7$braking.distance, ch7$speed)
  expect_equal(tmp$pval.exact.stat, 32.5)
  expect_equal(tmp$pval.exact, 0.0011111111)
  expect_equal(jonckheere.terpstra(ch7$braking.distance, ch7$speed,
                                   do.exact = FALSE,
                                   do.asymp = TRUE)$pval.asymp, 0.00219353184)
})

test_that("Exercise 7.6", {
  expect_equal(jonckheere.terpstra(ch7$platelet.count, ch7$spleen.size,
                                   alternative = "greater", do.exact = FALSE,
                                   do.asymp = TRUE)$pval.asymp, 0.0000194445645)
})

test_that("Exercise 8.10", {
  expect_equal(jonckheere.terpstra(ch8$silver.content, ch8$dynasty,
                                   alternative = "less", do.exact = FALSE,
                                   do.asymp = TRUE)$pval.asymp, 0.0084307191)
})

test_that("Example 12.7", {
  expect_equal(
    jonckheere.terpstra(ch12$dose.side.effect, ch12$dose,
                        alternative = "less", do.asymp = TRUE,
                        do.exact = FALSE,
                        do.asymp.ties.adjust = FALSE)$pval.asymp.stat,
    -0.74273281)
  expect_equal(
    jonckheere.terpstra(ch12$dose.side.effect, ch12$dose,
                        alternative = "less", do.asymp = TRUE,
                        do.exact = FALSE,
                        do.asymp.ties.adjust = FALSE)$pval.asymp, 0.22882173)
})

test_that("Example 12.8", {
  expect_equal(jonckheere.terpstra(ch12$dose.side.effect, ch12$dose,
                                   alternative = "less", do.asymp = TRUE,
                                   do.exact = FALSE)$pval.asymp.stat,
               -2.1248932)
  expect_equal(jonckheere.terpstra(ch12$dose.side.effect, ch12$dose,
                                   alternative = "less", do.asymp = TRUE,
                                   do.exact = FALSE)$pval.asymp, 0.0167977625)
})

test_that("Example 12.9", {
  expect_equal(jonckheere.terpstra(ch12$platelet.count, ch12$spleen.size,
                                   alternative = "greater", do.exact = FALSE,
                                   do.asymp = TRUE)$pval.asymp.stat,
               -3.40771709)
  expect_equal(jonckheere.terpstra(ch12$platelet.count, ch12$spleen.size,
                                   alternative = "greater", do.exact = FALSE,
                                   do.asymp = TRUE)$pval.asymp, 0.00032754388)
})

test_that("Exercise 12.6", {
  expect_equal(
    jonckheere.terpstra(ch12$ethnic.group, ch12$diabetes.status,
                        do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
    0.025386376)
})
