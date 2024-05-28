test_that("Example 4.2", {
  tmp <- suppressWarnings(kstest.ANSM(ch4$breaks, "punif", min = 0, max = 6,
                                       alternative = "two.sided"))
  expect_equal(tmp$pval.exact.stat, 0.183333333)
  expect_equal(tmp$pval.exact, 0.45829257)
  tmp <- suppressWarnings(kstest.ANSM(ch4$breaks, "punif", min = 0, max = 6,
                                      alternative = "two.sided",
                                      do.exact = FALSE, do.asymp = TRUE))
  expect_equal(tmp$pval.asymp.stat, 0.183333333)
  expect_equal(tmp$pval.asymp, 0.51214416)
})

test_that("Exercise 4.2", {
  expect_equal(suppressWarnings(
    kstest.ANSM(ch4$waiting.time, "pexp", rate = 1 / 20)$pval.exact),
    0.16521394)
})

test_that("Exercise 4.3", {
  expect_equal(kstest.ANSM(ch4$visiting.supporters, "pexp",
                            rate = 2600)$pval.exact, 0)
  expect_equal(kstest.ANSM(ch4$visiting.supporters, "punif", min = 1100,
                            max = 7700)$pval.exact, 0.000227738766)
})

test_that("Example 6.15", {
  expect_equal(kstest.ANSM(ch6$salivaF, ch6$salivaM,
                            alternative = "greater")$pval.exact.stat,
               0.52380952)
  expect_equal(kstest.ANSM(ch6$salivaF, ch6$salivaM,
                            alternative = "greater")$pval.exact, 0.046479004)
  expect_equal(kstest.ANSM(ch6$salivaF, ch6$salivaM,
                            alternative = "two.sided")$pval.exact, 0.092958008)
})

test_that("Exercise 6.13", {
  expect_equal(kstest.ANSM(ch6$time.withoutLD,
                            ch6$time.withLD)$pval.exact, 0.055944056)
  expect_equal(kstest.ANSM(ch6$time.withoutLD, ch6$time.withLD,
                            alternative = "greater")$pval.exact, 0.027972028)
})

test_that("Exercise 6.15", {
  expect_equal(kstest.ANSM(ch6$doseI, ch6$doseII)$pval.exact,
               0.033566434)
})

test_that("Exercise 6.16", {
  expect_equal(kstest.ANSM(ch6$travel, ch6$politics)$pval.exact,
               0.117214055)
})

test_that("Exercise 9.2", {
  expect_equal(kstest.ANSM(ch9$boys.toothtime, ch9$girls.toothtime)$pval.exact,
               0.994672)
})
