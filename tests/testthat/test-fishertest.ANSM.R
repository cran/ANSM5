test_that("Example 6.7", {
  expect_equal(fishertest.ANSM(ch6$males, ch6$females)$pval.exact,
               0.38454106)
})

test_that("Example 12.2", {
  expect_equal(fishertest.ANSM(ch12$infection.site, ch12$district)$pval.exact,
               0.026133508)
})

test_that("Example 12.5", {
  expect_equal(fishertest.ANSM(ch12$drugAB, ch12$side.effect.level)$pval.exact,
               0.0127213014)
})

test_that("Exercise 12.2", {
  expect_equal(fishertest.ANSM(ch12$welsh.language,
                                ch12$opportunities)$pval.exact, 0.042385254)
})

test_that("Exercise 12.3", {
  expect_equal(fishertest.ANSM(ch12$diagnosis,
                                ch12$position.played)$pval.exact, 0.077922078)
})

test_that("Exercise 12.5", {
  expect_equal(fishertest.ANSM(ch12$win.opinion, ch12$supporter)$pval.exact,
               0.06950145)
})

test_that("Exercise 12.6", {
  expect_equal(fishertest.ANSM(ch12$ethnic.group,
                                ch12$diabetes.status)$pval.exact, 0.071578042)
})

test_that("Example 13.1", {
  expect_equal(
    fishertest.ANSM(ch13$physical.activity[ch13$gender == "Boy"],
                     ch13$tv.viewing[ch13$gender == "Boy"])$pval.exact,
    0.044196753)
  expect_equal(
    fishertest.ANSM(ch13$physical.activity[ch13$gender == "Girl"],
                     ch13$tv.viewing[ch13$gender == "Girl"])$pval.exact,
    0.0059470492)
})

test_that("Section 13.3.1", {
  expect_equal(fishertest.ANSM(ch13$medicine[ch13$location == "Urban"],
                                ch13$response[ch13$location == "Urban"],
                                alternative = "greater")$pval.exact,
               0.0000068532897)
  expect_equal(fishertest.ANSM(ch13$medicine[ch13$location == "Rural"],
                                ch13$response[ch13$location == "Rural"],
                                alternative = "greater")$pval.exact,
               0.00169314846)
  expect_equal(fishertest.ANSM(ch13$medicine, ch13$response,
                                alternative = "less")$pval.exact, 0.094564075)
})

test_that("Section 13.4", {
  expect_equal(fishertest.ANSM(ch13$group[ch13$company == "A"],
                                ch13$promoted[ch13$company == "A"],
                                alternative = "less")$pval.exact, 0.076517329)
  expect_equal(fishertest.ANSM(ch13$group[ch13$company == "B"],
                                ch13$promoted[ch13$company == "B"],
                                alternative = "less")$pval.exact, 0.035089803)
  expect_equal(fishertest.ANSM(ch13$group, ch13$company)$pval.exact,
               0.030086355)
})

test_that("Exercise 13.3", {
  expect_equal(
    fishertest.ANSM(ch13$breakfast.eaten[ch13$boys.girls == "Boys"],
                     ch13$VEL[ch13$boys.girls == "Boys"])$pval.exact,
    0.6897803)
  expect_equal(
    fishertest.ANSM(ch13$breakfast.eaten[ch13$boys.girls == "Girls"],
                     ch13$VEL[ch13$boys.girls == "Girls"])$pval.exact,
    0.0300630155)
})

test_that("Exercise 13.7", {
  expect_equal(
    fishertest.ANSM(ch13$medicine[ch13$location == "Urban"],
                     ch13$response[ch13$location == "Urban"])$pval.exact,
    1.17710918e-05)
  expect_equal(
    fishertest.ANSM(ch13$medicine[ch13$location == "Rural"],
                     ch13$response[ch13$location == "Rural"])$pval.exact,
    0.0032792054)
})

test_that("Exercise 13.10", {
  expect_equal(fishertest.ANSM(ch13$laid.off, ch13$employee.ages,
                                alternative = "greater")$pval.exact, 0.008529354)
  expect_equal(fishertest.ANSM(ch13$laid.off, ch13$employee.ages)$pval.exact,
               0.016494163)
})
