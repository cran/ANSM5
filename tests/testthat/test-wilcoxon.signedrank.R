test_that("Example 3.4", {
  tmp <- wilcoxon.signedrank(ch3$heartrates1, 70, "greater")
  expect_equal(tmp$pval.exact.stat,
               "4 (sum of negative ranks), 24 (sum of positive ranks)")
  expect_equal(tmp$pval.exact, 0.0546875)
})

test_that("Example 3.5", {
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15, "less")$pval.exact.stat,
               "64 (sum of negative ranks), 14 (sum of positive ranks)")
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15,
                                   "less")$pval.exact, 0.026123047)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15)$pval.exact.stat,
               "64 (sum of negative ranks), 14 (sum of positive ranks)")
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15)$pval.exact, 0.052246094)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15, do.asymp = TRUE,
                                   cont.corr = FALSE)$pval.asymp, 0.049860204)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15, do.asymp = TRUE)$pval.asymp, 0.054613544)
})

test_that("Example 3.6", {
  expect_equal(wilcoxon.signedrank(ch3$heartrates2)$CI.exact.lower, 2.5)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2)$CI.exact.upper, 16)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2)$actualCIwidth.exact, 0.95751953)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, CI.width = 0.99)$CI.exact.lower, 0.5)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, CI.width = 0.99)$CI.exact.upper, 18)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, CI.width = 0.99)$actualCIwidth.exact, 0.99072266)
})

test_that("Example 3.7", {
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15, do.asymp = TRUE)$pval.asymp, 0.054613544)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2, 15, do.asymp = TRUE,
                                   cont.corr = FALSE)$pval.asymp, 0.049860204)
})

test_that("Example 3.8", {
  expect_equal(wilcoxon.signedrank(ch3$withties, 6, "greater")$pval.exact.stat,
               "5 (sum of negative ranks), 23 (sum of positive ranks)")
  expect_equal(wilcoxon.signedrank(ch3$withties, 6, "greater")$pval.exact, 0.0859375)
  expect_equal(wilcoxon.signedrank(ch3$tiedifrounded1, 6, "greater")$pval.exact.stat,
               "3 (sum of negative ranks), 25 (sum of positive ranks)")
  expect_equal(wilcoxon.signedrank(ch3$tiedifrounded1, 6, "greater")$pval.exact, 0.0390625)
  expect_equal(wilcoxon.signedrank(ch3$tiedifrounded2, 6, "greater")$pval.exact.stat,
               "7 (sum of negative ranks), 21 (sum of positive ranks)")
  expect_equal(wilcoxon.signedrank(ch3$tiedifrounded2, 6, "greater")$pval.exact, 0.1484375)
})

test_that("Example 3.9", {
  expect_equal(wilcoxon.signedrank(ch3$ages, 30, "greater")$pval.exact.stat,
               "16.5 (sum of negative ranks), 61.5 (sum of positive ranks)")
  expect_equal(wilcoxon.signedrank(ch3$ages, 30, "greater")$pval.exact, 0.040527344)
})

test_that("Example 3.12", {
  expect_equal(wilcoxon.signedrank(ch3$sampleA, 9)$pval.exact, 0.1015625)
  expect_equal(wilcoxon.signedrank(ch3$sampleA, 9, do.asymp = TRUE)$pval.asymp, 0.092388857)
  expect_equal(wilcoxon.signedrank(ch3$sampleB, 9)$pval.exact, 0.0148468018)
  expect_equal(wilcoxon.signedrank(ch3$sampleB, 9, do.asymp = TRUE)$pval.asymp, 0.0168522505)
})

test_that("Example 3.13", {
  expect_equal(wilcoxon.signedrank(ch3$sampleB, do.asymp = TRUE)$CI.asymp.lower, 7.15)
  expect_equal(wilcoxon.signedrank(ch3$sampleB, do.asymp = TRUE)$CI.asymp.upper, 8.45)
})

test_that("Example 3.14", {
  expect_equal(wilcoxon.signedrank(ch3$heartrates2a)$CI.exact.lower, 2.5)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2a)$CI.exact.upper, 17)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2a, CI.width = 0.99)$CI.exact.lower, 0.5)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2a, CI.width = 0.99)$CI.exact.upper, 20)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2b)$CI.exact.lower, 2.5)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2b)$CI.exact.upper, 18.5)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2b, CI.width = 0.99)$CI.exact.lower, 0.5)
  expect_equal(wilcoxon.signedrank(ch3$heartrates2b, CI.width = 0.99)$CI.exact.upper, 34.5)
})

test_that("Exercise 3.6", {
  expect_equal(wilcoxon.signedrank(ch3$children, 400)$pval.exact, 0.12939453)
  expect_equal(wilcoxon.signedrank(ch3$children, 400)$CI.exact.lower, 200.5)
  expect_equal(wilcoxon.signedrank(ch3$children, 400)$CI.exact.upper, 433)
})

test_that("Exercise 3.7", {
  expect_equal(wilcoxon.signedrank(ch3$fishlengths, 73.5)$pval.exact, 0.0093478277)
  expect_equal(wilcoxon.signedrank(ch3$fishlengths, 73.5, do.asymp = TRUE)$pval.asymp, 0.0104156196)
})

test_that("Exercise 3.8", {
  expect_equal(wilcoxon.signedrank(ch3$ages)$CI.exact.lower, 26.5)
  expect_equal(wilcoxon.signedrank(ch3$ages)$CI.exact.upper, 74)
  expect_equal(wilcoxon.signedrank(ch3$ages, CI.width = 0.99)$CI.exact.lower, 22)
  expect_equal(wilcoxon.signedrank(ch3$ages, CI.width = 0.99)$CI.exact.upper, 85)
})

test_that("Exercise 3.10", {
  expect_equal(wilcoxon.signedrank(ch3$sleeptime, 2, alternative = "greater")$pval.exact, 0.024902344)
})

test_that("Exercise 3.11", {
  expect_equal(wilcoxon.signedrank(ch3$fishlengths, do.asymp = TRUE)$CI.exact.lower, 70)
  expect_equal(wilcoxon.signedrank(ch3$fishlengths, do.asymp = TRUE)$CI.exact.upper, 73)
})

test_that("Exercise 3.12", {
  expect_equal(wilcoxon.signedrank(ch3$weightloss, 5)$pval.exact, 0.03137207)
})

test_that("Exercise 3.13", {
  expect_equal(wilcoxon.signedrank(ch3$plants, 50)$pval.exact, 0.65339544)
})

test_that("Exercise 3.14", {
  expect_equal(wilcoxon.signedrank(ch3$birthprops, CI.width = 0.99)$CI.exact.lower, 0.28250)
  expect_equal(wilcoxon.signedrank(ch3$birthprops, CI.width = 0.99)$CI.exact.upper, 0.29550)
})

test_that("Exercise 3.15", {
  expect_equal(wilcoxon.signedrank(ch3$assembly, 14.2, do.asymp = TRUE)$pval.exact, 0.0050659180)
  expect_equal(wilcoxon.signedrank(ch3$assembly, 14.2, do.asymp = TRUE)$CI.exact.lower, 12.25)
  expect_equal(wilcoxon.signedrank(ch3$assembly, 14.2, do.asymp = TRUE)$CI.exact.upper, 13.7)
})

test_that("Exercise 3.16", {
  expect_equal(wilcoxon.signedrank(ch3$weightchange, 0)$pval.exact, 0.00041832903)
  expect_equal(wilcoxon.signedrank(ch3$weightchange)$CI.exact.lower, -2.15)
  expect_equal(wilcoxon.signedrank(ch3$weightchange)$CI.exact.upper, -0.55)
})

test_that("Exercise 3.17", {
  expect_equal(wilcoxon.signedrank(ch3$sampleI, 110)$pval.exact, 0.044921875)
  expect_equal(wilcoxon.signedrank(ch3$sampleI, 110, do.asymp = TRUE)$pval.asymp, 0.046710479)
  expect_equal(wilcoxon.signedrank(ch3$sampleII, 110)$pval.exact, 0.0175762177)
  expect_equal(wilcoxon.signedrank(ch3$sampleII, 110, do.asymp = TRUE)$pval.asymp, 0.019622174)
})

test_that("Exercise 4.15", {
  expect_equal(wilcoxon.signedrank(ch4$arrow.angles, 145)$pval.exact,
               0.048828125)
})

test_that("Example 5.1", {
  expect_equal(wilcoxon.signedrank(ch5$LVF - ch5$RVF, 0)$pval.exact,
               0.0009765625)
  expect_equal(wilcoxon.signedrank(ch5$LVF - ch5$RVF)$CI.exact.lower,
               9.5)
  expect_equal(wilcoxon.signedrank(ch5$LVF - ch5$RVF)$CI.exact.upper,
               23.5)
})

test_that("Example 5.2", {
  expect_equal(wilcoxon.signedrank(ch5$arithmetic, 10)$pval.exact,
               0.248046875)
  expect_equal(wilcoxon.signedrank(ch5$arithmetic, 10,
                                   do.asymp = TRUE)$pval.asymp,
               0.22887085)
  expect_equal(wilcoxon.signedrank(ch5$arithmetic, 10,
                                   do.asymp = TRUE)$CI.asymp.lower, 5)
  expect_equal(wilcoxon.signedrank(ch5$arithmetic, 10,
                                   do.asymp = TRUE)$CI.asymp.upper, 12)
})

test_that("Example 5.3", {
  expect_equal(wilcoxon.signedrank(ch5$bp)$CI.exact.lower, 17.5)
  expect_equal(wilcoxon.signedrank(ch5$bp)$CI.exact.upper, 33.5)
  expect_equal(wilcoxon.signedrank(ch5$bp.incorrect)$CI.exact.lower, 16)
  expect_equal(wilcoxon.signedrank(ch5$bp.incorrect)$CI.exact.upper, 32.5)
})

test_that("Exercise 5.1", {
  expect_equal(wilcoxon.signedrank(ch5$bp.diff, 0)$pval.exact, 0.025390625)
})

test_that("Exercise 5.2", {
  expect_equal(wilcoxon.signedrank(ch5$LabI - ch5$LabII, 0)$pval.exact,
               0.044921875)
  expect_equal(wilcoxon.signedrank(
    ch5$LabI - ch5$LabII, 0)$CI.exact.lower, 0.05)
  expect_equal(wilcoxon.signedrank(
    ch5$LabI - ch5$LabII, 0)$CI.exact.upper, 0.85)
  expect_equal(wilcoxon.signedrank(ch5$LabI - ch5$LabII, 0,
                                   CI.width = 0.99)$CI.exact.lower, -0.25)
  expect_equal(wilcoxon.signedrank(ch5$LabI - ch5$LabII, 0,
                                   CI.width = 0.99)$CI.exact.upper, 1.05)
})

test_that("Exercise 5.5", {
  expect_equal(wilcoxon.signedrank(ch5$online - ch5$lectures,
                                   0)$pval.exact, 0.07421875)
})

test_that("Exercise 5.9", {
  expect_equal(wilcoxon.signedrank(ch5$additiveA - ch5$additiveB,
                                   0)$pval.exact, 0.12890625)
})

test_that("Exercise 5.10", {
  expect_equal(wilcoxon.signedrank(ch5$round3 - ch5$round2, 3,
                                   alternative = "greater")$pval.exact,
               0.0283203125)
})

test_that("Exercise 5.11", {
  expect_equal(wilcoxon.signedrank(ch5$pollA - ch5$pollB,
                                   0)$pval.exact, 0.27056885)
})

test_that("Exercise 5.12", {
  expect_equal(wilcoxon.signedrank(ch5$kHz0.125 - ch5$kHz0.25,
                                   0)$pval.exact, 0.03125)
})
