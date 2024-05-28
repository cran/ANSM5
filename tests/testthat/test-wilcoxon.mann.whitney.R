test_that("Examples 6.1 and 6.2", {
  expect_equal(
    wilcoxon.mann.whitney(ch6$groupA, ch6$groupB)$pval.exact.stat,
    paste0("\n","76 (rank sum from ch6$groupA), ",
           "155 (rank sum from ch6$groupB)", "\n",
           "21 (Mann-Whitney U from ch6$groupA), ",
           "89 (Mann-Whitney U from ch6$groupB)"))
  expect_equal(wilcoxon.mann.whitney(ch6$groupA, ch6$groupB)$pval.exact,
               0.015871126)
})

test_that("Example 6.3", {
  expect_equal(
    wilcoxon.mann.whitney(ch6$groupA, ch6$groupB)$CI.exact.lower, -13)
  expect_equal(
    wilcoxon.mann.whitney(ch6$groupA, ch6$groupB)$CI.exact.upper, -2)
  expect_equal(wilcoxon.mann.whitney(ch6$groupA, ch6$groupB,
                                     CI.width = 0.99)$CI.exact.lower, -16)
  expect_equal(wilcoxon.mann.whitney(ch6$groupA, ch6$groupB,
                                     CI.width = 0.99)$CI.exact.upper, 1)
})

test_that("Example 6.6", {
  expect_equal(wilcoxon.mann.whitney(ch6$McGamma, ch6$McBeta,
                                     do.exact = FALSE, do.asymp = TRUE,
                                     cont.corr = FALSE)$pval.asymp, 0.41896844)
  expect_equal(wilcoxon.mann.whitney(ch6$McGamma, ch6$McBeta,
                                     do.exact = FALSE,
                                     do.asymp = TRUE)$pval.asymp, 0.425550392)
})

test_that("Example 6.8", {
  expect_equal(wilcoxon.mann.whitney(ch6$males, ch6$females, do.exact = FALSE,
                                     do.asymp = TRUE)$CI.asymp.lower,
               -1.00002198)
  expect_equal(wilcoxon.mann.whitney(ch6$males, ch6$females, do.exact = FALSE,
                                     do.asymp = TRUE)$CI.asymp.upper, 11.99999)
})

test_that("Example 6.17", {
  tmp <- wilcoxon.mann.whitney(ch6$groupA, ch6$groupB)
  expect_equal(tmp$pval.exact.stat,
               paste0("\n","76 (rank sum from ch6$groupA), ",
                      "155 (rank sum from ch6$groupB)", "\n",
                      "21 (Mann-Whitney U from ch6$groupA), ",
                      "89 (Mann-Whitney U from ch6$groupB)"))
  expect_equal(tmp$pval.exact, 0.015871126)
})

test_that("Exercise 6.1", {
  expect_equal(wilcoxon.mann.whitney(ch6$temp.H, ch6$temp.L)$pval.exact,
               0.054895105)
})

test_that("Exercise 6.3", {
  expect_equal(
    wilcoxon.mann.whitney(ch6$DMF.M, ch6$DMF.F, do.exact = FALSE,
                          do.asymp = TRUE)$pval.asymp, 0.004461907)
})

test_that("Exercise 6.5", {
  tmp <- wilcoxon.mann.whitney(ch6$cooling.time.standard,
                               ch6$cooling.time.cheap)
  expect_equal(tmp$pval.exact, 0.13879062)
  expect_equal(tmp$CI.exact.lower, -0.2)
  expect_equal(tmp$CI.exact.upper, 1.9)
})

test_that("Exercise 6.6", {
  expect_equal(
    wilcoxon.mann.whitney(ch6$wait.1979, ch6$wait.1983,
                          do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
    0.0040175478)
})

test_that("Exercise 6.7", {
  expect_equal(
    wilcoxon.mann.whitney(ch6$activity.boys, ch6$activity.girls,
                          do.exact = FALSE, do.asymp = TRUE)$pval.asymp,
    3.139189e-22)
})

test_that("Exercise 6.15", {
  tmp <- wilcoxon.mann.whitney(ch6$doseI, ch6$doseII)
  expect_equal(tmp$pval.exact, 0.0106129165)
  expect_equal(tmp$CI.exact.lower, -17.42)
  expect_equal(tmp$CI.exact.upper, -2.4)
})

test_that("Example 8.9", {
  groupAB.sum <- (ch8$periodI + ch8$periodII)[ch8$sequence == "AB"]
  groupBA.sum <- (ch8$periodI + ch8$periodII)[ch8$sequence == "BA"]
  groupAB.diff <- (ch8$periodI - ch8$periodII)[ch8$sequence == "AB"]
  groupBA.diff <- (ch8$periodI - ch8$periodII)[ch8$sequence == "BA"]
  groupBA.diff2 <- (ch8$periodII - ch8$periodI)[ch8$sequence == "BA"]
  tmp <- wilcoxon.mann.whitney(groupAB.sum, groupBA.sum)
  expect_equal(tmp$pval.exact.stat,
               paste0("\n","22 (rank sum from groupAB.sum), ",
                      "33 (rank sum from groupBA.sum)", "\n",
                      "7 (Mann-Whitney U from groupAB.sum), ",
                      "18 (Mann-Whitney U from groupBA.sum)"))
  expect_equal(tmp$pval.exact, 0.30952381)
  tmp <- wilcoxon.mann.whitney(groupAB.diff, groupBA.diff)
  expect_equal(tmp$pval.exact.stat,
               paste0("\n","16 (rank sum from groupAB.diff), ",
                      "39 (rank sum from groupBA.diff)", "\n",
                      "1 (Mann-Whitney U from groupAB.diff), ",
                      "24 (Mann-Whitney U from groupBA.diff)"))
  expect_equal(tmp$pval.exact, 0.015873016)
  tmp <- wilcoxon.mann.whitney(groupAB.diff, groupBA.diff2)
  expect_equal(tmp$pval.exact.stat,
               paste0("\n","33 (rank sum from groupAB.diff), ",
                      "22 (rank sum from groupBA.diff2)", "\n",
                      "18 (Mann-Whitney U from groupAB.diff), ",
                      "7 (Mann-Whitney U from groupBA.diff2)"))
  expect_equal(tmp$pval.exact, 0.30952381)
})

test_that("Following Example 9.3", {
  expect_equal(
    wilcoxon.mann.whitney(ch9$sampleI.survtime, ch9$sampleII.survtime.2,
                          alternative = "less")$pval.exact, 0.21756022)
})

test_that("Exercise 9.2", {
  expect_equal(
    wilcoxon.mann.whitney(ch9$boys.toothtime, ch9$girls.toothtime)$pval.exact,
    0.85181485)
})

test_that("Example 12.5", {
  side.effect.level.A <- ch12$side.effect.level[ch12$drugAB == "Drug A"]
  side.effect.level.B <- ch12$side.effect.level[ch12$drugAB == "Drug B"]
  tmp <- wilcoxon.mann.whitney(side.effect.level.A, side.effect.level.B,
                                   do.CI = FALSE, do.exact = FALSE,
                                   do.asymp = TRUE)
  expect_equal(tmp$pval.asymp.stat,
               paste0("\n","2624 (rank sum from side.effect.level.A), ",
                      "2326 (rank sum from side.effect.level.B)", "\n",
                      "1589 (Mann-Whitney U from side.effect.level.A), ",
                      "841 (Mann-Whitney U from side.effect.level.B)"))
  expect_equal(tmp$pval.asymp, 0.00184697407)
})

test_that("Exercise 12.4", {
  feedback.satisfaction.Representative <-
    ch12$feedback.satisfaction[ch12$PPI.person.2 == "Representative"]
  feedback.satisfaction.Researcher <-
    ch12$feedback.satisfaction[ch12$PPI.person.2 == "Researcher"]
  expect_equal(
    wilcoxon.mann.whitney(feedback.satisfaction.Representative,
                          feedback.satisfaction.Researcher, do.exact = FALSE,
                          do.asymp = TRUE)$pval.asymp, 0.99598853)
})
