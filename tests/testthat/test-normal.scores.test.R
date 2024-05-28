test_that("Example 6.10", {
  tmp <- normal.scores.test(ch6$groupA, ch6$groupB)
  expect_equal(tmp$pval.exact.stat, -4.9861603)
  expect_equal(tmp$pval.exact, 0.0117034668)
  expect_equal(normal.scores.test(ch6$groupA, ch6$groupB, do.exact = FALSE,
                                  do.asymp = TRUE)$pval.asymp.stat, -4.9861603)
  expect_equal(normal.scores.test(ch6$groupA, ch6$groupB, do.exact = FALSE,
                                  do.asymp = TRUE)$pval.asymp, 0.010654119)
})

test_that("Exercise 6.15", {
  expect_equal(normal.scores.test(ch6$doseI, ch6$doseII)$pval.exact,
               0.0080213904)
})

test_that("Example 8.9", {
  groupAB.sum <- (ch8$periodI + ch8$periodII)[ch8$sequence == "AB"]
  groupBA.sum <- (ch8$periodI + ch8$periodII)[ch8$sequence == "BA"]
  tmp <- normal.scores.test(groupAB.sum, groupBA.sum)
  expect_equal(tmp$pval.exact.stat, -1.4942462)
  expect_equal(tmp$pval.exact, 0.27777778)
})
