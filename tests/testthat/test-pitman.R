test_that("Example 3.11", {
  tmp <- pitman(ch3$heartrates1, 70, "greater")
  expect_equal(tmp$pval.exact.stat, 95)
  expect_equal(tmp$pval.exact, 0.0546875)
  expect_equal(pitman(ch3$heartrates1, 70, "greater", do.exact = FALSE,
                      do.asymp = TRUE)$pval.asymp.stat, 1.63691848888252)
  expect_equal(pitman(ch3$heartrates1, 70, "greater", do.exact = FALSE,
                      do.asymp = TRUE)$pval.asymp, 0.050823751)
})

test_that("Example 3.12", {
  expect_equal(pitman(ch3$sampleA, 9)$pval.exact, 0.3515625)
  expect_equal(pitman(ch3$sampleA, 9, do.exact = FALSE,
                      do.asymp = TRUE)$pval.asymp, 0.3051186)
  expect_equal(pitman(ch3$sampleB, 9)$pval.exact, 0.9159565)
  expect_equal(pitman(ch3$sampleB, 9, do.exact = FALSE,
                      do.asymp = TRUE)$pval.asymp, 0.81653892)
  expect_equal(pitman(ch3$sampleA2, 9)$pval.exact, 0.123046875)
  expect_equal(pitman(ch3$sampleA3, 9)$pval.exact, 0.035497665)
})

test_that("Exercise 3.17", {
  expect_equal(pitman(ch3$sampleI, 110)$pval.exact, 0.05859375)
  expect_equal(pitman(ch3$sampleI, 110, do.exact = FALSE,
                      do.asymp = TRUE)$pval.asymp, 0.064739074)
  expect_equal(pitman(ch3$sampleII, 110)$pval.exact, 0.0220565796)
  expect_equal(pitman(ch3$sampleII, 110, do.exact = FALSE,
                      do.asymp = TRUE)$pval.asymp, 0.0241634455)
})
