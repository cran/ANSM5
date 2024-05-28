test_that("Example 8.11", {
  expect_equal(friedman.lsd(ch8$prey.preference, ch8$prey, ch8$larva,
                            c("Cyclops", "Anopheles"))$pval.asymp.stat,
               6.4649763)
  expect_equal(friedman.lsd(ch8$prey.preference, ch8$prey, ch8$larva,
                            c("Cyclops", "Anopheles"))$pval.asymp,
               1.9153983e-07)
})
