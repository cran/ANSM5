test_that("Example 8.10", {
  expect_equal(kruskal.wallis.lsd(ch8$sentences, ch8$authors,
                                  c("Vulliamy", "Queen"))$pval.asymp.stat,
               3.7294898)
  expect_equal(kruskal.wallis.lsd(ch8$sentences, ch8$authors,
                                  c("Vulliamy", "Queen"))$pval.asymp,
               0.00100687047)
})

test_that("Exercise 8.8", {
  expect_equal(kruskal.wallis.lsd(ch8$seizure.score, ch8$hospital,
                                  c("HospitalA", "HospitalC"))$pval.asymp.stat,
               2.4743409)
  expect_equal(kruskal.wallis.lsd(ch8$seizure.score, ch8$hospital,
                                  c("HospitalA", "HospitalC"))$pval.asymp,
               0.0104035113)
})
