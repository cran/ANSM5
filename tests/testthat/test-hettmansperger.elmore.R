test_that("Example 8.6", {
  expect_equal(hettmansperger.elmore(ch8$plant.weight.2, ch8$growth.hormone,
                                     ch8$undersoil.heating)$pval.asymp.stat,
               16.8033333)
  expect_equal(hettmansperger.elmore(ch8$plant.weight.2, ch8$growth.hormone,
               ch8$undersoil.heating)$pval.asymp, 4.1460411e-05)
})

test_that("Example 8.7", {
  expect_equal(hettmansperger.elmore(ch8$plant.weight.4, ch8$growth.hormone,
                                     ch8$undersoil.heating)$pval.asymp.stat,
               11.2133333)
  expect_equal(hettmansperger.elmore(ch8$plant.weight.4, ch8$growth.hormone,
                                     ch8$undersoil.heating)$pval.asymp,
               0.00081211713)
  expect_equal(hettmansperger.elmore(ch8$plant.weight.4, ch8$growth.hormone,
                                     ch8$undersoil.heating,
                                     median.polish = TRUE)$pval.asymp.stat,
               9.72)
  expect_equal(hettmansperger.elmore(ch8$plant.weight.4, ch8$growth.hormone,
                                     ch8$undersoil.heating,
                                     median.polish = TRUE)$pval.asymp,
               0.00182273517)
})

test_that("Exercise 8.3", {
  expect_equal(hettmansperger.elmore(ch8$game.time, ch8$experience,
                                     ch8$game)$pval.asymp.stat, 7.0520487)
  expect_equal(hettmansperger.elmore(ch8$game.time, ch8$experience,
                                     ch8$game)$pval.asymp, 0.0294216535)
})
