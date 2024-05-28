
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ANSM5

<!-- badges: start -->
<!-- badges: end -->

The ANSM5 package provides a number of functions and data to accompany
the book “Applied Nonparametric Statistical Methods”, 5th edition
(Smeeton, Spencer & Sprent, 2024, Taylor Francis), the revisions from
the 4th edition including a move from describing the output from a
miscellany of statistical software packages to using R. While the output
from many of the functions can also be obtained using a range of
pre-existing R functions, this package provides functions in a unified
setting and give output using both p-values and confidence intervals,
exemplifying the book’s approach of treating p-values as a guide to
statistical importance and not an end product in their own right.

## Installation

You can install the development version of ANSM5 like so:

``` r
# install.packages("devtools")
devtools::install_github("neilhspencer/ANSM5")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ANSM5)

#From Example 3.1 in "Applied Nonparametric Statistical Methods" (5th edition)
sgn.test(ch3$sampleI, 110)
#> 
#> Sign test for ch3$sampleI 
#> 
#> Null hypothesis: theta = 110 
#> Alternative hypothesis (2-sided): theta != 110 
#> 
#> Exact p-value: 0.10937
#> Exact 95% Confidence Interval (97.852% achieved) 
#> (11.00000, 151.00000)

#From Example 6.1 in "Applied Nonparametric Statistical Methods" (5th edition)
wilcoxon.mann.whitney(ch6$groupA, ch6$groupB)
#> 
#> Wilcoxon-Mann-Whitney test for ch6$groupA and ch6$groupB
#> 
#> H0: samples are from the same population
#> H1: samples differ in location
#> 
#> Statistic for exact test: 
#> 76 (rank sum from ch6$groupA), 155 (rank sum from ch6$groupB)
#> 21 (Mann-Whitney U from ch6$groupA), 89 (Mann-Whitney U from ch6$groupB) 
#> Exact p-value: 0.01587
#> Exact 95% Confidence Interval 
#> (-13.00000, -2.00000)

#From Example 7.6 in "Applied Nonparametric Statistical Methods" (5th edition)
friedman(ch7$pulse, ch7$time.period, ch7$student)
#> 
#> Friedman test for ch7$pulse and ch7$time.period (as groups) with ch7$student (as blocks)
#> 
#> H0: distributions of ch7$pulse are identical
#> H1: distributions differ
#> 
#> Statistic for exact test: 10.57143 
#> Exact p-value: 0.00272
```
