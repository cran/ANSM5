% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fishertest.ANSM.R
\name{fishertest.ANSM}
\alias{fishertest.ANSM}
\title{Perform Fisher exact test}
\usage{
fishertest.ANSM(
  x,
  y,
  H0 = NULL,
  alternative = c("two.sided", "less", "greater"),
  max.exact.cases = 10000,
  do.exact = TRUE
)
}
\arguments{
\item{x}{Numeric vector or factor}

\item{y}{Numeric vector or factor}

\item{H0}{Null hypothesis value (defaults to \code{NULL})}

\item{alternative}{Type of alternative hypothesis (defaults to \code{two.sided})}

\item{max.exact.cases}{Maximum number of cases allowed for exact calculations (defaults to \code{10000})}

\item{do.exact}{Boolean indicating whether or not to perform exact calculations (defaults to \code{TRUE})}
}
\value{
An ANSMtest object with the results from applying the function
}
\description{
\code{fishertest.ANSM()} is a wrapper for fisher.test() from the \code{stats} package - performs the Fisher exact test and is used in chapters 6, 12 and 13 of "Applied Nonparametric Statistical Methods" (5th edition)
}
\examples{
# Example 6.7 from "Applied Nonparametric Statistical Methods" (5th edition)
fishertest.ANSM(ch6$males, ch6$females)

# Exercise 13.10 from "Applied Nonparametric Statistical Methods" (5th edition)
fishertest.ANSM(ch13$laid.off, ch13$employee.ages)

}
