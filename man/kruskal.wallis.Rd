% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kruskal.wallis.R
\name{kruskal.wallis}
\alias{kruskal.wallis}
\title{Perform Kruskal-Wallis test}
\usage{
kruskal.wallis(
  x,
  g,
  max.exact.cases = 15,
  nsims.mc = 10000,
  seed = NULL,
  do.asymp = FALSE,
  do.exact = TRUE,
  do.mc = FALSE
)
}
\arguments{
\item{x}{Numeric vector or factor of same length as g}

\item{g}{Factor of same length as x}

\item{max.exact.cases}{Maximum number of cases allowed for exact calculations (defaults to \code{15})}

\item{nsims.mc}{Number of Monte Carlo simulations to be performed (defaults to \code{10000})}

\item{seed}{Random number seed to be used for Monte Carlo simulations (defaults to \code{NULL})}

\item{do.asymp}{Boolean indicating whether or not to perform asymptotic calculations (defaults to \code{FALSE})}

\item{do.exact}{Boolean indicating whether or not to perform exact calculations (defaults to \code{TRUE})}

\item{do.mc}{Boolean indicating whether or not to perform Monte Carlo calculations (defaults to \code{FALSE})}
}
\value{
An ANSMtest object with the results from applying the function
}
\description{
\code{kruskal.wallis()} performs the Kruskal-Wallis test and is used in chapters 7 and 12 of "Applied Nonparametric Statistical Methods" (5th edition)
}
\examples{
# Example 7.1 from "Applied Nonparametric Statistical Methods" (5th edition)
kruskal.wallis(ch7$affordability, ch7$regions, do.exact = FALSE, do.asymp = TRUE)

# Exercise 7.16 from "Applied Nonparametric Statistical Methods" (5th edition)
kruskal.wallis(ch7$affordability, ch7$regions)

}
