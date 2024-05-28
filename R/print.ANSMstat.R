#' Prints an ANSMstat object
#'
#' @description
#' `print.ANSMstat()` prints the output contained in an ANSMstat object
#'
#' @param x An ANSMstat object
#' @param ... Further arguments relevant to the default `print` function
#' @returns No return value, called to display results
#'
#' @export
print.ANSMstat <- function(x, ...) {
  #print space
  cat("\n")

  #print stat name
  if (is.null(x$varname2)){
    cat(paste0(x$title, " for ", x$varname1), "\n")
  }else if(is.null(x$varname3)){
    cat(paste0(x$title, " for ", x$varname1, " and ", x$varname2, "\n"))
  }else{
    cat(paste0(x$title, " for ", x$varname1, " and ", x$varname2,
               "\nby levels of ", x$varname3, "\n"))
  }

  #print space
  cat("\n")

  #print stat and CI
  if (!is.null(x$stat)){
    for (i in 1:length(x$stat)){
      #print stat
      if (is.character(x$stat[i])){
        cat(paste0(x$statlabel[i], ": ", x$stat[i]), "\n")
      }else{
        cat(paste0(x$statlabel[i], ": ", round(x$stat[i],5)), "\n")
      }
      #print exact CI
      if (!is.null(x$CI.exact.lower) && !is.null(x$CI.exact.upper)){
        if (is.na(x$CI.exact.lower) | is.na(x$CI.exact.upper)){
          cat(paste0("Exact ", 100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval cannot be calculated"), "\n")
        }else if (!is.null(x$actualCIwidth.exact)){
          cat(paste0("Exact ", 100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval (",
                     100 * round(x$actualCIwidth.exact, 5), "% achieved)"), "\n")
          cat(paste0("(", sprintf("%.5f", x$CI.exact.lower), ", ",
                     sprintf("%.5f", x$CI.exact.upper), ")"), "\n")
        }else{
          cat(paste0("Exact ", 100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval"), "\n")
          cat(paste0("(", sprintf("%.5f", x$CI.exact.lower), ", ",
                     sprintf("%.5f", x$CI.exact.upper), ")"), "\n")
        }
      }
      if (!is.null(x$CI.exact.note)){cat(x$CI.exact.note, "\n")}
      #print asymp CI
      if (!is.null(x$CI.asymp.lower) && !is.null(x$CI.asymp.upper)){
        if (is.na(x$CI.asymp.lower) | is.na(x$CI.asymp.upper)){
          cat(paste0("Asymptotic ", 100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval cannot be calculated"), "\n")
        }else{
          cat(paste0("Asymptotic ", 100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval"), "\n")
          cat(paste0("(", sprintf("%.5f", x$CI.asymp.lower), ", ",
                     sprintf("%.5f", x$CI.asymp.upper), ")"), "\n")
        }
      }
      if (!is.null(x$CI.asymp.note)){cat(x$CI.asymp.note, "\n")}
      #print Monte Carlo CI
      if (!is.null(x$CI.mc.lower) && !is.null(x$CI.mc.upper)){
        if (length(grep("bootstrap", x$CI.mc.note, ignore.case = TRUE)) != 0){
          CItype <- "Bootstrap"
        }else{
          CItype <- "Monte Carlo"
        }
        if (is.na(x$CI.mc.lower) | is.na(x$CI.mc.upper)){
          cat(paste0(CItype, " ", 100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval cannot be calculated"), "\n")
        }else{
          cat(paste0(CItype, " ", 100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval (", sprintf("%1$d",x$nsims.mc),
                     " simulations)"), "\n")
          cat(paste0("(", sprintf("%.5f", x$CI.mc.lower), ", ",
                     sprintf("%.5f", x$CI.mc.upper), ")"), "\n")
        }
      }
      if (!is.null(x$CI.mc.note)){cat(x$CI.mc.note, "\n")}
      #print sample-based CI that is neither exact nor asymptotic
      if (!is.null(x$CI.sample.lower) && !is.null(x$CI.sample.upper)){
        if (is.na(x$CI.sample.lower) | is.na(x$CI.sample.upper)){
          cat(paste0(100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval cannot be calculated"), "\n")
        }else{
          cat(paste0(100 * round(x$targetCIwidth, 5),
                     "% Confidence Interval (", sprintf("%1$d",x$nsims.mc),
                     " simulations)"), "\n")
          cat(paste0("(", sprintf("%.5f", x$CI.sample.lower), ", ",
                     sprintf("%.5f", x$CI.sample.upper), ")"), "\n")
        }
      }
      if (!is.null(x$CI.sample.note)){cat(x$CI.sample.note, "\n")}
    }
  }

  #print space
  if (!is.null(x$stat))
  {cat("\n")}

  #print hypotheses
  if (!is.null(x$H0)){
    if (is.character(x$H0)){
      cat(x$H0)
    }else{
      cat(paste0("Null hypothesis: theta = ", x$H0), "\n")
      if (x$alternative == "two.sided") {
        cat(paste0("Alternative hypothesis (2-sided): theta != ", x$H0), "\n")
      }else if (x$alternative == "less") {
        cat(paste0("Alternative hypothesis (1-sided): theta < ", x$H0), "\n")
      }else{
        cat(paste0("Alternative hypothesis (1-sided): theta > ", x$H0), "\n")
      }
    }
    #print space
    cat("\n")
  }

  #print pval.exact
  if (!is.null(x$pval.exact)){
    if (!is.null(x$pval.exact.stat)){
      if (is.character(x$pval.exact.stat)){
        cat(paste0("Statistic for exact test: ", x$pval.exact.stat), "\n")
      }else{
        cat(paste0("Statistic for exact test: ", round(x$pval.exact.stat,5)), "\n")
      }
    }
    if (is.na(x$pval.exact)){
      cat("Exact p-value cannot be calculated\n")
    }else if (x$pval.exact < 0.00001){
      cat("Exact p-value: < 0.00001\n")
    }else{
      cat(paste0("Exact p-value: ", sprintf("%.5f",x$pval.exact), "\n"))
    }
  }
  if (!is.null(x$pval.exact.note)){cat(x$pval.exact.note, "\n")}

  #print space
  if (!is.null(x$pval.exact) | !is.null(x$pval.exact.note))
  {cat("\n")}

  #print pval
  if (!is.null(x$pval)){
    if (!is.null(x$pval.stat)){
      if (is.character(x$pval.stat)){
        cat(paste0("Statistic for test: ", x$pval.stat), "\n")
      }else{
        cat(paste0("Statistic for test: ", round(x$pval.stat,5)), "\n")
      }
    }
    if (is.na(x$pval)){
      cat("p-value cannot be calculated\n")
    }else if (is.character(x$pval)){
      cat(x$pval, "\n")
    }else if (x$pval < 0.00001){
      cat("p-value: < 0.00001\n")
    }else{
      cat(paste0("p-value: ", sprintf("%.5f",x$pval), "\n"))
    }
  }
  if (!is.null(x$pval.note)){cat(x$pval.note, "\n")}

  #print space
  if (!is.null(x$pval) | !is.null(x$pval.note))
  {cat("\n")}

  #print pval.asymp
  if (!is.null(x$pval.asymp)){
    if (!is.null(x$pval.asymp.stat)){
      if (is.character(x$pval.asymp.stat)){
        cat(paste0("Statistic for asymptotic test: ", x$pval.asymp.stat), "\n")
      }else{
        cat(paste0("Statistic for asymptotic test: ", round(x$pval.asymp.stat,5)), "\n")
      }
    }
    if (!is.null(x$cont.corr)){
      if (x$cont.corr == TRUE){
        pval.asymp.label = "Asymptotic p-value (continuity correction used): "
      }else{
        pval.asymp.label = "Asymptotic p-value (continuity correction not used): "
      }
    }else{
      pval.asymp.label = "Asymptotic p-value: "
    }
    if (is.na(x$pval.asymp)){
      cat("Asymptotic p-value cannot be calculated\n")
    }else if (x$pval.asymp < 0.00001){
      cat(paste0(pval.asymp.label, "< 0.00001\n"))
    }else{
      cat(paste0(pval.asymp.label, sprintf("%.5f",x$pval.asymp), "\n"))
    }
  }
  if (!is.null(x$pval.asymp.note)){cat(x$pval.asymp.note, "\n")}

  #print space
  if (!is.null(x$pval.asymp) | !is.null(x$pval.asymp.note))
  {cat("\n")}

  #print pval.mc
  if (!is.null(x$pval.mc)){
    if (!is.null(x$pval.mc.stat)){
      if (is.character(x$pval.mc.stat)){
        cat(paste0("Statistic for Monte Carlo test: ", x$pval.mc.stat), "\n")
      }else{
        cat(paste0("Statistic for Monte Carlo test: ", round(x$pval.mc.stat,5)), "\n")
      }
    }
    pval.mc.label = paste0("Monte Carlo p-value (", sprintf("%1$d",x$nsims.mc),
                           " simulations): ")
    if (x$pval.mc < 0.00001){
      cat(paste0(pval.mc.label, "< 0.00001\n"))
    }else{
      if (x$nsims.mc >= 100000){
        cat(paste0(pval.mc.label, sprintf("%.5f",x$pval.mc), "\n"))
      }else if (x$nsims.mc >= 10000){
        cat(paste0(pval.mc.label, sprintf("%.4f",x$pval.mc), "\n"))
      }else{
        cat(paste0(pval.mc.label, sprintf("%.3f",x$pval.mc), "\n"))
      }
    }
  }
  if (!is.null(x$pval.mc.note)){cat(x$pval.mc.note, "\n")}

  #print space
  if (!is.null(x$pval.mc) | !is.null(x$pval.mc.note))
  {cat("\n")}

  #print stat note
  if (!is.null(x$stat.note)){
    cat(x$stat.note, "\n")
    cat("\n")
  }
}

