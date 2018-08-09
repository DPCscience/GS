#' Generate the dat.as for asreml-w
#' @param dat Data.frame structure
#' @return Return the file in default workong directory named:Model-dave.as
#' @examples
#' library(asreml)
#' data(harvey)
#' Generate_asreml_template(harvey)
#'
#'

Generate_asreml_template <- function(dat){
  a <- capture.output(paste(cat("!WORKSPACE 16 !RENAME !out !ARGS 1 2 3 4 5 6 // !DOPART $1 \nTitle:",substitute(dat),"\n\n"),
          cat(as_temp(dat),"\n"),
          cat(substitute(dat),"!skip 1\n
\n!part 1\n!nodisplay\ny1  ~ mu ,\n!r ","giv(ID)","\nresidual units\nvpredict !define\nF total ID Residual\nH h2 ID total\n
\n!part 2\n!nodisplay\ny2  ~ mu ,\n!r ","giv(ID)","\nresidual units\nvpredict !define\nF total ID Residual\nH h2 ID total\n
\n!part 3\n!nodisplay\ny3  ~ mu ,\n!r ","giv(ID)","\nresidual units\nvpredict !define\nF total ID Residual\nH h2 ID total\n
\n!part 4\n!nodisplay\ny4  ~ mu ,\n!r ","giv(ID)","\nresidual units\nvpredict !define\nF total ID Residual\nH h2 ID total\n
\n!part 5\n!nodisplay\ny5  ~ mu ,\n!r ","giv(ID)","\nresidual units\nvpredict !define\nF total ID Residual\nH h2 ID total\n
\n!part 6\n!nodisplay\ny  ~ mu ,\n!r ","giv(ID)","\nresidual units
              \nvpredict !define\nF total ID Residual\nH h2 ID total\n#"), sep = ""))
  write.table(a,"dat_model.as",col.names = F,row.names = F,quote = F)
}


