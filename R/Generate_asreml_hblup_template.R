#' Generate the dat.as for asreml-w
#' @param dat Data.frame structure: dat has tow column, named ID and y
#' @return Return the file in default workong directory named:hblup_asreml_templete.as
#' @examples
#' library(asreml)
#' data(harvey)
#' Generate_asreml_template(harvey)
#'
#'

Generate_asreml_hblup_template <- function(dat){
  a <- capture.output(paste(cat("!WORKSPACE 15 !RENAME !out !ARGS 1// !DOPART $1 \nTitle:",substitute(dat),"\n\n"),
                            cat(as_temp(dat),"\n"),
                            cat(substitute(dat),"!skip 1\n
\n!part 1\n!nodisplay\n!gdense\ny  ~ mu ,\n!r ","giv(ID)","\nresidual units\n\n

              \nvpredict !define\nF total ID Residual\nH h2 ID total\n#"), sep = ""))
  write.table(a,"hblup_asreml_templete.as",col.names = F,row.names = F,quote = F)
}


