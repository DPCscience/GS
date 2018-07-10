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
  a <- capture.output(paste(cat("WORKSPACE 4G !RENAME !ARGS // !DOPART $1 \nTitle:",substitute(dat),"\n\n"),
          cat(as_temp(dat),"\n"),
          cat(substitute(dat),"!skip 1 \n \ny  ~ mu ,\n!r ","ID","\nresidual units\n#"), sep = ""))
  write.table(a,"Model-dave.as",col.names = F,row.names = F,quote = F)
}


