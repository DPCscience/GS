#' Generate the data structure for asreml-w
#'
#' @param dat Data fram that used to be analysised.
#' @return The structure of the data.
#' @examples
#' library(asreml)
#' data(harvey)
#' as_temp(harvey)

as_temp <- function(dat){
  dat$ID <- as.factor(dat$ID)
  for(i in 1:dim(dat)[2]){
    # i =1
    if(is.factor(dat[,i])){cat(paste(names(dat)[i],"!A",nlevels(dat[,i]),"!L id.csv !skip 1","\n"))}
    else{
      cat(paste(names(dat)[i],"!M -999 \n"))
    }
  }
  cat("\nHinv.giv")
}
