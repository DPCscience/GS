#' Change the Siminute SNP data to 012 format
#'
#' @param geno A format that has row.name(ID), and first column is first allet of SNP(1,2)
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' gg <- data.frame(V1 = c(1,1,2,2,1,2,1,2),V2 = c(1,1,2,2,1,2,1,2),V3 = c(2,1,2,2,1,1,2,1),V4=round(runif(8,1,2)),V5=round(runif(8,1,2)),V6=round(runif(8,1,2)))
#' row.names(gg) <- paste("ID",1:8,sep="")
#' gg
#' snp_trans_012(gg)
#'
snp_trans_012 <- function(geno) {

  Time = proc.time() # begin
  cat("Begin... \n") # begin

  X <- matrix(0,dim(geno)[1],dim(geno)[2]/2)
  for(i in 1:dim(X)[2]){
    X[,i] <- paste(geno[,c(2*i-1)],geno[,2*i],sep = "")
  }
  row.names(X) <- row.names(geno)
  dim(X)
  colnames(X) <- paste("snp",1:dim(X)[2],sep = "")
  X[which(X=="11")] <- "AA"
  X[which(X=="12")] <- "AT"
  X[which(X=="21")] <- "AT"
  X[which(X=="22")] <- "TT"
  library(snpReady)
  geno.ready <- raw.data(data = as.matrix(X), frame = "wide", base = TRUE, sweep.sample = 0.5, call.rate = 0.95, maf = 0.01, imput = FALSE)
  M <- geno.ready$M.clean

  Time = as.matrix(proc.time() - Time) #end
  cat("\n", "Completed! Time =", Time[3]/60, " minutes \n") #end

  return(M)
}
