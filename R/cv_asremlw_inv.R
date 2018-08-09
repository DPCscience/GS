#' Get the cv-result that run asremlw in the working directory
#'
#' @param dat contins two columns: ID and y, that the ginv and the id.csv should
#' @param Hinv That is Hinv-matrix, is inverse of H-matrix
#' @param seed set the seed number, default is 123
#' @return result
#' @examples

cv_asremlw_inv <- function(dat,Hinv,seed=123){
  library(dplyr)
  library(data.table)
  library(caret)
  library(synbreed)
  dat <- unique(dat)
  ginv <- Hinv
  fwrite(as.data.frame(rownames(Hinv)),"id.csv")
  fwrite(as.data.table(ginv),"Hinv.giv",col.names = FALSE)

  set.seed(seed)
  # dat <- dat3
  names(dat) <- c("ID","y")
  tt <- dat
  dd <- dat

  w <- createFolds(1:length(dd$ID),k = 5)
  id1 <- as.character(dd$ID[w[[1]]])
  id2 <- as.character(dd$ID[w[[2]]])
  id3 <- as.character(dd$ID[w[[3]]])
  id4 <- as.character(dd$ID[w[[4]]])
  id5 <- as.character(dd$ID[w[[5]]])
  tt$y1 = tt$y2 = tt$y3 = tt$y4 = tt$y5 = tt$y
  tt <- as.data.frame(tt)
  row.names(tt) <- tt$ID
  tt[id1,]$y1 <- NA
  tt[id2,]$y2 <- NA
  tt[id3,]$y3 <- NA
  tt[id4,]$y4 <- NA
  tt[id5,]$y5 <- NA
  fwrite(tt,"dat.csv",na = "-999")
  dat.csv <- tt
  Generate_asreml_template(dat.csv)
}
