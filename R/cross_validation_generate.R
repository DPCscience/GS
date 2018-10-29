#' Get the cross validation, generate five data, 0 as missing value
#'
#' @param dat
#' @param trait the trait want to cross vation
#' @param seed set the seed number, default is 123
#' @return finve files in the working directory
#' @examples
#' 1+1

cross_validation_generate <- function(dat,trait,seed=123){
  library(caret)
  w = createFolds(1:dim(dat)[1],k=5)
  dat1 = dat2 = dat3 = dat4 = dat5 = dat

  dat1[w[[1]],trait] =0
  dat2[w[[2]],trait] =0
  dat3[w[[3]],trait] =0
  dat4[w[[4]],trait] =0
  dat5[w[[5]],trait] =0
  # dat5[w[[5]],]$y1 =0
  # dat3[w[[3]],]$y1 =0
  # dat4[w[[4]],]$y1 =0
  # dat5[w[[5]],]$y1 =0

  fwrite(dat1,"dat1.csv")
  fwrite(dat2,"dat2.csv")
  fwrite(dat3,"dat3.csv")
  fwrite(dat4,"dat4.csv")
  fwrite(dat5,"dat5.csv")
}
