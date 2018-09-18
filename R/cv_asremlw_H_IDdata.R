#' Generate the asreml template as and data file using Gmatrix and dat
#'
#' @param dat contins two columns: ID and y
#' @param Gmatrix That is G-matrix, build by SNP marker, can be A.mat in sommer
#' @param seed set the seed number, default is 123
#' @return result
#' @examples
#' 1+1

cv_asremlw_H_genotypeID <- function(dat_full,dat_genotype,Gmatrix,seed=123){
  library(dplyr)
  library(data.table)
  library(caret)
  library(synbreed)
  dat <- dat_full
  dat <- unique(dat)
  dat <- na.omit(dat)
  G <- Gmatrix
  fwrite(as.data.frame(row.names(G)),"id.csv")
  ginv <- write.relationshipMatrix(G,sorting = "ASReml",type="ginv")
  fwrite(ginv,"Hinv.giv",col.names = FALSE)

  set.seed(seed)
  # dat <- dat3
  tt <- dat
  dd <- dat_genotype

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
  fwrite(tt,"dat.csv",na = "*")
  dat.csv <- tt
  Generate_asreml_template(dat.csv)
  #system("asreml dat_model.as")
}
