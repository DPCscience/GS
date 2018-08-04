#' Get the cv-result that run asremlw in the working directory
#'
#' @param dat contins two columns: ID and y, that the ginv and the id.csv should
#' @param Gmatrix That is G-matrix, build by SNP marker, can be A.mat in sommer
#' @param seed set the seed number, default is 123
#' @return result
#' @examples

cv_asremlw_H <- function(dat,Gmatrix,seed=123){
  library(dplyr)
  library(data.table)
  library(caret)
  library(synbreed)
  dat <- unique(dat)
  G <- Gmatrix
  fwrite(as.data.frame(row.names(G)),"id.csv")
  ginv <- write.relationshipMatrix(G,sorting = "ASReml",type="ginv")
  fwrite(ginv,"Hinv.giv",col.names = FALSE)

  set.seed(seed)
  # dat <- dat3
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
  fwrite(tt,"dat.csv",na = "*")
  dat.csv <- tt
  Generate_asreml_template(dat.csv)

  system("asreml dat_model.as")
  # pv <- fread("w/tt6/tt.pvs",skip = 13)
  pv1 <- fread("dat_model1/dat_model.pvs",skip = 13)
  pv2 <- fread("dat_model2/dat_model.pvs",skip = 13)
  pv3 <- fread("dat_model3/dat_model.pvs",skip = 13)
  pv4 <- fread("dat_model4/dat_model.pvs",skip = 13)
  pv5 <- fread("dat_model5/dat_model.pvs",skip = 13)
  pv <- fread("dat_model6/dat_model.pvs",skip = 13)

  py1 <- pv1[pv1$ID %in% id1,][,1:2];py1$ID <- as.factor(py1$ID)
  py2 <- pv2[pv2$ID %in% id2,][,1:2];py2$ID <- as.factor(py2$ID)
  py3 <- pv3[pv3$ID %in% id3,][,1:2];py3$ID <- as.factor(py3$ID)
  py4 <- pv4[pv4$ID %in% id4,][,1:2];py4$ID <- as.factor(py4$ID)
  py5 <- pv5[pv5$ID %in% id5,][,1:2];py5$ID <- as.factor(py5$ID)

  library(dplyr)
  p_list <- list(py1,py2,py3,py4,py5)
  r_pearson <- NULL
  r_spearman <- NULL
  r_unbiased <- NULL
  r_MSD <- NULL
  # dat[dat$y==0,]$y <- NA
  for(i in 1:5){
    # dat_omit <- na.omit(dat)
    # x <- left_join(p_list[[i]],dat,"ID")
    x <- merge(p_list[[i]],dat,"ID")
    names(x) <- c("ID","pv","y")
    x <- na.omit(x)
    r_pearson[i] <- cor(x$pv,x$y,method = "pearson")
    r_spearman[i] <- cor(x$pv,x$y,method = "spearman")
    r_unbiased[i] <- summary(lm(x$y ~ x$pv -1))$coefficients[1]
    r_MSD[i] <- mean((x$y - x$pv)^2)
  }
  se <- function(x){
    sd(x)/sqrt(length(x))
  }
  r_pearson
  n1 <- mean(r_pearson);n1e <- se(r_pearson)
  n2<- mean(r_spearman);n2e <- se(r_spearman)
  n3 <- mean(r_unbiased);n3e <- se(r_unbiased)
  n4 <- mean(r_MSD);n4e <- se(r_MSD)
  ax <- data.frame(type=c("r_pearson","r_spearman","r_unbiased","r_MSD"),
                   value = c(n1,n2,n3,n4),
                   se = c(n1e,n2e,n3e,n4e))
  ax
  return(ax)
}
