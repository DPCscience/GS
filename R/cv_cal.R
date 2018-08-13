#' Build the H-matrix from G-inverse-matrix and pedigree-full and pedigree-genotype, It calculate alpha and bera based the function of the matrix G_22 and A_22
#' @param M012 It is the matrix which has rownames(ID), the first column is SNP infor, the SNP is 012 format.
#' @param ped_full It contains the full pedigree, it has three columns:ID,Sire,Dam
#' @return The H-matrix form the formula
#' @examples
#' library(MASS)
#' animal  <- 13:26

cv_cal <- function(dat.csv,h2=1){
  dat <- dat.csv
  library(dplyr)
  library(data.table)
  library(caret)
  library(synbreed)
  dat <- unique(dat)

  names(dat) <- c("ID","y")
  dat$ID <- as.factor(dat$ID)
  dd <- dat
  set.seed(123)
  w <- createFolds(1:length(dd$ID),k = 5)
  id <- as.character(dd$ID)
  id1 <- as.character(dd$ID[w[[1]]])
  id2 <- as.character(dd$ID[w[[2]]])
  id3 <- as.character(dd$ID[w[[3]]])
  id4 <- as.character(dd$ID[w[[4]]])
  id5 <- as.character(dd$ID[w[[5]]])

  names(dat) <- c("ID","y")
  options(warn = -1)
  pv1 <- fread("dat_model1/dat_model.sln");pv1=pv1[-1,-1];names(pv1)=c("ID","BLUP","SE")
  pv2 <- fread("dat_model2/dat_model.sln");pv2=pv2[-1,-1];names(pv2)=c("ID","BLUP","SE")
  pv3 <- fread("dat_model3/dat_model.sln");pv3=pv3[-1,-1];names(pv3)=c("ID","BLUP","SE")
  pv4 <- fread("dat_model4/dat_model.sln");pv4=pv4[-1,-1];names(pv4)=c("ID","BLUP","SE")
  pv5 <- fread("dat_model5/dat_model.sln");pv5=pv5[-1,-1];names(pv5)=c("ID","BLUP","SE")
  pv <- fread("dat_model6/dat_model.sln");pv=pv[-1,-1];names(pv)=c("ID","BLUP","SE")
  options(warn = 1)

  py1 <- pv1[pv1$ID %in% id1,][,1:2];py1$ID <- as.factor(py1$ID)
  py2 <- pv2[pv2$ID %in% id2,][,1:2];py2$ID <- as.factor(py2$ID)
  py3 <- pv3[pv3$ID %in% id3,][,1:2];py3$ID <- as.factor(py3$ID)
  py4 <- pv4[pv4$ID %in% id4,][,1:2];py4$ID <- as.factor(py4$ID)
  py5 <- pv5[pv5$ID %in% id5,][,1:2];py5$ID <- as.factor(py5$ID)
  py <- pv[pv5$ID %in% id,][,1:2];py$ID <- as.factor(py$ID)

  library(dplyr)
  p_list <- list(py1,py2,py3,py4,py5)
  r_pearson <- NULL
  # r_spearman <- NULL
  r_unbiased <- NULL
  # r_MSD <- NULL
  # dat[dat$y==0,]$y <- NA
  for(i in 1:5){
    # i=1
    x <- merge(p_list[[i]],dat,"ID")
    names(x) <- c("ID","pv","y")
    x <- na.omit(x)
    r_pearson[i] <- cor(x$pv,x$y,method = "pearson")/sqrt(h2)
    r_unbiased[i] <- lm(x$y ~ x$pv)$coefficients[2]
  }
  se <- function(x){
    sd(x)/sqrt(length(x))
  }
  r_pearson
  n1 <- mean(r_pearson);n1e <- se(r_pearson)
  n3 <- mean(r_unbiased);n3e <- se(r_unbiased)
  ax <- data.frame(type=c("r_pearson","r_unbiased"),
                   value = c(n1,n3),
                   se = c(n1e,n3e))
  ax
  return(ax)
}

