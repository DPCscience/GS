#' Build the H-inverse-matrix from G-inverse-matrix and pedigree-full and pedigree-genotype#'
#' @param G_mat It is the matrix which has rownames and colnames(ID)
#' @param ped_full It contains the full pedigree, it has three columns:ID,Sire,Dam
#' @param ped_geno It contains the pedigree that has the genotype, it is a part of the ped_full pedigree
#' @return The H-inverse-matrix form the formula
#' @examples
#' animal  <- 13:26
#' data.11.1 <- data.frame(animal,
#'                         sire  = c(0,0,13,15,15,14,14,14,1,14,14,14,14,14),
#'                         dam   = c(0,0,4,2,5,6,9,9,3,8,11,10,7,12),
#'                         mean  = rep(1,length(animal)),
#'                         EDC   = c(558,722,300,73,52,87,64,103,13,125,93,66,75,33),
#'                         fat_DYD = c(9.0,13.4,12.7,15.4,5.9,7.7,10.2,4.8,7.6,8.8,9.8,9.2,11.5,13.3),
#'                         SNP1  = c(2,1,1,0,0,1,0,0,2,0,0,1,0,1),
#'                         SNP2  = c(0,0,1,0,1,1,0,1,0,0,1,0,0,0),SNP3  = c(1,0,2,2,1,0,1,1,0,0,1,0,0,1),
#'                         SNP4  = c(1,0,1,1,2,1,1,0,0,1,0,0,1,1),
#'                         SNP5  = c(0,0,1,0,0,0,0,0,0,1,0,1,1,0),
#'                         SNP6  = c(0,2,0,1,0,2,2,1,1,2,1,1,2,2),
#'                         SNP7  = c(0,0,0,0,0,0,0,0,2,0,0,0,0,0),
#'                         SNP8  = c(2,2,2,2,2,2,2,2,2,2,2,2,2,1),
#'                         SNP9  = c(1,1,1,2,1,2,2,2,1,0,2,0,1,0),
#'                         SNP10 = c(2,0,2,1,2,1,0,0,2,0,1,0,0,0))
#'                         rm(list="animal")
#'                         animal <- 1:26
#'                         sire   <- c(rep(0,12), data.11.1$sire)
#'                         dam    <- c(rep(0,12), data.11.1$dam)
#'                         ped    <- data.frame(animal, sire, dam)
#'                         rm(list=c("animal","dam","sire"))
#'                         M <- data.11.1[6:14, c(1, 7:16)]
#'                         rownames(M) <- M[, 1]
#'                         M1 <- as.matrix(M[, -1])
#'                         round(hinv_matrix(M1,ped),2)
#'
#'

hinv_matrix <- function(M012,ped_full,wts=c(0.95,0.05,1,1)){
  if (length(wts) != 4) stop("You need 4 wts (alpha, beta, tau, omega) in that order")
  alpha = wts[1]
  beta  = wts[2]
  tau   = wts[3]
  omega = wts[4]
  library(asreml)
  library(MASS)
  library(sommer)
  names(ped_full) <- c("ID","Sire","Dam")

  # alpha=0.95;beta=0.05;tau=1;omega=1

  G <- A.mat(M012-1)
  iG <- ginv(G)
  rownames(iG) <- rownames(G)
  colnames(iG) <- rownames(G)

  A_inv <- asreml.Ainverse(ped_full)$ginv
  iA <- asreml.sparse2mat(A_inv)
  A <- ginv(iA)


  row.names(A) = colnames(A) = attr(A1_inv,"rowNames")
  row.names(iA) = colnames(iA) = attr(A1_inv,"rowNames")

  list.geno <- row.names(M012)
  list.nongeno <- setdiff(ped_full$ID,list.geno)

  A_11 <- A[rownames(A) %in% list.nongeno, colnames(A) %in% list.nongeno]
  A_22 <- A[rownames(A) %in% list.geno,    colnames(A) %in% list.geno]
  A_12 <- A[rownames(A) %in% list.nongeno, colnames(A) %in% list.geno]
  A_21 <- t(A_12)

  iA_11 <- iA[rownames(iA) %in% list.nongeno, colnames(iA) %in% list.nongeno]
  iA_12 <- iA[rownames(iA) %in% list.nongeno,    colnames(iA) %in% list.geno]
  iA_22 <- iA[rownames(iA) %in% list.geno,    colnames(iA) %in% list.geno]
  iA_21 <- t(iA_12)


  Gstar <- (alpha*G) + (beta*A_22[order(rownames(G)), order(colnames(G))])

  iGiA22 <- (tau*solve(Gstar)) -  (omega*iA_22[order(rownames(Gstar)), order(colnames(Gstar))])
  rownames(iGiA22) <- rownames(Gstar)
  colnames(iGiA22) <- colnames(Gstar)

  iH_11 <- iA_11
  iH_12 <- iA_12
  iH_21 <- t(iH_12)
  iH_22 <- iA_22[order(rownames(iGiA22)), order(colnames(iGiA22))] + iGiA22[order(rownames(iGiA22)), order(colnames(iGiA22))]

  topH <- cbind(iH_11, iH_12)
  botH <- cbind(iH_21, iH_22)
  iH   <- rbind(topH, botH)
  return(iH)
}
