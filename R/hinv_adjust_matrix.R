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

hinv_adjust_matrix <- function(M012,ped_full,wts=c(1,0,1,1)){
  if (length(wts) != 4) stop("You need 4 wts (alpha, beta, tau, omega) in that order")
  alpha = wts[1]
  beta  = wts[2]
  tau   = wts[3]
  omega = wts[4]

  library(MASS)
  library(sommer)
  G <- A.mat(M012-1)

  # library(asreml)
  # A_inv <- asreml.Ainverse(ped_full)$ginv
  # iA <- asreml.sparse2mat(A_inv)
  # A <- ginv(iA)
  # id <- attr(A_inv,"rowNames")
  # row.names(A) = colnames(A) = id
  # row.names(iA) = colnames(iA) = id

  library(nadiv)
  A <- as.matrix(makeA(prepPed(ped_full)))
  id <- row.names(A)
  iA <- ginv(A)
  row.names(iA) = colnames(iA) = id

  list.geno <- row.names(M012)
  list.nongeno <- setdiff(id,list.geno)

  iG <- ginv(G)
  rownames(iG) <- rownames(G)
  colnames(iG) <- rownames(G)

  list.geno <- row.names(G)
  list.nongeno <- setdiff(id,list.geno)

  A22 <- A[list.geno,list.geno]

  meanG=mean(G)
  meandiagG=mean(diag(G))
  meanAgg=mean(A22)
  meandiagAgg=mean(diag(A22))
  cat("Means G is:",meanG,";Means G diag is:",meandiagG,";Means A22 is:",meanAgg,";Means diag A22 is:",meandiagAgg,"\n")
  beta=(meandiagAgg-meanAgg)/(meandiagG-meanG)
  alpha=meandiagAgg-meandiagG*beta
  cat("alpha and beta value:",alpha,beta,"\n")
  G=alpha+beta*G # 调整后的G

  iA22 <- ginv(A22)
  rownames(iA22) <- colnames(iA22) <- row.names(A22)
  x22 <- tau*iG - omega*iA22
  iH11 <- iA[list.nongeno,list.nongeno]
  iH21 <- iA[list.geno,list.nongeno]
  iH12 <- t(iH21)
  iH22 <- x22
  Hinv <- cbind(rbind(iH11,iH21),rbind(iH12,x22))
  return(Hinv)
}
