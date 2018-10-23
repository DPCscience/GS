#' a Build the H-inverse-matrix from G-inverse-matrix and pedigree-full and pedigree-genotype#'
#' @param G_mat It is the matrix which has rownames and colnames(ID)
#' @param ped_full It contains the full pedigree, it has three columns:ID,Sire,Dam
#' @param w a,b,tau,omega, tau[a*G + b*A22]^{-1} - omega*A22^{-1}
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
#'                         library(asreml)
#'                         ainv = asreml.Ainverse(ped)$ginv
#'                         id_full = attr(ainv,"rowNames")
#'                         id_geno = row.names(M1)
#'                         round(hinv_ainv_adjust_matrix(M1,ainv,id_full,id_geno),2)
#'
#'

hinv_ainv_adjust_matrix <- function(M012,ainv,id_full,id_geno,wts=c(0.95,0.05,1,1)){
  if (length(wts) != 4) stop("You need 4 wts (a, b, tau, omega) in that order")
  a = wts[1]
  b  = wts[2]
  tau   = wts[3]
  omega = wts[4]

  library(MASS)
  library(sommer)
  library(GS)

  id_full = as.character(id_full)
  id_geno = as.character(id_geno)

  iA <- sparse2mat(ainv)
  rownames(iA) = colnames(iA) =  id_full

  G <- A.mat(M012-1)

  id1 = setdiff(id_full,id_geno)
  id2 = id_geno

  iA22 =  iA[id2,id2] - iA[id2,id1] %*% solve(iA[id1,id1]) %*% iA[id1,id2]
  A22 <- solve(iA22)
  diagG <- diag(G)
  diagA22 <- diag(A22)
  offdiagA22 <- c(A22[upper.tri(A22)],A22[lower.tri(A22)])
  diagG <- diag(G)
  offdiagG <- c(G[upper.tri(G)],G[lower.tri(G)])
  meanoffdiagG=mean(offdiagG)
  meandiagG=mean(diagG)
  meanoffdiagA22=mean(offdiagA22)
  meandiagA22=mean(diagA22)

  beta=(meandiagA22 - meanoffdiagA22)/(meandiagG - meanoffdiagG)
  alpha=meandiagA22-meandiagG*beta
  cat("Adjust G, and the value of alpha and beta is:",alpha,beta,"\n\n\n")
  G=alpha+beta*G #
  G = a*G + b*A22
  iG <- solve(G)
  rownames(iG) = colnames(iG) = id2

  x22 <- tau*iG - omega*iA22

  iH11 <- iA[id1,id1]
  iH21 <- iA[id2,id1]
  iH12 <- t(iH21)
  iH22 <-  iA[id2,id2] + x22
  Hinv <- cbind(rbind(iH11,iH21),rbind(iH12,iH22))
  return(Hinv)
}
