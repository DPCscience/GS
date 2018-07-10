#' Build the H-inverse-matrix from G-inverse-matrix and pedigree-full and pedigree-genotype#'
#' @param G_mat It is the matrix which has rownames and colnames(ID)
#' @param ped_full It contains the full pedigree, it has three columns:ID,Sire,Dam
#' @param ped_geno It contains the pedigree that has the genotype, it is a part of the ped_full pedigree
#' @return The H-inverse-matrix form the formula
#' @examples
#' 1+1
#'

hinv_matrix <- function(G_mat,ped_full,ped_geno){
  library(asreml)
  library(MASS)
  diag(G_mat) <- diag(G_mat) + 0.01
  G_inv <- ginv(G_mat)
  names(ped_geno) <- c("Progeny","Sire","Dam")
  A1_inv <- asreml.Ainverse(ped_full)$ginv
  A1_mat <- ginv(asreml.sparse2mat(A1_inv))
  dim(A1_mat)
  A11 <- ginv(A1_mat)
  row.names(A1_mat) = colnames(A1_mat) = attr(A1_inv,"rowNames")
  id_geno <- ped_geno$Progeny
  A22_inv <- ginv(A1_mat[rownames(A1_mat) %in% id_geno,colnames(A1_mat) %in% id_geno])
  x22 <- G_inv - A22_inv
  m <- dim(A11)[1]
  n <- dim(x22)[1]
  Hinv <- A11 + cbind(rbind(matrix(0,(m-n),(m-n)),matrix(0,n,(m-n))),rbind(matrix(0,(m-n),n),x22))
  return(Hinv)
}
