#' Build the H-inverse-matrix from G-inverse-matrix and pedigree-full and pedigree-genotype#'
#' @param G_mat It is the matrix which has rownames and colnames(ID)
#' @param ped_full It contains the full pedigree, it has three columns:ID,Sire,Dam
#' @param ped_geno It contains the pedigree that has the genotype, it is a part of the ped_full pedigree
#' @return The H-inverse-matrix form the formula
#' @examples
#' # 这里将iA22是A矩阵提取的A22，进行求逆


hinv_matrix <- function(M012,ped_full){
  library(sommer)
  library(MASS)

  library(nadiv)
  A <- as.matrix(makeA(prepPed(ped_full)))
  id <- row.names(A)
  iA <- ginv(A)
  row.names(iA) = colnames(iA) = id

  G <- A.mat(M012-1)
  diag(G) <- diag(G) + 0.01
  iG <- ginv(G)
  rownames(iG) <- colnames(iG) <- rownames(G)


  # library(asreml)
  # A_inv <- asreml.Ainverse(ped_full)$ginv
  # iA <- asreml.sparse2mat(A_inv)
  # A <- ginv(iA)
  # id <- attr(A_inv,"rowNames")
  # row.names(A) = colnames(A) = id
  # row.names(iA) = colnames(iA) = id

  list.geno <- row.names(G)
  list.nongeno <- setdiff(id,list.geno)
  iG1 <- iG[list.geno,list.geno]

  A22 <- A[list.geno,list.geno]
  iA22 <- ginv(A22)
  rownames(iA22) = colnames(iA22) <- row.names(A22)
  x22 <- iG1 - iA22
  iH11 <- iA[list.nongeno,list.nongeno]
  iH21 <- iA[list.geno,list.nongeno]
  iH12 <- t(iH21)
  iH22 <- x22
  Hinv <- cbind(rbind(iH11,iH21),rbind(iH12,x22))
  return(Hinv)
}
