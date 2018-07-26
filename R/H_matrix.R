#' Build the H-inverse-matrix from G-inverse-matrix and pedigree-full and pedigree-genotype#'
#' @param M012 It is the matrix which has rownames(ID), the first column is SNP infor, the SNP is 012 format.
#' @param ped_full It contains the full pedigree, it has three columns:ID,Sire,Dam
#' @return The H-matrix form the formula
#' @examples
#' library(MASS)
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
#'                         round(ginv(H_matrix(M1,ped)),2)

H_matrix <- function(M012,ped_full){
  library(MASS)
  library(sommer)

  library(nadiv)
  A <- as.matrix(makeA(prepPed(ped_full)))
  id <- row.names(A)
  iA <- ginv(A)
  row.names(iA) = colnames(iA) = id

  G <- A.mat(M012-1)
  diag(G) <- diag(G) + 0.01
  iG <- ginv(G)
  rownames(iG) <- colnames(iG) <- rownames(G)
  genotyped=colnames(G)

  # library(asreml)
  # A_inv <- asreml.Ainverse(ped_full)$ginv
  # iA <- asreml.sparse2mat(A_inv)
  # A <- ginv(iA)
  # id <- attr(A_inv,"rowNames")
  # row.names(A) = colnames(A) = id
  # row.names(iA) = colnames(iA) = id

  genotyped=colnames(G)
  inpedigree=colnames(A)
  nongenotyped=setdiff(inpedigree,genotyped)

  cat("In pedigree nongenotyped length: ",length(nongenotyped),"\n")

  genotypednotinpedigree=setdiff(genotyped,inpedigree)
  cat("genotyped not in pedigree",length(genotypednotinpedigree),"\n")

  genotypedinpedigree=intersect(genotyped,inpedigree)
  cat("genotyped in pedigree",length(genotypedinpedigree),"\n")
  G=G[genotypedinpedigree,genotypedinpedigree]
  genotyped=genotypedinpedigree

  Agg=matrix(NA,ncol(G),nrow(G))
  Agg=A[genotyped,genotyped]
  Aggi <- ginv(Agg)

  meanG=mean(G)
  meandiagG=mean(diag(G))
  meanAgg=mean(Agg)
  meandiagAgg=mean(diag(Agg))
  cat("Means G is:",meanG,";Means G diag is:",meandiagG,";Means Agg is:",meanAgg,";Means diag Agg is:",meandiagAgg,"\n")
  beta=(meandiagAgg-meanAgg)/(meandiagG-meanG)
  alpha=meandiagAgg-meandiagG*beta
  cat("alpha and beta value:",alpha,beta,"\n")
  # G=alpha+beta*G

  Aggi=solve(Agg)
  H=matrix(NA,ncol(A),nrow(A))
  colnames(H)=colnames(A)
  rownames(H)=rownames(A)
  H[genotyped,genotyped]=G #H22
  H[nongenotyped,genotyped]=A[nongenotyped,genotyped]%*%(Aggi%*%G) # H12
  H[genotyped,nongenotyped]=t(H[nongenotyped,genotyped]) # H21
  H[nongenotyped,nongenotyped]=A[nongenotyped,nongenotyped] +A[nongenotyped,genotyped]%*%(Aggi%*%(G-Agg)%*%Aggi)%*%A[genotyped,nongenotyped] # H11
  cat("H diag means is :",mean(diag(H)),"H means is:",mean(H),"\n")
  return(H)
}


