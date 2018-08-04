#' Build the H-inverse-matrix from G-inverse-matrix and pedigree-full and pedigree-genotype#'
#' @param G_mat It is the matrix which has rownames and colnames(ID)
#' @param ped_full It contains the full pedigree, it has three columns:ID,Sire,Dam
#' @param ped_geno It contains the pedigree that has the genotype, it is a part of the ped_full pedigree
#' @return The H-inverse-matrix form the formula
#' @examples
#' # 这里将iA22是A矩阵提取的A22，进行求逆
#' # paper:Legarra A, Aguilar I, Misztal I. A relationship matrix including full pedigree and genomic information.[J]. Journal of Dairy Science, 2009, 92(9):4656-63.
#'ped_full <- data.frame(ID=9:17,Sire=c(1,3,5,7,9,11,11,13,13),Dam=c(2,4,6,8,10,12,4,15,14))
#'ped_full
#'G <- matrix(0.7,4,4)
#'diag(G) <- 1
#'rownames(G) <- colnames(G) <- 9:12
#'G
#'#another example
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
#'                         ped_full    <- data.frame(animal, sire, dam)
#'                         rm(list=c("animal","dam","sire"))
#'                         M <- data.11.1[6:14, c(1, 7:16)]
#'                         rownames(M) <- M[, 1]
#'                         M012 <- as.matrix(M[, -1])
#'                         hinv_matrix(M012,ped_full) %>% round(3)



hinv_matrix_asreml <- function(M012,ped_full){

  Time = proc.time() # begin
  cat("Begin to build the Hinv matrix... \n\n") # begin

  library(MASS)
  library(sommer)
  library(asreml)
  library(nadiv)

  Timex = proc.time() # begin
  cat("Begin to build A matrix... \n\n") # begin

  A <- as.matrix(makeA(prepPed(ped_full)))

  Timex = as.matrix(proc.time() - Timex) #end
  cat("\n", "A matrix takes time =", Timex[3]/60, " minutes \n") #end

  ainv <- asreml.Ainverse(ped_full)$ginv
  iA <- asreml.sparse2mat(ainv)

  inpedigree <- attr(ainv,"rowNames")
  row.names(iA) = colnames(iA) = inpedigree

  Timex = proc.time() # begin
  cat("Begin to build G matrix... \n\n") # begin
  G <- A.mat(M012-1)

  Timex = as.matrix(proc.time() - Timex) #end
  cat("\n", "G matrix takes time =", Timex[3]/60, " minutes \n") #end

  diag(G) <- diag(G) + 0.01
  genotyped=rownames(G)

  inpedigree=colnames(A)
  nongenotyped=setdiff(inpedigree,genotyped)

  cat("In pedigree nongenotyped length: ",length(nongenotyped),"\n")

  genotypednotinpedigree=setdiff(genotyped,inpedigree)
  cat("genotyped not in pedigree",length(genotypednotinpedigree),"\n")

  genotypedinpedigree=intersect(genotyped,inpedigree)
  cat("genotyped in pedigree",length(genotypedinpedigree),"\n")

  genotypedinpedigree <- as.character(genotypedinpedigree)

  G=G[genotypedinpedigree,genotypedinpedigree]
  genotyped=genotypedinpedigree


  genotype <- as.character(genotyped)
  nongenotyped <- as.character(nongenotyped)

  Timex = proc.time() # begin
  cat("Begin to inverse G matrix... \n\n") # begin

  iG <- solve(G)
  rownames(iG) = colnames(iG) = genotyped

  Timex = as.matrix(proc.time() - Timex) #end
  cat("\n", "Inverse G matrix takes time =", Timex[3]/60, " minutes \n") #end


  iG1 <- iG[genotype,genotype]
  A22 <- A[genotype,genotype]
  iA22 <- solve(A22)
  # iA22 <- iA[genotype,genotype]

  meanG=mean(G)
  meandiagG=mean(diag(G))
  meanAgg=mean(A22)
  meandiagAgg=mean(diag(A22))
  cat("Means G is:",meanG,";Means G diag is:",meandiagG,";Means A22 is:",meanAgg,";Means diag A22 is:",meandiagAgg,"\n")
  beta=(meandiagAgg-meanAgg)/(meandiagG-meanG)
  alpha=meandiagAgg-meandiagG*beta
  cat("Not adjust, but the value of alpha and beta is:",alpha,beta,"\n")
  # G=alpha+beta*G

  rownames(iA22) = colnames(iA22) <- row.names(A22)
  x22 <- iG1 - iA22
  iH11 <- iA[nongenotyped,nongenotyped]
  iH21 <- iA[genotype,nongenotyped]
  iH12 <- t(iH21)
  iH22 <-  iA[genotype,genotype] + x22
  Hinv <- cbind(rbind(iH11,iH21),rbind(iH12,iH22))

  Time = as.matrix(proc.time() - Time) #end
  cat("\n", "Completed! Time =", Time[3]/60, " minutes \n") #end
  Hinv <- as.matrix(Hinv)
  return(Hinv)
}
