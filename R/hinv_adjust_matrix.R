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
#'                         round(hinv_adjust_matrix(M1,ped),2)
#'
#'

hinv_adjust_matrix <- function(M012,ped_full,wts=c(0.95,0.05,1,1)){
  if (length(wts) != 4) stop("You need 4 wts (a, b, tau, omega) in that order")
  a = wts[1]
  b  = wts[2]
  tau   = wts[3]
  omega = wts[4]

  Time = proc.time()

  cat(rep("******",10),"\n")
  cat("G* = a*G + b*A22, a is ",a,"; b is ",b,"\n")
  cat("iG = tau*G - omega*A22, tau is ",tau,"; omega is ",omega,"\n\n\n")

  library(MASS)
  library(sommer)
  library(nadiv)

  Timex = proc.time() # begin
  cat("Begin to build A matrix... \n") # begin

  A <- as.matrix(makeA(prepPed(ped_full)))
  iA <- solve(A)
  rownames(iA) <- colnames(iA) <- rownames(A)

  A = makeA(prepPed(ped_full))
  id=colnames(A)
  A <- as.matrix(A)
  iA <- solve(A)
  rownames(iA) = colnames(iA) = rownames(A) = colnames(A) = id

  Timex = as.matrix(proc.time() - Timex) #end
  cat("\n", "A matrix takes time =", Timex[3]/60, " minutes \n\n\n") #end

  Timex = proc.time() # begin
  cat("Begin to build G matrix... \n") # begin
  G <- A.mat(M012-1)

  Timex = as.matrix(proc.time() - Timex) #end
  cat("\n", "G matrix takes time =", Timex[3]/60, " minutes \n\n\n") #end
  # diag(G) <- diag(G) + 0.01
  genotyped=rownames(G)

  diagG <- diag(G)

  cat(rep("******",10),"\n")
  options(digits = 2)
  cat("Frequency - Diagonal of G\n","\tN:\t\t",length(diagG),"\n",
      "\tMean:\t\t",mean(diagG),"\n",
      "\tMin:\t\t",min(diagG),"\n",
      "\tMax:\t\t",max(diagG),"\n",
      "\tRange:\t\t",range(diagG),"\n")
  options(digits = 7)

  inpedigree=colnames(A)
  nongenotyped=setdiff(inpedigree,genotyped)

  cat("\n\nIn pedigree nongenotyped length: ",length(nongenotyped),"\n")

  genotypednotinpedigree=setdiff(genotyped,inpedigree)
  cat("genotyped not in pedigree",length(genotypednotinpedigree),"\n")

  genotypedinpedigree=intersect(genotyped,inpedigree)
  cat("genotyped in pedigree",length(genotypedinpedigree),"\n")

  genotypedinpedigree <- as.character(genotypedinpedigree)

  G=G[genotypedinpedigree,genotypedinpedigree]
  genotyped=genotypedinpedigree


  genotype <- as.character(genotyped)
  nongenotyped <- as.character(nongenotyped)

  A22 <- A[genotype,genotype]
  diagA22 <- diag(A22)
  offdiagA22 <- c(A22[upper.tri(A22)],A22[lower.tri(A22)])
  options(digits = 2)
  cat("\n\nStatistic of Rel Matrix A22\n","\t\t","N","\t","Mean","\t","Min","\t","Max","\t","Var","\n",
      "Diagonal\t",length(diagA22),"\t",mean(diagA22),"\t",min(diagA22),"\t",max(diagA22),"\t",var(diagA22),
      "\nOff-diagonal\t",length(offdiagA22),"\t",mean(offdiagA22),"\t",min(offdiagA22),"\t",max(offdiagA22),"\t",var(offdiagA22),"\n\n")
  options(digits = 7)


  diagG <- diag(G)
  offdiagG <- c(G[upper.tri(G)],G[lower.tri(G)])
  options(digits = 2)
  cat("Statistic of Genomic Matrix\n","\t\t","N","\t","Mean","\t","Min","\t","Max","\t","Var","\n",
      "Diagonal\t",length(diagG),"\t",mean(diagG),"\t",min(diagG),"\t",max(diagG),"\t",var(diagG),
      "\nOff-diagonal\t",length(offdiagG),"\t",mean(offdiagG),"\t",min(offdiagG),"\t",max(offdiagG),"\t",var(offdiagG),"\n\n")
  options(digits = 7)

  cat("Correlation of Genomic Inbreeding and Pedigree InbreedingA22\n",
      "\tCorrelation:",cor(diagA22,diagG),"\n\n\n")
  cat("Correlation Off-diagonal G and A22\n",
      "\tCorrelation:",cor(offdiagA22,offdiagG),"\n\n\n")

  iA22 <- solve(A22)

  cat(rep("******",10),"\n")
  cat("G and A22 Inv matrix\n")
  cat(rep("******",10),"\n")

  meanoffdiagG=mean(offdiagG)
  meandiagG=mean(diagG)
  meanoffdiagA22=mean(offdiagA22)
  meandiagA22=mean(diagA22)
  cat("Means_off_diag\t","Means_diag:\n","G\t",meanoffdiagG,"\t",meandiagG,
      "\nA22\t",meanoffdiagA22,"\t",meandiagA22,"\n")
  beta=(meandiagA22 - meanoffdiagA22)/(meandiagG - meanoffdiagG)
  alpha=meandiagA22-meandiagG*beta
  cat("Adjust G, and the value of alpha and beta is:",alpha,beta,"\n\n\n")
  G=alpha+beta*G #

  G = a*G + b*A22

  Timex = proc.time() # begin
  cat("Begin to inverse G matrix... \n\n") # begin

  iG <- solve(G)
  rownames(iG) = colnames(iG) = genotyped

  Timex = as.matrix(proc.time() - Timex) #end
  cat("\n", "Inverse G matrix takes time =", Timex[3]/60, " minutes \n\n\n") #end

  iG1 <- iG[genotype,genotype]
  rownames(iA22) = colnames(iA22) <- row.names(A22)
  x22 <- tau*iG1 - omega*iA22
  diagx22 <- diag(x22)
  offdiagx22 <- c(x22[upper.tri(x22)],x22[lower.tri(x22)])
  options(digits = 2)
  cat("Statistic of iG - iA22 Matrix\n","\t\t","N","\t","Mean","\t","Min","\t","Max","\t","Var","\n",
      "Diagonal\t",length(diagx22),"\t",mean(diagx22),"\t",min(diagx22),"\t",max(diagx22),"\t",var(diagx22),
      "\nOff-diagonal\t",length(offdiagx22),"\t",mean(offdiagx22),"\t",min(offdiagx22),"\t",max(offdiagx22),"\t",var(offdiagx22),"\n\n")
  options(digits = 7)


  iH11 <- iA[nongenotyped,nongenotyped]
  iH21 <- iA[genotype,nongenotyped]
  iH12 <- t(iH21)
  iH22 <-  iA[genotype,genotype] + x22
  Hinv <- cbind(rbind(iH11,iH21),rbind(iH12,iH22))

  Time = as.matrix(proc.time() - Time) #end
  cat("\n", "hinv_matrix completed! total time =", Time[3]/60, " minutes \n\n\n") #end
  return(Hinv)
}
