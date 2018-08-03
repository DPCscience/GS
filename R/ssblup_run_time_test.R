#' Change the Siminute SNP data to 012 format
#'
#' @param geno A format that has row.name(ID), and first column is first allet of SNP(1,2)
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' gg <- data.frame(V1 = c(1,1,2,2,1,2,1,2),V2 = c(1,1,2,2,1,2,1,2),V3 = c(2,1,2,2,1,1,2,1),V4=round(runif(8,1,2)),V5=round(runif(8,1,2)),V6=round(runif(8,1,2)))
#' row.names(gg) <- paste("ID",1:8,sep="")
#' gg
#' snp_trans_012(gg)

ssblup_run_time_test<- function(){
#! /opt/R-3.4.4/bin/Rscript
library(data.table)
library(GS)
dat <- fread("p1_data_001.txt")
dat <- dat[,c(1:3,10)]
ped_full <- dat[,1:3]
dat <- dat[,c(1,4)]
names(dat) <- c("ID","y")
dat$ID <- as.factor(dat$ID)
cat("\n", "Phenotype Data total rows:", dim(dat)[1], "\n\n") #end
Mak <- data.frame(fread("p1_mrk_001.txt",skip = 1),row.names = 1)
cat("\n", "Genotype total rows and 2*Numbewr columns:", "\n\n")
dim(Mak)

cat("\n", "Look at the Markers: 5 rows and 5 columns", "\n")
Mak[1:5,1:5]

library(GS)
M <- GS::snp_trans_012(Mak) # 4000 24000 # Completed! Time = 3.180167  minutes
H <- H_matrix(M,ped_full)
Hinv <- hinv_matrix(M,ped_full)
}
