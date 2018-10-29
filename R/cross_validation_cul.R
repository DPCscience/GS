#' Get the cross validation, generate five data, 0 as missing value
#'
#' @param dat
#' @param trait the trait want to cross vation
#' @param seed set the seed number, default is 123
#' @return finve files in the working directory
#' @examples
#' 1+1

cross_validation_cul <- function(h2=1,trait){
  se <- function(x){
    sd(x)/sqrt(length(x))
  }
  # h2 =1
  library(dplyr)
  dat = as.data.frame(fread("dat.csv"))
  accur = NULL;unbia=NULL
  names(dat)[1] = "ID"
  dat$ID = as.character(dat$ID)

  t1 = as.data.frame(fread("dat1/dat.sln"))
  t1$n = 1:dim(t1)[1]
  hang = grep("animal",t1$Model_Term)[1]
  t1 = t1[-c(1:hang),]
  select(dat,Level=ID,trait) %>% merge(t1,by="Level") ->a
  cor(a[,2],a[,4])
  dat1 = as.data.frame(fread("dat1.csv"))
  id1 = as.character(dat1[dat1[,trait]==0,]$animal)
  b = a[a$Level %in% id1,]
  accur[1] = cor(b[,trait],b$Effect)/sqrt(h2)
  unbia[1] = lm(b[,trait] ~ b$Effect)$coefficients[[2]]

  t2 = as.data.frame(fread("dat2/dat.sln"))
  t2$n = 1:dim(t2)[1]
  ang = grep("animal",t2$Model_Term)[1]
  t2 = t2[-c(1:hang),]
  select(dat,Level=ID,trait) %>% merge(t2,by="Level") ->a
  dat2 = as.data.frame(fread("dat2.csv"))
  id2 = as.character(dat2[dat2[,trait]==0,]$animal)
  b = a[a$Level %in% id2,]
  accur[2] = cor(b[,trait],b$Effect)/sqrt(h2)
  unbia[2] = lm(b[,trait] ~ b$Effect)$coefficients[[2]]

  t3 = as.data.frame(fread("dat3/dat.sln"))
  t3$n = 1:dim(t3)[1]
  ang = grep("animal",t3$Model_Term)[1]
  t3 = t3[-c(1:hang),]
  select(dat,Level=ID,trait) %>% merge(t3,by="Level") ->a
  dat3 = as.data.frame(fread("dat3.csv"))
  id3 = as.character(dat3[dat3[,trait]==0,]$animal)
  b = a[a$Level %in% id3,]
  accur[3] = cor(b[,trait],b$Effect)/sqrt(h2)
  unbia[3] = lm(b[,trait] ~ b$Effect)$coefficients[[2]]


  t4 = as.data.frame(fread("dat4/dat.sln"))
  t4$n = 1:dim(t4)[1]
  ang = grep("animal",t4$Model_Term)[1]
  t4 = t4[-c(1:hang),]
  select(dat,Level=ID,trait) %>% merge(t4,by="Level") ->a
  dat4 = as.data.frame(fread("dat4.csv"))
  id4 = as.character(dat4[dat4[,trait]==0,]$animal)
  b = a[a$Level %in% id4,]
  accur[4] = cor(b[,trait],b$Effect)/sqrt(h2)
  unbia[4] = lm(b[,trait] ~ b$Effect)$coefficients[[2]]

  t5 = as.data.frame(fread("dat5/dat.sln"))
  t5$n = 1:dim(t5)[1]
  ang = grep("animal",t5$Model_Term)[1]
  t5 = t5[-c(1:hang),]
  select(dat,Level=ID,trait) %>% merge(t5,by="Level") ->a
  dat5 = as.data.frame(fread("dat5.csv"))
  id5 = as.character(dat5[dat5[,trait]==0,]$animal)
  b = a[a$Level %in% id5,]
  accur[5] = cor(b[,trait],b$Effect)/sqrt(h2)
  unbia[5] = lm(b[,trait] ~ b$Effect)$coefficients[[2]]

  n1 <- mean(accur);n1e <- se(accur)
  n3 <- mean(unbia);n3e <- se(unbia)
  ax <- data.frame(type=c("r_pearson","r_unbiased"),
                   value = c(n1,n3),
                   se = c(n1e,n3e))
  ax
  return(ax)
}
