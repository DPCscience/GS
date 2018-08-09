cv_result <- function(dat_six_cloums){
  library(data.table)
  dd <- dat_six_cloums
  head(dd)
  id1 <- dd[dd$y1=="*",]$ID
  id2 <- dd[dd$y2=="*",]$ID
  id3 <- dd[dd$y3=="*",]$ID
  id4 <- dd[dd$y4=="*",]$ID
  id5 <- dd[dd$y5=="*",]$ID

  # pv <- fread("w/tt6/tt.pvs",skip = 13)
  options(warn=-1)
  pv1 <- fread("dat_model1/dat_model.pvs")
  pv2 <- fread("dat_model2/dat_model.pvs")
  pv3 <- fread("dat_model3/dat_model.pvs")
  pv4 <- fread("dat_model4/dat_model.pvs")
  pv5 <- fread("dat_model5/dat_model.pvs")
  pv <- fread("dat_model6/dat_model.pvs")
  options(warn=1)

  py1 <- pv1[pv1$ID %in% id1,][,1:2];py1$ID <- as.factor(py1$ID)
  py2 <- pv2[pv2$ID %in% id2,][,1:2];py2$ID <- as.factor(py2$ID)
  py3 <- pv3[pv3$ID %in% id3,][,1:2];py3$ID <- as.factor(py3$ID)
  py4 <- pv4[pv4$ID %in% id4,][,1:2];py4$ID <- as.factor(py4$ID)
  py5 <- pv5[pv5$ID %in% id5,][,1:2];py5$ID <- as.factor(py5$ID)

  p_list <- list(py1,py2,py3,py4,py5)
  r_pearson <- NULL
  r_unbiased <- NULL
  head(dd)
  dat <- dd[,1:2]
  dat$ID <- as.factor(dat$ID)
  head(dat)
  for(i in 1:5){
    # i=1
    x <- merge(p_list[[i]],dat,"ID")
    names(x) <- c("ID","pv","y")
    x <- na.omit(x)
    r_pearson[i] <- cor(x$pv,x$y,method = "pearson")
    # r_spearman[i] <- cor(x$pv,x$y,method = "spearman")
    r_unbiased[i] <- summary(lm(x$y ~ x$pv -1))$coefficients[1]
    # r_MSD[i] <- mean((x$y - x$pv)^2)
  }
  se <- function(x){
    sd(x)/sqrt(length(x))
  }
  r_pearson
  n1 <- mean(r_pearson);n1e <- se(r_pearson)
  # n2<- mean(r_spearman);n2e <- se(r_spearman)
  n3 <- mean(r_unbiased);n3e <- se(r_unbiased)
  # n4 <- mean(r_MSD);n4e <- se(r_MSD)
  ax <- data.frame(type=c("r_pearson","r_unbiased"),
                   value = c(n1,n3),
                   se = c(n1e,n3e))
  ax
  return(ax)
}
