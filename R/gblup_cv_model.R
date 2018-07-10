#' Get the arrucacy and other values when run crossvalidation
#'
#' @param pheno_data contins two columns: ID and y
#' @param inv contains three columns that is ginv or ainv: row, col, inverse_values
#' @return result
#' @examples
#'

gblup_cv_model <- function(pheno_data,inv){
  library(asreml)
  dd <- pheno_data
  ainv <- inv
  names(dd) <- c("ID","y")
  dd$ID <- as.factor(as.character(dd$ID))
  row.names(dd) <- dd$ID
  library(caret)
  w <- createFolds(1:length(dd$ID),k = 5)
  r_pearson <- NULL
  r_spearman <- NULL
  r_unbiased <- NULL
  r_MSD <- NULL
  for(i in 1:5){
    vv <- dd$ID[w[[i]]]
    vv <- as.character(vv)
    tt <- dd
    tt[vv,]$y <- NA #将测试集观测值为NA
    mod <- asreml(y ~ 1, random=~ giv(ID),ginverse= list(ID=ainv),data=tt,workspace=1e8)
    gblup <- predict(mod,"ID")$predictions$pval
    rownames(gblup) <- gblup$ID
    r_pearson[i] <- cor(gblup[vv,]$predicted.value,dd[vv,]$y,method = "pearson")
    r_spearman[i] <- cor(gblup[vv,]$predicted.value,dd[vv,]$y,method = "spearman")
    r_unbiased[i] <- summary(lm(dd[vv,]$y~gblup[vv,]$predicted.value-1))$coefficients[1]
    r_MSD[i] <- mean((dd[vv,]$y-gblup[vv,]$predicted.value)^2)
  }
  se <- function(x){
    sd(x)/sqrt(length(x))
  }
  n1 <- mean(r_pearson);n1e <- se(r_pearson)
  n2<- mean(r_spearman);n2e <- se(r_spearman)
  n3 <- mean(r_unbiased);n3e <- se(r_unbiased)
  n4 <- mean(r_MSD);n4e <- se(r_MSD)
  ax <- data.frame(type=c("Accuracy","Rank_corr","Unbiased","Error_means_squre"),
                   value = c(n1,n2,n3,n4),
                   se = c(n1e,n2e,n3e,n4e))
  return(ax)
}
