#' Generate the K-fold number in cross-validation
#'
#' @param y The length of the y variate
#' @param k the k folder, is a number
#' @return K-subset of the length
#' @examples
#' dat <- data.frame(ID=paste("ID",151:160,sep=""),y=rnorm(10))
#' dat
#' index <- cv_generater(1:length(dat$y),k=3)
#' dat$y1 = dat$y2 = dat$y3 =dat$y;
#' index
#' for(i in 1:length(index)){
#'   vv <- dat$ID[index[[i]]]
#'   vv <- as.character(vv)
#'   dat[match(vv,dat$ID),i+2] <- NA
#' }
#' dat
#' index
#'

cv_generater <- function (y, k = n, list = TRUE, returnTrain = FALSE) {
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}



