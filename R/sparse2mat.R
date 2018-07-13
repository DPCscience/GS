#' Generate inv to mat, use DMU formart to asreml format
#'
#' @param x sparse inv that have three column:id1 id2 value
#' @return matrix
#' @examples
#' library(GS)
#' dmu_inv <- data.frame(ID1 = rep(1:3,each=3),ID2 = rep(1:3,3),value=rnorm(9)+5)
#' dmu_inv
#' asr_mat <- sparse2mat(dmu_inv)
#' asr_inv <- write_mat_to_inv(asr_mat,"none")
#' asr_inv

sparse2mat <- function (x) {
  nrow <- max(x[, 1])
  ncol <- max(x[, 2])
  y <- rep(0, nrow * ncol)
  y[(x[, 2] - 1) * nrow + x[, 1]] <- x[, 3]
  y[(x[, 1] - 1) * nrow + x[, 2]] <- x[, 3]
  matrix(y, nrow = nrow, ncol = ncol, byrow = FALSE)
}
