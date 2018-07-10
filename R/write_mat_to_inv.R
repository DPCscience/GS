#' Change the G-matrix ginv that contains three columns format
#'
#' @param x A matrix that is G-matrix
#' @return Three column that contains row, column, and the value of a matrix
#' @examples
#' rMinv <- matrix(1:25,5,5)
#' rMinv
#' res <- write_mat_to_inv(rMinv,type = "none")
#' res


write_mat_to_inv <- function (x, type = c("ginv", "inv", "none"), digits = 10){
  library(MASS)
  type <- match.arg(type)
  if (type == "ginv") rMinv <- ginv(x)
  if (type == "none") rMinv <- x
  rMinv <- round(rMinv, digits)
  res <- data.frame(Row = rep(1:nrow(rMinv), nrow(rMinv)),
                    Column = rep(1:nrow(rMinv), each = nrow(rMinv)),
                    coeff = as.numeric(rMinv),
                    lower = as.logical(lower.tri(rMinv, diag = TRUE)))
  rm(rMinv)
  res <- res[res$lower == TRUE, c("Row", "Column", "coeff")]
  res <- res[order(res$Row, res$Column), ]
  res <- res[res$coeff != 0, ]
    attr(res, "rowNames") <- rownames(x)
    rm(x)
    return(res)
}
