#' Build a matrix base the column vector
#'
#' @param column Data with one column
#' @return The matrix
#' @examples
#' 1+1


build_matrix <- function(column){
  v = column
  if(is.numeric(v)){
    va=v
    mrow=length(va)
    mcol=max(va)
  }
  if(is.character(v)){
    vf = factor(v)
    va = as.numeric(vf)
    mrow = length(va)
    mcol = length(levels(vf))
  }
  X = matrix(data=c(0),nrow = mrow, ncol=mcol)

  for(i in 1:mrow){
    ic = va[i]
    X[i,ic] =1
  }
  return(X)
}
