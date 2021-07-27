#'@title Function to perform traditional bootstrap for changepoint
#'
#'@description This function select a bootstrap sample of indices and compute the Cramer-von Mises and Kolmogorov-Smirnov test statistics for changepoint
#'@param X        n x d matrix of pseudo-observations;
#'@param n        length of the series.
#
#'
#'
#'@return \item{cvm}{simulated value of the Cramer-von Mises statistic}
#'@return \item{ks}{simulated value of the Kolmogorov-Smirnov statistic}
#'
#'
#'@keywords internal
#'
#'@export
#'
bootstrapfun <- function(X,n){

 ind <- sample(c(1:n),n,replace=TRUE)

 X1  <- X[ind,]
out0 <- cpCopStats(X1)
 simKS  <- max(out0$statT)
 simCVM <- max(out0$statS)


 out = list(cvm=simCVM,ks=simKS)
 return(out)


}
