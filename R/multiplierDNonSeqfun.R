#'@title Function to perform multiplier bootstrap for changepoint
#'
#'@description This function simulates a random sample of Gaussian multipliers  null hypothesis of a Gaussian HMM
#'and compute the Cramer-von Mises and Kolmogorov-Smirnov test statistics.
#'
#'@param MC      matrix needed for multipliers
#'@param MC1     matrices needed for multipliers
#'@param grad   gradient of the copula
#'@param s       sequence of normalized values in (0,1)
#'@param n       length of the series
#'@param d      number of variables
#
#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{cvm}{simulated value of the Cramer-von Mises statistic}
#'@return \item{ks}{simulated value of the Kolmogorov-Smirnov statistic}
#'
#'
#'@references Nasri, B. R. Remillard, B., & Bahraoui, T. (2022).
#   Change-point problems for multivariate time series using pseudo-observations
#'
#'@keywords internal
#'
#'@export
#'
multiplierDNonSeqfun <- function(MC,MC1,grad,s,n,d){
  xi <- rnorm(n)

 out0 <- cpCopMultStatsBKRSNonSeq(MC,MC1,grad,xi,s,n,d)

  simKS  <- max(out0$statT)
  simCVM <- max(out0$statS)


 out = list(cvm=simCVM,ks=simKS)
  return(out)

}
