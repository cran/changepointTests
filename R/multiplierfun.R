#'@title Function to perform multiplier bootstrap for changepoint
#'
#'@description This function simulates a random sample of Gaussian multipliers  null hypothesis of a Gaussian HMM
#'and compute the Cramer-von Mises and Kolmogorov-Smirnov test statistics.
#'
#'@param MC         n x n matrix = MM  - C, with MM[i,j] = 1(Xi <= Xj) and C=mean(M[,j]);
#'@param n          length of the series.
#
#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{cvm}{simulated value of the Cramer-von Mises statistic}
#'@return \item{ks}{simulated value of the Kolmogorov-Smirnov statistic}
#'
#'@references Chapter 8 of B. Remillard (2013). Statistical Methods for Financial Engineering,
#'Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
#'
#'
#'@keywords internal
#'
#'@export
#'
multiplierfun <- function(MC,s,n){
  xi <- rnorm(n)

 out0 <- cpCopMultStats(MC,xi,s,n)

  simKS  <- max(out0$statT)
  simCVM <- max(out0$statS)


 out = list(cvm=simCVM,ks=simKS)
  return(out)

}
