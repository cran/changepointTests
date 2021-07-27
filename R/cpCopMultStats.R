#'@title Function to compute bootstrapped changepoint statistics
#'
#'@description This function compute bootstrapped multipliers values for the Cramer-von Mises and Kolmogorov-Smirnov test statistics for changepoint
#'@param MC      matrix needed for multipliers
#'@param xi      multipliers
#'@param s       sequence of normalized values in (0,1)
#'@param n       length of the series
#
#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{statS}{Simulated values of the Cramer-von Mises statistics}
#'@return \item{statT}{Simulated values of the Kolmogorov-Smirnov statistics}
#'

#'
#'@keywords internal
#'
#'@export
#'
cpCopMultStats = function(MC,xi,s,n)
{


  out0 = .C("cpCopulaStatsMult",
            as.double(MC),
            as.double(xi),
            as.double(s),
            as.integer(n),
            statS = double(n),
            statT = double(n),
            PACKAGE = "changepointTests"
  )
}
