#'@title Function to compute the bootstrapped statistics for the BKRS process
#'
#'@description This function is used to compute bootstrapped statistics from the BKRS process in the non-sequential case
#'@param MC      matrix needed for multipliers
#'@param MC1     matrices needed for multipliers
#'@param grad   gradient of the copula
#'@param xi      multipliers
#'@param s       sequence of normalized values in (0,1)
#'@param n       length of the series
#'@param d        number of variables
#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{statS}{Values of the Cramer-von Mises statistics}
#'@return \item{statT}{Values of the Kolmogorov-Smirnov statistics}
#'

#'
#'@keywords internal
#'
#'@export
#'


cpCopMultStatsBKRSNonSeq = function(MC,MC1,grad,xi,s,n,d)
{



  out0 = .C("cpCopulaStatsMultBucherNonSeq",
            as.double(MC),
            as.double(MC1),
            as.double(grad),
            as.double(xi),
            as.double(s),
            as.integer(n),
            as.integer(d),
            statS = double(n),
            statT = double(n),
            PACKAGE = "changepointTests"
  )


  #out


}
