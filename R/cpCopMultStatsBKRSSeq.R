#'@title Function to compute the statistics for the BKRS process
#'
#'
#'@description This function is used to compute bootstrapped statistics from the BKRS process in the sequential case
#'@param U      matrix needed for multipliers
#'@param grad   gradient of the copula
#'@param xi      multipliers
#'@param n       length of the series
#'@param d      number of variables
#
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


cpCopMultStatsBKRSSeq = function(U,grad,xi,n,d)
{



  out0 = .C("cpCopulaStatsMultBucherSeq",
            as.double(U),
            as.double(grad),
            as.double(xi),
            as.integer(n),
            as.integer(d),
            statS = double(n),
            statT = double(n),
            PACKAGE = "changepointTests"
  )


  #out


}
# dyn.unload('changepointCopR.dll')
