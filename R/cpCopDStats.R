#'@title Function to compute changepoint statistics
#'
#'@description This function is used to compute the Cramer-von Mises and Kolmogorov-Smirnov test statistics for changepoint
#'@param U        n x d matrix of pseudo-observations
#'@param n        length of the series
#'@param d        number of variables
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
cpCopDStats = function(U,n,d)
{



  out0 = .C("cpChangePointDStat",
            as.double(U),
            as.integer(n),
            as.integer(d),
            statS = double(n),
            statT = double(n),
            PACKAGE = "changepointTests"
  )


  #out


}
