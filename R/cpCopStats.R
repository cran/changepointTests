#'@title Function to compute the statistics for the traditional empirical process
#'
#'@description This function is used to compute the empirical process used in the changepoint tests
#'@param x        n x d matrix of pseudo-observations
#
#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{M}{Matrix for the cdf}
#'@return \item{statS}{Values of the Cramer-von Mises statistics}
#'@return \item{statT}{Values of the Kolmogorov-Smirnov statistics}
#'

#'
#'@keywords internal
#'
#'@export
#'


cpCopStats = function(x)
{
  n = nrow(x)
  d = ncol(x)


  out0 = .C("cpCopulaStats",
            as.double(x),
            as.integer(n),
            as.integer(d),
            M = double(n*n),
            statS = double(n),
            statT = double(n),
            PACKAGE = "changepointTests"
  )


  #out


}
