#'@title Function to compute the statistics for the BKRS process
#'
#'@description This function is used to compute the empirical process used in the changepoint tests
#'@param x        n x d matrix of pseudo-observations
#'@param n        length of the series
#'@param d        number of variables
#
#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{MC}{Matrix needed for multipliers}
#'@return \item{MC1}{Matrices needed for multipliers}
#'@return \item{grad1}{Esatimated gradient of the copula}
#'

#'
#'@keywords internal
#'
#'@export
#'

cpCopStatsBKRS = function(x,n,d)
{



  out0 = .C("cpCopulaStatsBucher",
            as.double(x),
            as.integer(n),
            as.integer(d),
            MC = double(n*n),
            MC1 = double(n*n*d),
            grad = double(n*d),
            PACKAGE = "changepointTests"
  )


  #out


}
