#'@title Function to compute the empirical cdf at given points
#'
#'
#'@description This function is used to compute the empirical process used in the changepoint tests
#'@param x        n x d matrix of pseudo-observations
#'@param u        evaluation points of the cdf
#
#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{M}{Values of the empirical process}
#'@return \item{cumsum}{Associated cumulativ sums}
#'

#'
#'@keywords internal
#'
#'@export
#'
emp.cdf = function(x,u)
{
  n = nrow(x)
  d = ncol(x)


  out0 = .C("empcdf",
            as.double(x),
            as.integer(n),
            as.integer(d),
            as.double(u),
            M = double(n),
            cumsum = double(n),
            PACKAGE = "changepointTests"
)


  #out


}





