#'@title Function to perform changepoint tests with multiplier bootstrap using the usual sequential process
#'
#'@description This function compute the Cramer-von Mises and Kolmogorov-Smirnov test statistics based on the new sequential process of Bucher et al (2014), using multipliers and parallel computing.
#'
#'@param x  (n x d) matrix of data (observations or pseudo-observations, including residuals), d>=1
#'@param N   number of multipliers samples to compute the P-value
#'@param n_cores number of cores for parallel computing (default = 2)
#'@param boot.method bootstrapping method: 'multipliers' (default, fastest) or 'bootstrap'
#'@param est   if TRUE, tau is estimated (default = FALSE)

#'
#'@author Bouchra R Nasri  and Bruno N Remillard, August 6, 2020
#'
#'@return \item{CVM}{Cramer-von Mises statistic}
#'@return \item{KS}{Kolmogorov-Smirnov statistic}
#'@return \item{pvalueCVM}{Pvalue for the Cramer-von Mises statistic}
#'@return \item{pvalueKS}{Pvalue for theKolmogorov-Smirnov statistic}
#'@return \item{tauCVM}{Estimated changepoint using the Cramer-von Mises statistic}
#'@return \item{tauKS}{Estimated changepoint using the Kolmogorov-Smirnov statistic}
#'
#'@references Nasri, B. R. Remillard, B., & Bahraoui, T. (2022). Change-point problems for multivariate time series using pseudo-observations, J. Multivariate Anal., 187, 104857.
#'
#'@examples
#'x=matrix(rnorm(600),ncol=3)
#'out = test.change.point(x)
#'
#'

test.change.point <- function(x, N=1000, n_cores=2, boot.method = 'multipliers', est = FALSE){

  tauKS <- NA
  tauCVM <- NA
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)

  out0 <- cpCopStats(x)

 MC <- out0$M
 statT <- out0$statT
 statS <- out0$statS
  KS <-  max(statT)
  CVM <- max(statS)
  ##
  if(est){
    tauKS  = 1 + order(statT,decreasing=T)[1]
    tauCVM = 1 + order(statS,decreasing=T)[1]
  }


  ## Parallel bootstrapping
   #n_cores <- max(2,detectCores()-2);
   cl <- makeCluster(n_cores)
   registerDoParallel(cl)
   cvm_sim <- rep(0,N)
   ks_sim  <- rep(0,N)
   if(boot.method== 'multipliers')
     {
                fun <- c('multiplierfun','cpCopMultStats')

                s <- c(1:n)/n

             result <- foreach(i=1:N, .export=fun, .packages = "changepointTests") %dopar% multiplierfun(MC,s,n)

   }else
     {
     fun <- c('bootstrapfun','cpCopStats')

     result <- foreach(i=1:N, .export=fun, .packages = "changepointTests") %dopar% bootstrapfun(x,n)

     }

     stopCluster(cl)

     for (i in 1:N)
     {
       cvm_sim[i] <-  result[[i]]$cvm
       ks_sim[i]  <-  result[[i]]$ks

     }

  pvalueKS  <- 100*mean (ks_sim  > KS)
  pvalueCVM <- 100*mean (cvm_sim > CVM)

  return(list(CVM=CVM, KS=KS, pvalueCVM = pvalueCVM, pvalueKS=pvalueKS, tauCVM=tauCVM, tauKS=tauKS, method= boot.method,NBoot=N))
}

