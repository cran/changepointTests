#'@title Function toperform changepoint test for the copula with multiplier bootstrap using for changepoint the new sequential process of Bucher et al (2014)
#'
#'@description This function compute the Cramer-von Mises and Kolmogorov-Smirnov test statistics based on the new sequential process of Bucher et al (2014), using multipliers and parallel computing.
#'Two methods of bootstrapping are used: non-sequential (fastest) and sequential. Both methods yields basically the same P-valueas.
#'
#'@param x  (n x d) matrix of data (observations or pseudo-observations, including residuals), d >=2
#'@param N number of multipliers samples to compute the P-value
#'@param n_cores number of cores for parallel computing (default = 2)
#'@param method 'nonseq' (default) or 'seq'
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
#'@references Nasri, B. R. Remillard, B., & Bahraoui, T. (2021). Change-point problems for multivariate time series using pseudo-observations
#'@references Bucher, A., Kojadinovic, I., Rohmer, T., & Segers, J. (2014). Detecting changes in cross-sectional dependence in multivariate time series, J. Multiv. Anal., 132, 111--128.
#'
#'@examples
#'x<-matrix(rnorm(100),ncol=2)
#'out = test.change.point.copula.BKRS(x)
#'
#'
#'

test.change.point.copula.BKRS <- function(x, N=1000, n_cores=2, method = 'nonseq',est = FALSE){



  tauKS <- NA
  tauCVM <- NA
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)

  U = apply(x,2,rank)/n

  # Computes the D stat in BKRS (2014)

  out0 <- cpCopDStats(U,n,d)

 statT <- out0$statT
 statS <- out0$statS
  KS <-  max(statT)
  CVM <- max(statS)
  ##
  if(est){
    tauKS  = 1 + order(statT,decreasing=T)[1]
    tauCVM = 1 + order(statS,decreasing=T)[1]
  }
  # Generates centered indicator functions and gradient
  out1 <-cpCopStatsBKRS(U,n,d)

  grad = out1$grad
  if(method == 'nonseq'){
    MC  = out1$MC
    MC1 = out1$MC1
  }
  ## Parallel bootstrapping
  # n_cores <- max(2,detectCores()-2);
   cl <- makeCluster(n_cores)
   registerDoParallel(cl)

   fun <- c('multiplierDNonSeqfun','cpCopMultStatsBKRSNonSeq','multiplierDSeqfun','cpCopMultStatsBKRSSeq')

  s <- c(1:n)/n
if(method == 'nonseq'){
  result <- foreach(i=1:N, .export=fun) %dopar% multiplierDNonSeqfun(MC,MC1,grad,s,n,d)
}else{
  result <- foreach(i=1:N, .export=fun) %dopar% multiplierDSeqfun(U,grad,n,d)
  }

  cvm_sim <- rep(0,N)
  ks_sim  <- rep(0,N)
  for (i in 1:N){
    cvm_sim[i] <- result[[i]]$cvm
    ks_sim[i]  <- result[[i]]$ks
   }
  stopCluster(cl)






  pvalueKS  <- 100*mean (ks_sim  > KS)
  pvalueCVM <- 100*mean (cvm_sim > CVM)

  return(list(CVM=CVM, KS=KS, pvalueCVM = pvalueCVM, pvalueKS=pvalueKS, tauCVM=tauCVM, tauKS=tauKS, method=method,NBoot=N, process = 'BKRS with multipliers'))
}
