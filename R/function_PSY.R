#' Data generating process for simulations: random walk with drift 1/n
#'
#' This function simulates a random walk with drift 1/n, which is used when computing
#' the tests for bubbles by Phillips et al. (2015)
#'
#' @param n is the number of simulated time steps
#' @param niter is the number of replications for the simulated DGP
#' @return y a n x niter matrix of simulated data, where each colum represents a single simulated process
#' @details
#' This function simulates a n x 1 vector following a a random walk with drift 1/n,
#' and this operation is repeated niter times
#' @export
#' @importFrom MASS mvrnorm
#-------------------------------------------------------------------
DGP <- function (n, niter){ #random walk with drift 1/n@
  u0<-1/n
  set.seed(niter)
  rn <- MASS::mvrnorm (n, rep(u0, niter), diag( niter))
  y <- apply (rn, 2, cumsum)
  return(y)
}





#' Critical values for the sup ADF (SADF) statistic.
#'
#' This function computes the critical values for the sup ADF statistic by Phillips et al. (2015)
#'
#' @param qe is a q x 1 vector of quantiles
#' @param m is the number of replications for the simulated DGP
#' @param T is the number of simulated time steps
#' @param r0 is the minimum window size fraction
#' @param address is a character object containing the address where the plot density of SADF test stat will be saved in png format
#' @param clust_number is the number of clusters for parallel computation
#' @return quantile_sadf a q x 1 vector of simulated critical values for the sup ADF statistic
#' @details
#' This function computes the critical values for the sup ADF statistic by Phillips et al. (2015) and saves the plot density of SADF test statistic as a png file in a specified address
#' @export
#' @importFrom pbapply pbapply
#' @importFrom urca ur.df
#' @importFrom grDevices dev.copy dev.off png
#' @importFrom graphics plot
#' @importFrom stats density quantile
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel stopCluster
#' @importFrom grDevices dev.copy
#' @importFrom grDevices dev.off
#'
#' @examples
#'
#'  \dontrun{
#'  # Replicate finite sample critical values in Table 1 for SADF by PSY(2015)
#'  qe <-c(0.90,0.95,0.99)
#'  address <- "D:/bubbletest/"
#'  m <- 1000           # number of replications
#'  table1_sadf <- matrix(0,5,length(qe))
#'
#'  table1_sadf[1,] <- CV_SADF(qe,m,100,0.190,address,clust_number=8)
#'  table1_sadf[2,] <- CV_SADF(qe,m,200,0.137,address,clust_number=8)
#'  table1_sadf[3,] <- CV_SADF(qe,m,400,0.100,address,clust_number=8)
#'  table1_sadf[4,] <- CV_SADF(qe,m,800,0.074,address,clust_number=8)
#'  table1_sadf[5,] <- CV_SADF(qe,m,1600,0.055,address,clust_number=8)
#'
#'  write.table( table1_sadf,paste0(address, "table1_sadf.csv"), row.names = F, col.names =F)
#'    }

CV_SADF <- function(qe, m, T, r0, address, clust_number) {
  # Attention: PSY(2015) suggest to use swindow0=floor((0.01+1.8/sqrt(T))*T),
  # but we leave it as an input, because we want to replicate their code, where they used
  # different minimum window sizes
  swindow0<-floor(r0*T)
  set.seed( 1 )
  y <- DGP(T , m)
  cl <- parallel::makeCluster(clust_number)
  parallel::clusterEvalQ(cl, library(urca))

  badfs <- pbapply::pbapply(y, 2, function(col_y) {
    sapply(swindow0 : T, function(i) {
      x <- urca::ur.df(col_y[1:i] , type = "drift", lags = 0, selectlags = "Fixed")
      return(x@teststat[1])
    })
  }, cl= cl)
  parallel::stopCluster(cl)
  sadf <- apply(badfs, 2, max)

  quantile_sadf <- quantile(sadf, probs = qe)

  #plot density of test stat
  d <- density(as.numeric(sadf))
  name<-paste("SADF_","T", T ,"_r", round(r0, digits = 3),".png", sep="")
  name <- paste0(address, name)

  plot(d, main = "SADF", type = "l", col = "blue")
  grDevices::dev.copy(png,filename = name)
  grDevices::dev.off()

  return(quantile_sadf)
}


#' Critical value sequences for the backward ADF (BADF) statistic.
#'
#' This function computes the critical values  for the backward ADF statistic  by Phillips et al. (2015)
#'
#' @param qe is a q x 1 vector of quantiles
#' @param m is the number of replications for the simulated DGP
#' @param T is the number of simulated time steps
#' @param r0 is the minimum window size fraction
#' @param address is a character object containing the address where the plot density of SADF test stat will be saved in png format
#' @param clust_number is the number of clusters for parallel computation
#' @return quantile_badfs a q x 1 vector of simulated critical values for the sup ADF statistic
#' @details
#' This function computes the critical values for the backward ADF statistic by Phillips et al. (2015) and saves the plot density of the test statistic as a png file in a specified address
#' @export
#' @importFrom urca ur.df
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#'
#' @examples
#'
#' \dontrun{
#' qe <-c(0.90,0.95,0.99)
#' address <- "D:/bubbletest/"
#' m <- 1000     # number of replications
#' system.time(aa<-CV_BADF(qe,m,100,0.190,address,clust_number=8))
#' }

CV_BADF <- function(qe, m, T, r0, address, clust_number){
  # Attention: PSY(2015) suggest to use swindow0=floor((0.01+1.8/sqrt(T))*T),
  # but we leave it as an input, because we want to replicate their code, where they used
  # different minimum window sizes
  swindow0=floor(r0*T)
  adfs <- matrix(0, m, (T - swindow0 + 1))
  cl <- parallel::makeCluster(clust_number)
  registerDoParallel(cl)
  r <- foreach ( r2 = swindow0:T, .packages = c('urca', 'bubble')) %dopar%
  {
    set.seed(1)
    y<-DGP(r2,m)
    adf_r<-apply(y, 2, function(x){
      x1 <- urca::ur.df(x, type = "drift", lags = 0, selectlags = "Fixed")
      x1 <- x1@teststat[1]
    })
  }
  parallel::stopCluster(cl)
  for (r2 in 1:(T - swindow0 + 1)) {adfs[,r2] <- r[[r2]][]}

  #plot density of test stat
  d <- density(as.numeric(adfs))
  name<-paste("BADF_","T",T,"_r",round(r0, digits = 3),".png", sep="")
  name <- paste0(address, name)
  plot(d, main ="BADF", type = "l", col = "red")
  dev.copy(png,filename = name)
  dev.off()

  quantile_badfs <- quantile(adfs,qe)
  return(quantile_badfs)
}


#'  Critical values for the backward SADF (BSADF) statistic
#'
#' This function computes the critical values for the backward SADF (BSADF) statistic  by Phillips et al. (2015)
#'
#' @param qe is a q x 1 vector of quantiles
#' @param m is the number of replications for the simulated DGP
#' @param T is the number of simulated time steps
#' @param r0 is the minimum window size fraction
#' @param address is a character object containing the address where the plot density of SADF test stat will be saved in png format
#' @param clust_number is the number of clusters for parallel computation
#' @param beta is the significance level for the critical values for the BSADF sequence
#' @return results is a list containing quantile_bsadf (a q x 1 vector of simulated critical values for the BSADF statistic) and bsadf_cv (a (T-swindow0+1) vector containing the critical values for the BSADF sequence)
#' @details
#' This function computes the simulated critical values for the backward SADF (BSADF) statistic by Phillips et al. (2015) and saves the plot density of the test statistic as a png file in a specified address.
#' The main difference between the SADF and GSADF tests on one side and BADF BSADF tests on the other side, lies in the data generating process:
#' the SADF and GSADF tests use a random walk with an overall drift value of 1/T, while the
#' the BADF BSADF tests use a random walk with an drift value of 1/r2 for each r2
#' @export
#' @importFrom urca ur.df

CV_BSADF <- function(qe, m, T, r0, address, clust_number, beta=0.95){
  # Attention: PSY(2015) suggest to use swindow0=floor((0.01+1.8/sqrt(T))*T),
  # but we leave it as an input, because we want to replicate their code, where they used
  # different minimum window sizes
  swindow0=floor(r0*T)
  Msadfs <- matrix(0 ,m ,(T - swindow0 + 1))

  for (r2 in swindow0 : T)
  {
    set.seed(1)
    y <- DGP( r2 , m)

    badfs <- matrix(0,r2-swindow0+1,m);

    cl <- makeCluster(clust_number) #
    registerDoParallel(cl) #
    r <- foreach ( j = 1 : m, .packages = 'urca' ) %dopar%
    {
      badfs_r <- rep(0,r2-swindow0+1)
      for (r1 in 1:(r2-swindow0+1))
      {
        x<-ur.df(y[r1:r2,j], type = "drift", lags = 0, selectlags = "Fixed")
        badfs_r[r1] <- x@teststat[1]
      }
      badfs[,j] <- badfs_r ##
    }
    stopCluster(cl) #
    for (j in 1:m) {badfs[,j] <- r[[j]][]} ##

    if (r2 == swindow0) {sadfs <- badfs}
    else  {sadfs <- apply( badfs,2,max)}
    Msadfs[,(r2-swindow0+1)] <- sadfs
  }

  # I create the critical values for the BSADF sequence
  bsadf_cv <- apply(Msadfs, 2, quantile, probs = beta,  na.rm = TRUE) ##

  #plot density of test stat
  d <- density(as.numeric(Msadfs))
  name <- paste("BSADF_","T",T,"_r",round(r0, digits = 3),".png", sep="")
  name <- paste0(address, name)
  plot(d, main ="BSADF", type = "l", col = "green")
  dev.copy(png,filename = name)
  dev.off()

  quantile_bsadf <- quantile(Msadfs,qe)
  results <-list(quantile_bsadf=quantile_bsadf, bsadf_cv=bsadf_cv)
  return(results)
}


#'  Critical values for the generalized sup ADF (GSADF) statistic.
#'
#' This function computes the critical values for the generalized sup ADF (GSADF) statistic by Phillips et al. (2015)
#'
#' @param qe is a q x 1 vector of quantiles
#' @param m is the number of replications for the simulated DGP
#' @param T is the number of simulated time steps
#' @param r0 is the minimum window size fraction
#' @param address is a character object containing the address where the plot density of SADF test stat will be saved in png format
#' @param clust_number is the number of clusters for parallel computation
#' @param beta is the significance level for the critical values for the BSADF sequence
#' @return results is a list containing quantile_gsadf (a q x 1 vector of simulated critical values for the GSADF statistic) and gsadf_cv (a (T-swindow0+1) vector containing the critical values for the BSADF sequence)
#' @details
#' This function computes the simulated critical values for the GSADF statistic,  the critical value sequence for the backward SADF (BSADF) statistic sequence by Phillips et al. (2015) and saves the plot density of BADF test statistic as a png file in a specified address
#' The main difference between the SADF and GSADF tests on one side and BADF BSADF tests on the other side, lies in the data generating process:
#' the SADF and GSADF tests use a random walk with an overall drift value of 1/T, while the
#' the BADF BSADF tests use a random walk with an drift value of 1/r2 for each r2
#' @export
#' @importFrom urca ur.df
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#'
#' @examples
#'
#'  \dontrun{
#'  # Replicate Finite sample critical values in Table 1 for GSADF by PSY(2015)
#'  qe <-c(0.90,0.95,0.99)
#'  m <- 1000           # number of replications
#'  address <- "D:/bubbletest/"
#'  clust_number <- 7
#'  beta=0.95
#'
#'  table1_gsadf <- matrix(0,5,length(qe))
#'
#'  table1_gsadf[1,] <- CV_GSADF(qe,m,100,0.190,address, clust_number, beta)$quantile_gsadf
#'  table1_gsadf[2,] <- CV_GSADF(qe,m,200,0.137,address, clust_number, beta)$quantile_gsadf
#'  table1_gsadf[3,] <- CV_GSADF(qe,m,400,0.100,address, clust_number, beta)$quantile_gsadf
#'  table1_gsadf[4,] <- CV_GSADF(qe,m,800,0.074,address, clust_number, beta)$quantile_gsadf
#'  table1_gsadf[5,] <- CV_GSADF(qe,m,1600,0.055,address, clust_number, beta)$quantile_gsadf
#'
#'  write.table( table1_gsadf,paste0(address, "table1_gsadf.csv"), row.names = F, col.names =F)
#'    }

CV_GSADF <- function(qe, m, T, r0, address, clust_number, beta=0.95){
  # Attention: PSY(2015) suggest to use swindow0=floor((0.01+1.8/sqrt(T))*T),
  # but we leave it as an input, because we want to replicate their code, where they used
  # different minimum window sizes
  swindow0=floor(r0*T)

  cl <- makeCluster(clust_number)
  y <- DGP(T,m)
  Msadfs <- matrix(0,m,(T-swindow0+1))##

  registerDoParallel(cl)

  r <- foreach ( j = 1 : m, .packages = 'urca' ) %dopar%
  {
    #####################################
    # This part of code creates a log.txt file
    # in your current working directory and write
    # the number of current iteration in it
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",j,"\n"))
    sink()
    #####################################
    sadfs <- rep(0,T-swindow0+1)
    for (r2 in swindow0:T)
    {
      dim0 <- r2 - swindow0 + 1
      rwadft <- rep( 0 , dim0)
      for (r1 in 1:dim0)
      {
        x <- ur.df(y[(r1:r2),j],type = "drift",lags = 0, selectlags = "Fixed")  #% two tail 5% significant level
        rwadft[r1] <- x@teststat[1]
      }
      sadfs[r2 - swindow0 + 1] <- max(rwadft)
    }
    Msadfs[j,] <- sadfs ##
  }
  stopCluster(cl)

  for (j in 1:m) {Msadfs[j,] <- r[[j]][]} ##

  # I create the critical values for the BSADF sequence
  gsadf_cv <- apply(Msadfs, 2, quantile, probs = beta,  na.rm = TRUE) ##
  # I create the distribution for the GSADF test
  gsadf <- apply(Msadfs, 1, max, na.rm = TRUE) ##

  #plot density of test stat
  d <- density(as.numeric(gsadf))
  name<-paste("GSADF_","T",T,"_r",round(r0, digits = 3),".png", sep="")
  name <- paste0(address, name)
  plot(d, main ="GSADF", type = "l", col = "orchid")
  dev.copy(png,filename = name)
  dev.off()

  quantile_gsadf <- quantile(gsadf,qe)
  results <-list(quantile_gsadf=quantile_gsadf, gsadf_cv=gsadf_cv)
  return(results)
}




#'  Size of the  sup ADF test
#'
#' This function computes the size of the sup ADF statistic to replicate the results in Table 2 by Phillips et al. (2015)
#'
#' @param qe is a q x 1 vector of quantiles
#' @param m is the number of replications for the simulated DGP
#' @param T is the number of simulated time steps
#' @param r0 is the minimum window size fraction
#' @param lag is the number of lags in the ADF regression
#' @param select is a character object choosing either a "fixed" number of lags in the ADF regression or by using the the "BIC"/"AIC" criteria
#' @param cv is the critical value of the SADF test
#' @param clust_number is the number of clusters for parallel computation
#' @return size is the size of the test
#' @details
#' This function computes the size of the sup ADF statistic to replicate the results in Table 2 by Phillips et al. (2015),
#' @export
#' @importFrom urca ur.df
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel stopCluster
#' @importFrom pbapply pbapply
#'
#' @examples
#'
#'  \dontrun{
#'  # Replicate sizes of the SADF test in Table 2 by PSY(2015)
#'  qe <-c(0.90,0.95,0.99)
#'  m <- 1000           # number of replications
#'
#'  table2_SADF <- matrix(0, 5, 3)
#'  table2_SADF[1,1] <- size_SADF(qe, m, 100, 0.190, 0 , "Fixed", 1.30, clust_number=8)
#'  table2_SADF[2,1] <- size_SADF(qe, m, 200, 0.137, 0 , "Fixed", 1.40, clust_number=8)
#'  table2_SADF[3,1] <- size_SADF(qe, m, 400, 0.100, 0 , "Fixed", 1.49, clust_number=8)
#'  table2_SADF[4,1] <- size_SADF(qe, m, 800, 0.074, 0 , "Fixed", 1.53, clust_number=8)
#'  table2_SADF[5,1] <- size_SADF(qe, m, 1600,0.055, 0 , "Fixed", 1.57, clust_number=8)
#'
#'  table2_SADF[1,2] <- size_SADF(qe, m, 100, 0.190, 3 , "Fixed", 1.30, clust_number=8)
#'  table2_SADF[2,2] <- size_SADF(qe, m, 200, 0.137, 3 , "Fixed", 1.40, clust_number=8)
#'  table2_SADF[3,2] <- size_SADF(qe, m, 400, 0.100, 3 , "Fixed", 1.49, clust_number=8)
#'  table2_SADF[4,2] <- size_SADF(qe, m, 800, 0.074, 3 , "Fixed", 1.53, clust_number=8)
#'  table2_SADF[5,2] <- size_SADF(qe, m, 1600,0.055, 3 , "Fixed", 1.57, clust_number=8)
#'
#'  table2_SADF[1,3] <- size_SADF(qe, m, 100, 0.190, 6 , "BIC", 1.30, clust_number=8)
#'  table2_SADF[2,3] <- size_SADF(qe, m, 200, 0.137, 6 , "BIC", 1.40, clust_number=8)
#'  table2_SADF[3,3] <- size_SADF(qe, m, 400, 0.100, 6 , "BIC", 1.49, clust_number=8)
#'  table2_SADF[4,3] <- size_SADF(qe, m, 800, 0.074, 6 , "BIC", 1.53, clust_number=8)
#'  table2_SADF[5,3] <- size_SADF(qe, m, 1600,0.055, 6 , "BIC", 1.57, clust_number=8)
#'  }

size_SADF <- function(qe, m, T, r0, lag ,select, cv, clust_number) {
  swindow0=floor(r0*T)
  set.seed(1)
  y<-DGP(T,m)
  cl <- parallel::makeCluster(clust_number)
  parallel::clusterEvalQ(cl, library(urca))

  badfs <- pbapply::pbapply(y, 2, function(col_y) {
    sapply(swindow0:T, function(i) {
      x <- urca::ur.df(col_y[1:i], type = "drift", lags = lag, selectlags = select)
      return(x@teststat[1])
    })
  }, cl= cl)
  parallel::stopCluster(cl)
  sadf <- apply(badfs, 2, max)
  size <- sum(sadf > cv)/m
  return(size)
}

#' Size of the generalized sup ADF (GSADF) statistic.
#'
#' This function computes the size of the generalized sup ADF (GSADF) statistic to replicate the results in Table 2 by Phillips et al. (2015)
#'
#' @param qe is a q x 1 vector of quantiles
#' @param m is the number of replications for the simulated DGP
#' @param T is the number of simulated time steps
#' @param r0 is the minimum window size fraction
#' @param lag is the number of lags in the ADF regression
#' @param select is a character object choosing either a "fixed" number of lags in the ADF regression or by using the the "BIC"/"AIC" criteria
#' @param cv is the critical value of the GSADF test
#' @param clust_number is the number of clusters for parallel computation
#' @return size is the size of the test
#' @details
#' This function computes the size of the generalized sup ADF (GSADF) statistic to replicate the results in Table 2 by Phillips et al. (2015),
#' @export
#' @importFrom urca ur.df
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#'
#' @examples
#'
#' \dontrun{
#' # Replicate sizes of the GSADF test in Table 2 by PSY(2015)
#' qe <-c(0.90,0.95,0.99)
#' m <- 1000           # number of replications
#' clust_number <- 7
#'
#' table2_GSADF <- matrix(0, 5, 3)
#' table2_GSADF[1,1] <- size_GSADF(qe, m, 100, 0.190, 0 , "Fixed", 2.00, clust_number)
#' table2_GSADF[2,1] <- size_GSADF(qe, m, 200, 0.137, 0 , "Fixed", 2.08, clust_number)
#' table2_GSADF[3,1] <- size_GSADF(qe, m, 400, 0.100, 0 , "Fixed", 2.20, clust_number)
#' table2_GSADF[4,1] <- size_GSADF(qe, m, 800, 0.074, 0 , "Fixed", 2.34, clust_number)
#' table2_GSADF[5,1] <- size_GSADF(qe, m, 1600,0.055, 0 , "Fixed", 2.41, clust_number)
#'
#' table2_GSADF[1,2] <- size_GSADF(qe, m, 100, 0.190, 3 , "Fixed", 2.00, clust_number)
#' table2_GSADF[2,2] <- size_GSADF(qe, m, 200, 0.137, 3 , "Fixed", 2.08, clust_number)
#' table2_GSADF[3,2] <- size_GSADF(qe, m, 400, 0.100, 3 , "Fixed", 2.20, clust_number)
#' table2_GSADF[4,2] <- size_GSADF(qe, m, 800, 0.074, 3 , "Fixed", 2.34, clust_number)
#' table2_GSADF[5,2] <- size_GSADF(qe, m, 1600,0.055, 3 , "Fixed", 2.41, clust_number)
#'
#' table2_GSADF[1,3] <- size_GSADF(qe, m, 100, 0.190, 6 , "BIC", 2.00, clust_number)
#' table2_GSADF[2,3] <- size_GSADF(qe, m, 200, 0.137, 6 , "BIC", 2.08, clust_number)
#' table2_GSADF[3,3] <- size_GSADF(qe, m, 400, 0.100, 6 , "BIC", 2.20, clust_number)
#' table2_GSADF[4,3] <- size_GSADF(qe, m, 800, 0.074, 6 , "BIC", 2.34, clust_number)
#' table2_GSADF[5,3] <- size_GSADF(qe, m, 1600,0.055, 6 , "BIC", 2.41, clust_number)
#' }

size_GSADF <- function(qe, m, T, r0, lag ,select, cv, clust_number){
  swindow0=floor(r0*T)
  cl <- makeCluster(clust_number)
  y <- DGP(T,m)
  gsadf<- rep(1,m)

  registerDoParallel(cl)

  r <- foreach ( j = 1 : m, .packages = 'urca' ) %dopar%
  {
    #####################################
    # This part of code creates a log.txt file
    # in your current working directory and write
    # the number of current iteration in it
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",j,"\n"))
    sink()
    #####################################
    sadfs <- rep(0,T-swindow0+1)
    for (r2 in swindow0:T)
    {
      dim0 <- r2 - swindow0 + 1
      rwadft <- rep( 0 , dim0)
      for (r1 in 1:dim0)
      {
        x <- ur.df(y[(r1:r2),j],type = "drift",lags = lag, selectlags = select)  #% two tail 5% significant level
        rwadft[r1] <- x@teststat[1]
      }
      sadfs[r2 - swindow0 + 1] <- max(rwadft)
    }
    gsadf[j] <- max(sadfs)
  }
  stopCluster(cl)
  for (j in 1:m) {gsadf[j] <- r[[j]][]}

  size <- sum(gsadf > cv) / m
  return(size)
}



#'  SADF test statistic and the BADF sequence for a real time series data
#'
#' This function computes the SADF test statistic and the BADF sequence for a real time series data
#'
#' @param y is a T x 1 numeric vector
#' @param r0 is the minimum window size fraction
#' @param lag is the number of lags in the ADF regression
#' @param select is a character object choosing either a "Fixed" number of lags in the ADF regression or by using the the "BIC"/"AIC" criteria
#' @param clust_number is the number of clusters for parallel computation
#' @return results is a list containing the SADF test statistic and the BADF sequence
#' @details
#' This function computes the SADF test statistic and BADF sequence for a real time series data. The minimum window size fraction suggested by PSY(2015) is r0=0.01+1.8/sqrt(length(y))
#' @export
#' @importFrom urca ur.df
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @examples
#'
#' \dontrun{
#' #Replicate the empirical analysis in Table 8 by PSY(2015)
#'
#' path.bit <- system.file("extdata", "SP_DV.csv", package = "bubble")
#' SP_DV <- read.table(path.bit,  head = TRUE, sep = ";", fill=T, stringsAsFactors = FALSE)
#' SP_DV <- na.omit( SP_DV[,"PV.ratio", drop=F])
#' r0=0.01+1.8/sqrt(length(SP_DV[,1]))
#' sadf_test <-SADF_Y(SP_DV[,1], r0, 0 ,"Fixed",clust_number=8)
#'
#' qe <-c(0.90,0.95,0.99)    #quantiles
#' m <- 1000
#' T=length(SP_DV[,1])
#' cv_sadf<- CV_SADF(qe,m,T,r0, address)
#' cat(sadf_test$sadf, cv_sadf, sep = "\t")
#' }

SADF_Y <- function(y, r0, lag ,select,clust_number){
  T <- length(y)  # the number of observations
  swindow0=floor(r0*T)
  dim <- T - swindow0 + 1
  badfs <- rep(0,dim)
  cl <- parallel::makeCluster(clust_number)
  parallel::clusterEvalQ(cl, library(urca))
  registerDoParallel(cl)
  r <- foreach ( j = swindow0:T, .packages = 'urca' ) %dopar%
  {
    x <- ur.df(y[1:j], type = "drift", lags = lag, selectlags = select)
    teststat <- x@teststat[1]
  }
  stopCluster(cl)
  for (j in 1:dim) {badfs[j] <- r[[j]][]}
  sadf <- max(badfs)
  results <-list(sadf=sadf, badfs=badfs)
  return(results)
}



#'  GSADF test statistic and the BSADF sequence for a real time series data
#'
#' This function computes the GSADF test statistic and the BSADF sequence for a real time series data
#'
#' @param y is a T x 1 numeric vector
#' @param r0 is the minimum window size fraction
#' @param lag is the number of lags in the ADF regression
#' @param select is a character object choosing either a "Fixed" number of lags in the ADF regression or by using the the "BIC"/"AIC" criteria
#' @param clust_number is the number of clusters for parallel computation
#' @return results is a list containing the GSADF test statistic and the BSADF sequence
#' @details
#' This function computes the GSADF test statistic and BSADF sequence for a real time series data. The minimum window size fraction suggested by PSY(2015) is r0=0.01+1.8/sqrt(length(y))
#' @export
#' @importFrom urca ur.df
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @examples
#'
#' \dontrun{
#' #Replicate the empirical analysis in Table 8 by PSY(2015)
#'
#' # Read file
#' path.bit <- system.file("extdata", "SP_DV.csv", package = "bubble")
#' SP_DV <- read.table(path.bit,  head = TRUE, sep = ";", fill=T, stringsAsFactors = FALSE)
#' SP_DV <- na.omit( SP_DV[,"PV.ratio", drop=F])
#' r0=0.01+1.8/sqrt(length(SP_DV[,1]))
#' gsadf_test <-GSADF_Y(SP_DV[,1], r0, 6 , "BIC", 8)
#'
#' qe <-c(0.90,0.95,0.99)    #quantiles
#' m <- 1000
#' T=length(SP_DV[,1])
#' clust_number=7
#' cv_gsadf<- CV_GSADF(qe,m,T,r0, address)
#' cat(gsadf_test$gsadf, cv_gsadf, sep = "\t")
#' }
GSADF_Y <- function(y, r0, lag ,select,clust_number){
  T <- length(y)  # the number of observations
  swindow0=floor(r0*T)
  dim <- T - swindow0 + 1
  r2 <- c(swindow0:T)
  rw <- r2 - swindow0 + 1
  bsadfs <- rep(0,dim)
  cl <- parallel::makeCluster(clust_number)
  parallel::clusterEvalQ(cl, library(urca))
  registerDoParallel(cl)

  r <- foreach ( j = 1:length(r2), .packages = 'urca' ) %dopar%
  {
    #####################################
    # This creates a log.txt file in your current working directory
    # and write the number of current iteration in it
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",j,"\n"))
    sink()
    #####################################
    swindow <- c(swindow0:r2[j])
    r1 <- r2[j]-swindow+1
    rwadft <- rep(0 , length ( swindow ))
    for (i in 1:length(swindow))
    {
      x <- ur.df(y[r1[i]:r2[j]],lags = lag, type = "drift", selectlags = select)
      rwadft[i] <- x@teststat[1]
    }
    bsadfs_r <- max(rwadft)
  }
  stopCluster(cl)
  for (j in 1:dim) {bsadfs[j] <- r[[j]][]}

  gsadf=max(bsadfs)
  results <-list(gsadf=gsadf, bsadfs=bsadfs)
  return(results)
}
