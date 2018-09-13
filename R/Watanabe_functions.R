#' Estimate the particular AR(1) model by Watanabe et al.(2007a,b) using a simplified method
#'
#' This function estimates the particular AR(1) model by Watanabe et al.(2007a,b) using a simplified method
#'
#' @param x is a T x 1 numeric data vector
#' @return vec.1r.1 is a 1 x 1 numeric scalar, which is the estimated omega(i,ti) parameter
#' @details
#' This function estimates the particular AR(1) model by Watanabe et al.(2007a,b) using a simplified method
#' which does not required nonlinear optimization and returns the estimated omega(i,ti) parameter
#' @export
#' @importFrom stats ar.ols
est.war1<- function(x){
  P.null <- as.double(x[1])
  p.t <- x - P.null
  ar.1 <- stats::ar.ols(p.t, aic = FALSE, order.max = 1, intercept = FALSE)
  vec.1r.1 <- as.double(ar.1[2])
  return(vec.1r.1)
}

#' Estimate the particular AR(1) model by Watanabe et al.(2007a,b) using the original method
#'
#' This function estimates the particular AR(1) model by Watanabe et al.(2007a,b) using the original method
#'
#' @param x is a T x 1 numeric data vector
#' @return vec.1r.1 is a 1 x 1 numeric scalar, which is the estimated omega(i,ti) parameter
#' @details
#' This function estimates the particular AR(1) model by Watanabe et al.(2007a,b) using the original method
#' and returns the estimated omega(i,ti) parameter
#' @export
#' @importFrom stats nls coef lag ar.ols
est.war1.original<- function(x){
  P.null.start <- as.double(x[1])
  #Compute a good starting value for omega(i,ti)
  p.t <- x - P.null.start
  ar.1 <- stats::ar.ols(p.t, aic = FALSE, order.max = 1, intercept = FALSE)
  omega.start <- as.double(ar.1[2])

  dat<-as.data.frame(cbind(p.diff=diff(x), p.lag=lag(x,k=-1)))
  colnames(dat)<-c("p.diff", "p.lag")
  dat<-dat[-nrow(dat),]
  out <- tryCatch(
    {
      war1<-stats::nls(p.diff ~ (omega-1)*(p.lag-p0), data=dat,  start = list(omega = omega.start, p0=P.null.start))
      vec.1r.1 <- as.double(coef(war1)[1])
    },
    error = function(e){warning("NUmerical optimization failed"); -9999}
  )
  return(out)
}



#' Simulate data from a AR(p) model and estimate the particular AR(1) model by Watanabe et al.(2007a,b) using a simplified method
#'
#' This function simulates data from a AR(p) model and estimate the particular AR(1) model by Watanabe et al.(2007a,b) using a simplified method
#'
#' @param ar.par is a a p x 1 vector containing the AR parameters for the simulation of the AR(p)
#' @param step.win is the number of simulated time steps
#' @return vec.1r.1 is a 1 x 1 numeric scalar, which is the estimated omega(i,ti) parameter
#' @details
#' This function simulates data from a AR(p) model and estimate the particular AR(1) model by Watanabe et al.(2007a,b) using a simplified
#' method which does not required nonlinear optimization,
#'  and returns the estimated omega(i,ti) parameter
#' @export
#' @importFrom stats ar.ols
#' @importFrom stats arima.sim
est.war1.from.arp <- function(ar.par, step.win){
  Price.sim <- arima.sim(model=list(ar.par), n=step.win)
  P.null <- as.double(Price.sim[1])
  p.t <- Price.sim - P.null
  ar.1 <- stats::ar.ols(p.t, aic = FALSE, order.max = 1, intercept = FALSE)
  vec.1r.1 <- as.double(ar.1[2])
  return(vec.1r.1)
}

#' Simulate data from a AR(p) model and estimate the particular AR(1) model by Watanabe et al.(2007a,b) using the original method
#'
#' This function simulates data from a AR(p) model and estimate the particular AR(1) model by Watanabe et al.(2007a,b) using the original method
#'
#' @param ar.par is a a p x 1 vector containing the AR parameters for the simulation of the AR(p)
#' @param step.win is the number of simulated time steps
#' @return vec.1r.1 is a 1 x 1 numeric scalar, which is the estimated omega(i,ti) parameter
#' @details
#' This function simulates data from a AR(p) model and estimate the particular AR(1) model by Watanabe et al.(2007a,b) using the original method,
#'  and returns the estimated omega(i,ti) parameter
#' @export
#' @importFrom stats nls coef lag
#' @importFrom stats arima.sim ar.ols
est.war1.from.arp.original <- function(ar.par, step.win){
  Price.sim <- arima.sim(model=list(ar.par), n=step.win)
  P.null.start <- as.double(Price.sim[1])
  #Compute a good starting value for omega(i,ti)
  p.t <- Price.sim - P.null.start
  ar.1 <- stats::ar.ols(p.t, aic = FALSE, order.max = 1, intercept = FALSE)
  omega.start <- as.double(ar.1[2])

  dat<-as.data.frame(cbind(p.diff=diff(Price.sim), p.lag=lag(Price.sim,k=-1)))
  colnames(dat)<-c("p.diff", "p.lag")
  dat<-dat[-nrow(dat),]
  out <- tryCatch(
         {
           war1<-stats::nls(p.diff ~ (omega-1)*(p.lag-p0), data=dat,  start = list(omega = omega.start, p0=P.null.start))
           vec.1r.1 <- as.double(coef(war1)[1])
         },
         error = function(e){warning("NUmerical optimization failed"); -9999}
         )
  return(out)
}


#' Best time-scale according to the approach by Watanabe et al. (2007a,b)
#'
#' This function looks for the best time-scale according to the approach by Watanabe et al. (2007a,b)
#'
#' @param xdat is a T x 1 numeric data vector
#' @param maxlag is the maximum lag to be used to select the optimal AR lag order with AIC
#' @param min.win is the minimum window size to be tried during the simulations
#' @param max.win is the maximum window size to be tried during the simulations
#' @param n.replications is the number of simulated AR(p) processes (for each window size)
#' @param clust_number is the number of clusters for parallel computation
#' @param original if TRUE the original method by Watanabe et al. (2007a,b) is used, otherwise
#' a simplified approach which does not required nonlinear optimization.
#' @details
#' This function looks for the best time-scale according to the approach by Watanabe et al. (2007a,b), using either the original method
#' or a simplified approach which does not required nonlinear optimization.
#' The optimal windown size according to Watanabe et al. (2007a,b) is the minimum window size for which always omega(i, Ti)<=1 holds.
#' If an estimated parameter omega(i,ti)>1, the following message is printed to the screen,
#'  "series is divergent, at least for one iteration = [", i, "],  omega(i,Ti) > 1"
#'  and the function stops checking further.
#' @export
#' @importFrom stats ar.ols
#' @importFrom pbapply pbreplicate
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel stopCluster
#' @examples
#'
#'  \dontrun{
#'  x= arima.sim(list(order = c(3,0,0), ar = c(1.2, -0.2, -0.1)), n = 500)
#'  optimal.time.scale(xdat=x,maxlag=10,min.win=5,max.win=20,n.replications=2000,
#'  clust_number=8,original=FALSE)
#'  }

optimal.time.scale <- function(xdat, maxlag=10, min.win=10, max.win=20, n.replications=2000, clust_number=8, original=TRUE){
  # Select AR(p) model (without constant) and with AIC
  ar.X <- ar.ols(xdat, aic = TRUE, order.max = maxlag, intercept = FALSE)
  arpar<-as.vector(ar.X[2])
  step.win <- min.win
  diff.win <- max.win - min.win
  for (i in 1:diff.win){
    step.win <- min.win + i
    message("Window size = ", step.win, " ")
    cl <- parallel::makeCluster(clust_number)
    parallel::clusterEvalQ(cl, library(bubble))
    parallel::clusterEvalQ(cl, c("arpar", "step.win", "original"))
    if (original==TRUE){
      output <- as.vector(pbapply::pbreplicate(n.replications,est.war1.from.arp.original(arpar, step.win), cl=cl))
    }else{
      output <- as.vector(pbapply::pbreplicate(n.replications,est.war1.from.arp(arpar, step.win), cl=cl))
    }
    parallel::stopCluster(cl)
    for (i in 1:length(output)){
      if (output[i] > 1){
        message(" series is divergent, at least for one iteration = [", i, "],  omega(i,Ti) > 1" )
        break()
      }else if ( output[i] == -9999) {
        message(" Numerical optimization failure at iteration = [", i, "]" )
        break()
      }
    }
  }
}

#' Estimate the particular Watanabe et al.(2007a,b) model for detecting financial bubbles
#'
#' This function estimates the particular Watanabe et al.(2007a,b) model for detecting financial bubbles using a rolling regression with window size equal to time.scale and data given by x
#'
#' @param x is a T x 1 numeric data vector (or xts object)
#' @param time.scale is the optimal window size determined using the optimal.time.scale function
#' @param original if TRUE the original method by Watanabe et al. (2007a,b) is used, otherwise
#' a simplified approach which does not required nonlinear optimization.
#' @return results is a list containing the rolling estimated AR(1) parameters and the bubble index
#' @details
#' This function estimates the particular Watanabe et al.(2007a,b) model for detecting financial bubbles using a rolling AR(1) regression with window size equal to time.scale and data given by x
#' The estimated parameters with the rolling AR(1) regression are saved in the xts object war1.par, while a bubble index which is 1 in case of a bubble according to the method by Watanabe et al.(2007a,b) in that specific date and zero otherwise is saved in the xts object bubble.index .
#' Both the rolling estimated AR(1) parameters and the bubble index are plotted, with the bubble period highlighted in light red on the original xts object/numeric vector.
#' @export
#' @importFrom zoo rollapply
#' @importFrom xts xts
#' @importFrom xts as.xts
#' @importFrom zoo as.zoo
#' @importFrom zoo index
#' @importFrom zoo xblocks
#' @importFrom grDevices hcl
#' @importFrom graphics par
#' @importFrom graphics abline
#' @examples
#'
#'  \dontrun{
#'  #Example 1: Replicate the analysis of plot (a) in Fig.1 in Watanabe et al.(2007b)
#'  library(quantmod)
#'  R <- getSymbols("AABA", src = "yahoo", from = as.Date("1998-01-01"), to = as.Date("2001-12-31"))
#'  head(AABA)
#'  yahoo.adj <- AABA[,"AABA.Adjusted", drop = F]
#'  Watanabe.bubble.index(yahoo.adj, 100)
#'
#'  ### Example 2: Compute the EXponential Curve Fitting (EXCF) method  with bitcoin prices
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  path.bit <- system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat <- xts(dat[,2], order.by=as.Date(dat[,1]))
#'  Watanabe.bubble.index(dat, 100)
#'  }

Watanabe.bubble.index<- function(x, time.scale, original= FALSE){
  if (original==TRUE){
  war1.par <- rollapply(as.xts(x), width = time.scale,
                        FUN = est.war1.original,
                        by.column = FALSE, align = "right")
  }else{
  war1.par <- rollapply(as.xts(x), width = time.scale,
                          FUN = est.war1,
                          by.column = FALSE, align = "right")
  }

  bubble.index=xts(rep(0,length(war1.par)), index(war1.par))
  for (i in time.scale:length(war1.par))
  {
    if (war1.par[i] > 1) {
      bubble.index[(i-time.scale+1):i] <-1
    } else {
      bubble.index[i]<-0
    }
  }

  z <- as.zoo(war1.par)[, 1]
  zx <- as.zoo(x)[, 1]
  zbubble <- as.zoo(bubble.index)[, 1]
  ## this draws on top using semi-transparent colors.
  rgb <- hcl(c(0, 0, 260), c = c(100, 0, 100), l = c(50, 90, 50), alpha = 0.5)
  par(mar=c(4,3,3,2)) #for Rstudio users -otherwise the plot dimensions are too big

  par(mfcol = c(2,1))
  plot(zx, log="y", xlab="", ylab=colnames(x)) #plot in log scale
  xblocks(zbubble > 0, col = rgb[1])
  plot(z, main="Rolling estimates of omega(i,Ti)", xlab="", ylab="omega(i,Ti)")
  abline(h =1)
  par(mfcol = c(1,1))
  results <-list(bubble.index=bubble.index, war1.par=war1.par)
  return(results)
}
