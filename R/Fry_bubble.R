#' Log-likelihood for a Gaussian random walk model
#'
#' This function computes the log-likelihood for a Gaussian random walk model
#'
#' @param x is a n x 1 numeric data vector
#' @return y a 1 x 1 scalar which is the log-likelihood
#' @details
#' This function computes the log-likelihood for the Gaussian random walk model. It is used as a benchmark in Fry (2014)
#'
#' @export
#' @importFrom graphics plot
#' @importFrom PerformanceAnalytics CalculateReturns
#' @importFrom stats dnorm var
#' @examples
#'
#'  \dontrun{
#'  #Create data set with log-returns
#'  library(PerformanceAnalytics)
#'  path.bit=system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat2 <- subset(dat, as.numeric(Date) > 898 & as.numeric(Date) < 1233)
#'  row.names(dat2) <- levels(dat2$Date)[dat2$Date]
#'  dat3 <- dat2[, 'Close', drop=FALSE]
#'
#'  #compute log-returns
#'  bitret <- CalculateReturns(dat3, method="log")
#'  bitret = bitret[-1,]
#'  plot(bitret)
#'  plot(cumsum(bitret))
#'
#'  #Compute the log-likelihood of the Gaussian random walk model as a benchmark
#'  normlklhd(bitret)
#'  }
#'

normlklhd <- function(x){
  n<-length(x)
  mu<-mean(x)
  mu<-rep(mu, n)
  sigma<-sqrt(((n-1)/n)*var(x))
  sigma<-rep(sigma, n)
  normlklhd<-sum(log(dnorm(x, mu, sigma)))
  normlklhd
  }


#' Log-likelihood for the univariate bubble model by Fry (2014)
#'
#' This function computes the log-likelihood for the univariate bubble model by Fry (2014)
#'
#' @param x is a 4 x 1 parameter vector
#' @param data is a (n-1) x 1 numeric data vector
#' @return y a 1 x 1 scalar which is the log-likelihood
#' @details
#' This function computes the log-likelihood for the univariate bubble model by Fry (2014).
#' The parameters vector consists of the following parameters:
#'
#'   x[1]:      mean
#'
#'   x[2]:      variance (sigma)
#'
#'   x[3]:      alpha (from hazard function)
#'
#'   x[4]:      beta (from hazard function)
#'
#' @export
#' @importFrom graphics plot
#' @importFrom PerformanceAnalytics CalculateReturns
#' @importFrom numDeriv hessian
#' @importFrom stats dnorm
#'
#' @examples
#'
#'  \dontrun{
#'  #Create data set with log-returns
#'  library(PerformanceAnalytics)
#'  library(numDeriv)
#'  path.bit=system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat2 <- subset(dat, as.numeric(Date) > 898 & as.numeric(Date) < 1233)
#'  row.names(dat2) <- levels(dat2$Date)[dat2$Date]
#'  dat3 <- dat2[, 'Close', drop=FALSE]
#'
#'  #compute log-returns
#'  bitret <- CalculateReturns(dat3, method="log")
#'  bitret = bitret[-1,]
#'  plot(bitret)
#'  plot(cumsum(bitret))
#'
#'  #Set (good) starting values and optimize
#'  mu <- mean(bitret)
#'  n <- length(bitret)
#'  sigma <- sqrt(((n-1)/n)*var(bitret))
#'  startx.1 <- c(mu, sigma, 4, 4)
#'  result.1 <- optim(startx.1, jmf1, data=bitret,control=list(maxit=1000, fnscale=-1), hessian = TRUE)
#'
#'  #Compute hessian and variance-covariance matrix
#'   hess <- hessian(x=result.1$par, data=bitret, func=jmf1)
#'   hess2 <- -solve(hess) # variance-covariance matrix
#'   hess2
#'
#'  #Compute v and mu.tilde from estimated parameters and replicate
#'  #the parameter estimates of bubble model reported in Table 2, p.35, by Cheah and Fry (2015)
#'   mu.est <- result.1$par[1]; mu.est;
#'   sigma.est <- result.1$par[2]; sigma.est;
#'   alpha.est <- result.1$par[3]; alpha.est;
#'   beta.est <- result.1$par[4]; beta.est;
#'
#'   v <- round((sqrt(sigma.est*alpha.est*(beta.est-1)^((1/beta.est)-1))),3)
#'   v
#'   mu.tilda <- round((mu.est + (1/2)*sigma.est), 5)
#'   mu.tilda
#'
#'  }
#'

jmf1<-function(x, data){
  x1<-x[1]      #mean
  x2<-x[2]      #sigma
  x4<-x[3]      #alpha
  x5<-x[4]      #beta
  x3<-sqrt(x2*x4*(x5-1)^((1/x5)-1)) #v

  yobs <- data

  #need to make sure this matches the series you are fitting the model to
  n<-length(yobs)+1
  t2<-seq(1, n-1)
  t1<-t2-1

  mu<-x1+x3*log((x4^x5+t2^x5)/(x4^x5+t1^x5))
  sigma<-x2*(1-x4*(x5-1)^((1/x5)-1)*log((x4^x5+t2^x5)/(x4^x5+t1^x5)))
  loglklhd<-sum(log(dnorm(yobs, mu, sqrt(sigma))))
  loglklhd
}


#' Log-likelihood for the univariate bubble model by Fry (2014) - second variant
#'
#' This function computes the second variant of the log-likelihood for the univariate bubble model by Fry (2014),
#' expressed in terms of the parameter "v"
#'
#' @param x is a 4 x 1 parameter vector
#' @param data is a (n-1) x 1 numeric data vector
#' @return y a 1 x 1 scalar which is the log-likelihood
#' @details
#' This function computes the second variant of the log-likelihood for the univariate bubble model by Fry (2014).
#' The parameters vector consists of the following parameters:
#'
#'   x[1]:      mean
#'
#'   x[2]:      v = - ln[(1 - k)] > 0, where k percent is automatically wiped off
#'              the value of the asset in case of market crash
#'
#'   x[3]:      sigma (variance)
#'
#'   x[4]:      beta (from hazard function)
#'
#' @export
#' @importFrom graphics plot
#' @importFrom PerformanceAnalytics CalculateReturns
#' @importFrom numDeriv hessian
#' @importFrom stats dnorm
#'
#' @examples
#'
#'  \dontrun{
#'  #Create data set with log-returns
#'  library(PerformanceAnalytics)
#'  library(numDeriv)
#'  path.bit=system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat2 <- subset(dat, as.numeric(Date) > 898 & as.numeric(Date) < 1233)
#'  row.names(dat2) <- levels(dat2$Date)[dat2$Date]
#'  dat3 <- dat2[, 'Close', drop=FALSE]
#'
#'  #compute log-returns
#'  bitret <- CalculateReturns(dat3, method="log")
#'  bitret = bitret[-1,]
#'  plot(bitret)
#'  plot(cumsum(bitret))
#'
#' #Optimize using good starting values
#'  mu <- mean(bitret)
#' startx.2 <- c(mu, 0.546, 0.007, 1.136)
#' result.2 <- optim(startx.2, jmf1MK, data=bitret, control=list(maxit=1000, fnscale=-1), hessian=TRUE)
#'
#' #Compute Hessian and var/covar matrix
#' hess3 <- hessian(x=result.2$par, data=bitret, func=jmf1MK)
#' hess4 <- -solve(hess3) # variance-covariance matrix
#' hess4
#'
#' #Estimates
#' mu.2.est <- result.2$par[1]
#' v.2.est <- result.2$par[2]
#' sigma.2.est <- result.2$par[3]
#' beta.2.est <- result.2$par[4]
#' alpha.2.est <- (v.2.est^2)/(sigma.2.est*(beta.2.est-1)^((1/beta.2.est)-1))
#'
#' mu.tilda.2 <- round((mu.2.est + (1/2)*sigma.2.est), 5)
#' v.2 <- round(v.2.est,5)
#'
#' #Standard deviations
#'
#' stdev.mu.tilda.2 <- round(sqrt(hess4[1,1] + (1/4)*hess4[3,3] + hess4[1,3]),5)
#' stdev.v.2 <- round(sqrt(hess4[2,2]),3)
#' stdev.mu.tilda.2
#' stdev.v.2
#'
#' #T-stats and p-values
#' t.mu.tilda <- (mu.2.est + (1/2)*sigma.2.est)/(sqrt(hess4[1,1] + (1/4)*hess4[3,3] + hess4[1,3]))
#' t.v <- v.2.est/(sqrt(hess4[2,2]))
#' t.mu.tilda.r <- round(t.mu.tilda,3)
#' t.v.r <- round(t.v,3)
#'
#' pvalue.mu.tilda <- round(2*pt(t.mu.tilda, 329, lower.tail = FALSE),3)
#' pvalue.v <- round(2*pt(t.v, 329, lower.tail = FALSE),3)
#' pvalue.mu.tilda
#' pvalue.v
#'
#' # Reproduce Table 2, p.35, by Cheah and Fry (2015)
#'
#' DF.1 <- data.frame(Parameter=c("nu", "mu.tilda"),
#'                  Estimate=c(v.2, mu.tilda.2),
#'                   E.S.E.=c(stdev.v.2, stdev.mu.tilda.2),
#'                   "t-value"=c(t.v.r, t.mu.tilda.r),
#'                   "p-value"=c(pvalue.v, pvalue.mu.tilda)
#'                   )
#'
#' DF.1
#'
#' }

jmf1MK<-function(x, data){
  x1<-x[1]      #mean
  x3<-x[2]      #v
  x2<-x[3]      #sigma
  x5<-x[4]      #beta

  x4<-(x3^2)/(x2*(x5-1)^((1/x5)-1))      #alpha

  yobs <- data

  #need to make sure this matches the series you are fitting the model to
  n<-length(yobs)+1
  t2<-seq(1, n-1)
  t1<-t2-1

  mu<-x1+x3*log((x4^x5+t2^x5)/(x4^x5+t1^x5))
  sigma<-x2*(1-x4*(x5-1)^((1/x5)-1)*log((x4^x5+t2^x5)/(x4^x5+t1^x5)))
  loglklhd<-sum(log(dnorm(yobs, mu, sqrt(sigma))))
  loglklhd
}



#' Bubble component in percentage
#'
#' This function computes the bubble component (in percentage) using eq.(16) on page 34 of Cheah and Fry (2015)
#'
#' @param alpha.est is a parameter from hazard function
#' @param beta.est is a parameter from hazard function
#' @param v is equal to - ln[(1 - k)] > 0 and it is a key parameter of the Fry (2014) model
#' @param num.obs is the number of time observations
#' @return bubble.component is the bubble component (in percentage)
#' @details
#' This function computes the bubble component (in percentage) using eq.(16) on page 34 of Cheah and Fry (2015)
#'
#' @export
#' @importFrom stats integrate
#'
#' @references Cheah, E.T., and Fry,J. (2015). Speculative bubbles in Bitcoin markets? An empirical investigation into the fundamental value of Bitcoin. Economics Letters, 130, 32-36.
#' @references Fry, J. (2014). Multivariate bubbles and antibubbles. The European Physical Journal B, 87(8), 174.
#'
#' @examples
#' bubble.comp.perc(alpha.est=33.84404, beta.est=1.135961, v=0.546, num.obs=333)
#'

bubble.comp.perc<-function(alpha.est, beta.est, v, num.obs){
  f <- function(t) { (1+ (t^beta.est)/(alpha.est^beta.est))^(-(v-(v^2)/2)) }
  bubble.component<- ( 1-(1/num.obs)*stats::integrate(f, lower = 0, upper = num.obs)$value )*100
  return(bubble.component)
}
