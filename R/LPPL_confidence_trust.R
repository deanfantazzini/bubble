#' Compute the LPPLS Confidence and Trust indicators which satisfy the filtering condition 1 and 2, respectively, in Table 1 by Sornette et al. (2015)
#'
#' This function computes the LPPLS Confidence and Trust indicators which satisfy the filtering condition 1 and 2, respectively,  in Table 1 by Sornette et al. (2015)
#'
#' @param x is a T x 1 numeric data vector (or xts object)
#' @param time.scale is the (maximum) rolling window size
#' @param step.shrinking.window is the amount by which the time.scale is decreased at each step
#' @param max.window is the maximum amount to be substracted from time.scale to create shrinking estimation windows. The initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by the amount given by step.shrinking.window  and a new LPPL estimation is performed with window size (t2-t1+step.shrinking.window);
#' this process continues till the estimation window size is (t2-t1+max.window). As a consequence, there are [max.window/step.shrinking.window+1] estimations with endpoint given by t2,
#' whereas the starting dates are given by t1+i, where i=0,step.shrinking.window,2*step.shrinking.window...,max.window.
#' @param boot.rep is the number of bootstrp replications used to compute the LPPL trust indicator at each ending point t2
#' @param clust_number is the number of clusters for parallel computation
#' @return lppl.roll is a xts object containing the following columns:
#'
#' - bet: the average estimate of the LPPL parameter beta over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - ome: the average estimate of the LPPL parameter omega over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - phi: the average estimate of the LPPL parameter phi over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - A: the average estimate of the LPPL parameter A over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - B: the average estimate of the LPPL parameter B over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - C: the average estimate of the LPPL parameter C over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - tc: the average estimate of the LPPL parameter tc (= the critical time) over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T. It is reported as number of rows.
#'
#' - length: the average length of the estimation window over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - num.osc: the average estimate of the LPPL indicator called 'number of oscillations' (=omega/(2*pi) x log(abs((tc-t1)/(tc-t2))) )over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - damping: the average estimate of the LPPL indicator called 'Damping' (= (beta x abs(B))/(omega x abs(C)) ) over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - rel.err: the average estimate of the LPPL indicator called 'relative error' (= abs((Y-Yfit)/Yfit) ) over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T.
#'
#' - lppl.Confidence: the average fraction of fitting windows for which the LPPLS calibrations satisfy the filtering condition 1 in Table 1 by Sornette et al. (2015). This is called the LPPL Confidence indicator.
#'
#' - lppl.trust: the median level over the [max.window/step.shrinking.window+1] estimations windows of the fraction among the boot.rep synthetic time series that satisfy the filtering condition 2 in Table 1 by Sornette et al. (2015).
#'   The synthetic time series are created by resampling boot.rep times the LPPL residuals and adding them to the calibrated LPPLS structure.
#'
#' - col.index: numerical sequential index from 1 to the total number of rows of the final xts object.
#'
#' - crash lock-in: same information as in tc, but in date format.
#'
#' The function also compute a plot with 4 components:
#'
#' - the original price series (with the bubble period highlighted in light red)
#'
#' - the fraction of fitted LPPL satisfying the filtering condition 1 (over time), that is the DS LPPL Confidence indicator.
#'
#' - the DS LPPL Trust indicator.
#'
#' - the crash lock-in plot, that is the rolling estimates for the critical time tc
#'
#' @details
#' This function computes the LPPLS Confidence and Trust indicators which satisfy the filtering condition 1 and 2, respectively, in Table 1 by Sornette et al. (2015). More specifically, it estimates the LPPL model using the 3-step procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016) where, for each endpoint t2, the initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by the amount given by step.shrinking.window and a new LPPL estimation is performed with window size (t2-t1+step.shrinking.window);
#' this process continues till the estimation window size is (t2-t1+max.window).  As a consequence, there are [max.window/step.shrinking.window+1] estimations for each endpoint t2, which are then used to compute
#' the filtering conditions 1 in Table 1 by Sornette et al. (2015), that is the DS LPPL Confidence indicator. For each endpoint t2, the DS LPPL Trust indicator is also computed:
#' this is the median level over the [max.window/step.shrinking.window+1] estimations windows of the fraction among the boot.rep synthetic time series that satisfy the filtering condition 2 in Table 1 by Sornette et al. (2015).
#' The synthetic time series are created by resampling boot.rep times the LPPL residuals and adding them to the calibrated LPPLS structure.
#'
#' The estimated parameters and the DS LPPL Confidence and Trust indicators computed with the rolling LPPL regression are saved in the xts object lppl.roll.
#' The original price series (with the bubble period highlighted in light red), the DS LPPL Confidence and Trust indicators, and the crash lock-in plot are all reported in a single plot with 4 rows.
#' @export
#' @importFrom zoo rollapply
#' @importFrom xts xts
#' @importFrom xts as.xts
#' @importFrom zoo as.zoo
#' @importFrom zoo index
#' @importFrom zoo xblocks
#' @importFrom zoo as.yearmon
#' @importFrom grDevices hcl
#' @importFrom graphics par
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @importFrom stats median
#' @importFrom graphics axis.Date
#'
#' @references Fantazzini, D. (2016). The oil price crash in 2014/15: Was there a (negative) financial bubble?. Energy Policy, 96, 383-396.
#' @references Geraskin, P., & Fantazzini, D. (2013). Everything you always wanted to know about log-periodic power laws for bubble modeling but were afraid to ask. The European Journal of Finance, 19(5), 366-391.
#' @references Sornette, D., Demos, G., Qun, Z., Cauwels, P.,  Zhang, Q. (2015). Real-time prediction and post-mortem analysis of the Shanghai 2015 stock market bubble and crash. Journal of Investment Strategies, 4(4), 77-95.
#' @references Zhang, Q., Sornette, D., Balcilar, M., Gupta, R., Ozdemir, Z. A., Yetkiner, H. (2016). LPPLS bubble indicators over two centuries of the SP500 index. Physica A: Statistical Mechanics and its Applications, 458, 126-139.
#'
#' @examples
#'
#'  \dontrun{
#'  # Example1: Compute the LPPL confidence indicator  with the SP500 prices
#'  from 1791 till 2015 (adjusted for dividends and splits), using a 269-month time scale
#'  data(sp500)
#'  ct=lppl.confidence.trust.A(x=sp500,time.scale=269,max.window=120,
#'  step.shrinking.window=10,boot.rep=100,clust_number=7)
#'  ct
#'
#'  ## Example 2:  replicate the example in Sornette et al. (2015)
#'  library(quantmod)
#'  SSE<-getSymbols("000001.SS",auto.assign=FALSE,from=as.Date("2012-07-01"), to=as.Date("2015-07-31"))
#'  sse.adj <- SSE[,"000001.SS.Adjusted", drop = F]
#'  head(sse.adj)
#'  ct=lppl.confidence.trust.A(x=sse.adj,time.scale=250,max.window=120,
#'  step.shrinking.window=10,boot.rep=10,clust_number=7)
#'  ct
#'
#'  ### Example 3: Compute the LPPL trust indicator  with bitcoin prices
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  path.bit <- system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat <- xts(dat[,2], order.by=as.Date(dat[,1]))
#'  a1<-lppl.confidence.trust.A(x=dat,time.scale=750,max.window=600,step.shrinking.window=150,
#'  boot.rep=10,clust_number=8)
#'  a1
#'  }


lppl.confidence.trust.A<- function(x, time.scale, max.window, step.shrinking.window, boot.rep, clust_number)
{
  lppl.roll <- rollapply(as.xts(x), width = time.scale,
                FUN = function(z) {
                  boot.rep=boot.rep
                  stepf=seq(0, max.window, by=step.shrinking.window)
                  a1 <- matrix(0,13,length(stepf))
                  cl <- makeCluster(clust_number) #
                  registerDoParallel(cl) #
                  r <- foreach ( i = stepf, .packages=c('zoo', 'tseries', 'bubble')) %dopar%
                  {
                      lppl.par <- lppl_estimate_rob_3all(z[(1+i):length(z)], max.win.tc=0.2)
                      num.osc <-  (lppl.par$param.JLS[2]/(2*pi))*log(abs( lppl.par$param.JLS[7] / (lppl.par$param.JLS[7] - length(z[(1+i):length(z)]))  ))
                      damping <- (lppl.par$param.JLS[1]*abs(lppl.par$param.JLS[5]))/(lppl.par$param.JLS[2]*abs(lppl.par$param.JLS[6]))
                      conditions1 <- 0
                      if (lppl.par$param.JLS[1]>0.01 & lppl.par$param.JLS[1]< 1.2 & lppl.par$param.JLS[2]>2 & lppl.par$param.JLS[2]< 25 & lppl.par$param.JLS[7]>0.95*(length(z[(1+i):length(z)])) &  lppl.par$param.JLS[7]<1.1*(length(z[(1+i):length(z)])) & num.osc >2.5 & damping >0.8 & lppl.par$relative.error < 0.05) {
                        conditions1 <-1
                      }

                      a.boot <- matrix(0,1,boot.rep)
                      for(j in 1:boot.rep)
                      {
                        z.est<- lppl_simulate_boot(length(lppl.par$resids), lppl.par$param.JLS, lppl.par$resids)
                        z.boot.est <- lppl_estimate_rob_3all(z.est)
                        num.osc <- (z.boot.est$param.JLS[2]/(2*pi))*log(abs( z.boot.est$param.JLS[7] / (lppl.par$param.JLS[7] - length(z.est))  ))
                        damping <- (z.boot.est$param.JLS[1]*abs(z.boot.est$param.JLS[5]))/(z.boot.est$param.JLS[2]*abs(z.boot.est$param.JLS[6]))
                        conditions2.boot <- 0
                        if (z.boot.est$param.JLS[1]>0.01 & z.boot.est$param.JLS[1]< 0.99 & z.boot.est$param.JLS[2]>2 & z.boot.est$param.JLS[2]< 25 & z.boot.est$param.JLS[7]>0.95*(length(z[(1+i):length(z)])) &  z.boot.est$param.JLS[7]<1.1*(length(z[(1+i):length(z)])) & num.osc >2.5 & damping >1 & z.boot.est$relative.error < 0.2) {
                          conditions2.boot <-1
                        }
                        a.boot[,j] <- conditions2.boot
                      }
                      rownames(a.boot)<-c("cond2")
                      a2.boot=apply(a.boot,1,mean)
                      all.data=c(lppl.par$param.JLS, length(z[(1+i):length(z)]), num.osc, damping, lppl.par$relative.error,  conditions1, a2.boot )
                    }
                    stopCluster(cl) #
                    for (i in 1 : length(stepf)) {a1[,i] <- r[[i]][]} ##
                    rownames(a1)<-c("bet","ome","phi","A","B","C", "tc", "avg.length", "num.osc", "damping", "rel.err", "lppl.Confidence", "lppl.Trust")
                    a1[7,] <- (a1[7,]-a1[8,])+length(z)
                    a2a=apply(a1[1:12,],1,mean)
                    a2b=median(a1[13,]); names(a2b)<-"lppl.Trust"
                    a2=c(a2a,a2b)
                    return(a2)
                  }, by.column = FALSE, align = "right")

  ## PLOT
  zx <- as.zoo(x)[, 1]
  zbubble1 <- as.zoo(lppl.roll[,12])[, 1]
  zbubble2 <- as.zoo(lppl.roll[,13])[, 1]
  ## this draws on top using semi-transparent colors.
  rgb <- hcl(c(0, 0, 260), c = c(100, 0, 100), l = c(50, 90, 50), alpha = 0.5)
  par(mar=c(4,3,3,2)) #for Rstudio users -otherwise the plot dimensions are too big
  par(mfcol = c(4,1))
  # 1) plot the price
  plot(zx, log="y", xlab="", main="Price series (in log scale)")
  xblocks(zbubble1 > 0.05, col = rgb[1])
  # 2) plot the DS LPPLS Confidence
  plot(coredata(lppl.roll[,12]), type="l", main="DS LPPLS Confidence")
  # 3) plot the DS LPPLS Trust
  plot(coredata(lppl.roll[,13]), type="l", main="DS LPPL Trust")

  col.index<-seq(1,nrow(lppl.roll),1)
  #crash dates in xts time format
  crash.lockin<-coredata(lppl.roll[,7]-time.scale)+index(lppl.roll)
  #crash dates in rows index format
  lppl.roll[,7]<-coredata(lppl.roll[,7]-time.scale)+ col.index
  # 4) plot the crash lock-in plot
  plot(index(lppl.roll), crash.lockin, xaxt="n", yaxt="n", type="l", main="Crash lock-in plot (rolling estimates for the critical time tc)")
  axis.Date(1,at=index(lppl.roll), format="%m-%Y")
  axis.Date(2,at=crash.lockin, format="%m-%Y")
  par(mfcol = c(1,1))
  lppl.roll<-cbind.data.frame(lppl.roll, col.index, crash.lockin)
  return(lppl.roll)
}



#' Compute the LPPLS Confidence and Trust indicators  which satisfy the filtering conditions used by Fantazzini (2016)
#'
#' This function the LPPLS Confidence and Trust indicators  which satisfy the filtering conditions used by Fantazzini (2016)
#'
#' @param x is a T x 1 numeric data vector (or xts object)
#' @param time.scale is the (maximum) rolling window size
#' @param step.shrinking.window is the amount by which the time.scale is decreased at each step
#' @param max.window is the maximum amount to be substracted from time.scale to create shrinking estimation windows. The initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by the amount given by step.shrinking.window  and a new LPPL estimation is performed with window size (t2-t1+step.shrinking.window);
#' this process continues till the estimation window size is (t2-t1+max.window). As a consequence, there are [max.window/step.shrinking.window+1] estimations with endpoint given by t2,
#' whereas the starting dates are given by t1+i, where i=0,step.shrinking.window,2*step.shrinking.window...,max.window.
#' @param boot.rep is the number of bootstrp replications used to compute the LPPL trust indicator at each ending point t2
#' @param clust_number is the number of clusters for parallel computation
#' @return lppl.roll is a xts object containing the following columns:
#'
#' - bet: the average estimate of the LPPL parameter beta over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - ome: the average estimate of the LPPL parameter omega over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - phi: the average estimate of the LPPL parameter phi over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - A: the average estimate of the LPPL parameter A over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - B: the average estimate of the LPPL parameter B over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - C: the average estimate of the LPPL parameter C over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - tc: the average estimate of the LPPL parameter tc (= the critical time) over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T. It is reported as number of rows.
#'
#' - length: the average length of the estimation window over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - confidence.pos.bub: the average fraction of fitting windows for which the LPPLS calibrations satisfy the filtering conditions used by Fantazzini (2016) for a positive bubble. This is called the LPPL Confidence indicator for positive bubbles.
#'
#' - confidence.neg.bub: the average fraction of fitting windows for which the LPPLS calibrations satisfy the filtering conditions used by Fantazzini (2016) for a negative bubble. This is called the LPPL Confidence indicator for negative bubbles.
#'
#' - Trust.pos.bub: the median level over the [max.window/step.shrinking.window+1] estimations windows of the fraction among the boot.rep synthetic time series that satisfy the filtering conditions used by Fantazzini (2016) for a positive bubble.
#'   The synthetic time series are created by resampling boot.rep times the LPPL residuals and adding them to the calibrated LPPLS structure.
#'
#' - Trust.neg.bub: the median level over the [max.window/step.shrinking.window+1] estimations windows of the fraction among the boot.rep synthetic time series that satisfy the filtering conditions used by Fantazzini (2016) for a negative bubble.
#'   The synthetic time series are created by resampling boot.rep times the LPPL residuals and adding them to the calibrated LPPLS structure.
#'
#' - col.index: numerical sequential index from 1 to the total number of rows of the final xts object.
#'
#' - crash lock-in: same information as in tc, but in date format.
#'
#' The function also compute a plot with 4 components:
#'
#' - the original price series (with a positive bubble period highlighted in light red, while a negative bubble in light green)
#'
#' - the fraction of fitted LPPL satisfying the filtering conditions for a positive bubble and a negative bubble (over time), that is the DS LPPL Confidence indicator.
#'
#' - the DS LPPL Trust indicator.
#'
#' - the crash lock-in plot, that is the rolling estimates for the critical time tc
#'
#' @details
#' This function computes the LPPLS Confidence indicator which satisfies the filtering conditions used by Fantazzini (2016). More specifically, it estimates the LPPL model using the 3-step procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016) where, for each endpoint t2, the initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by 1 step and a new LPPL estimation is performed with window size (t2-t1+1);
#' this process continues till the estimation window size is (t2-t1+max.window).  As a consequence, there are (max.window+1) estimations for each endpoint t2, which are then used to compute
#' the filtering conditions  used by Fantazzini (2016), which is comparable to the  DS LPPL Confidence indicator.
#' For each endpoint t2, the DS LPPL Trust indicator is also computed:
#' this is the median level over the [max.window/step.shrinking.window+1] estimations windows of the fraction among the boot.rep synthetic time series that satisfy the filtering conditions used by Fantazzini (2016).
#' The synthetic time series are created by resampling boot.rep times the LPPL residuals and adding them to the calibrated LPPLS structure.
#' The estimated parameters and indicators with the rolling LPPL regression are saved in the xts object lppl.roll.
#' The original price series (with a positive bubble period highlighted in light red, while a negative bubble in light green), the DS LPPL Confidence and Trust indicators for a positive bubble and a negative bubble, and the crash lock-in plot are all reported in a single plot with 4 rows.
#' The filtering conditions used by Fantazzini (2016) are reported below:
#'
#' - Positive bubble: 0<beta<1,  B<0, b=[- B x beta - |C| x sqrt(beta^2+omega^2)]>0 (hazard rate), LPPL residuals stationary at the 5\% level (using the KPSS test statistic)
#'
#' - Negative bubble: 0<beta<1,  B>0, b=[- B x beta - |C| x sqrt(beta^2+omega^2)]<0 (hazard rate), LPPL residuals stationary at the 5\% level (using the KPSS test statistic)
#'
#' @export
#' @importFrom zoo rollapply
#' @importFrom xts xts
#' @importFrom xts as.xts
#' @importFrom zoo as.zoo
#' @importFrom zoo index
#' @importFrom zoo xblocks
#' @importFrom zoo as.yearmon
#' @importFrom grDevices hcl
#' @importFrom graphics par
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @importFrom stats median
#'
#' @references Fantazzini, D. (2016). The oil price crash in 2014/15: Was there a (negative) financial bubble?. Energy Policy, 96, 383-396.
#' @references Geraskin, P., & Fantazzini, D. (2013). Everything you always wanted to know about log-periodic power laws for bubble modeling but were afraid to ask. The European Journal of Finance, 19(5), 366-391.
#'
#' @examples
#'
#'  \dontrun{
#'  # Example 1: Compute the LPPL confidence indicator  with the SP500 prices
#'  from 1791 till 2015 (adjusted for dividends and splits), using a 269-month time scale
#'  data(sp500)
#'  ctb=lppl.confidence.trust.B(x=sp500,time.scale=269,max.window=120,
#'  step.shrinking.window=10,boot.rep=100,clust_number=7)
#'  ctb
#'
#'  ## Example 2: replicate the example in Sornette et al. (2015) but using the criteria proposed
#'  # by Fantazzini (2016) and Geraskin and Fantazzini (2013)
#'  library(quantmod)
#'  SSE<-getSymbols("000001.SS",auto.assign=FALSE,from=as.Date("2012-07-01"), to=as.Date("2015-07-31"))
#'  sse.adj <- SSE[,"000001.SS.Adjusted", drop = F]
#'  head(sse.adj)
#'  ctb=lppl.confidence.trust.A(x=sse.adj,time.scale=250,max.window=120,
#'  step.shrinking.window=10,boot.rep=10,clust_number=7)
#'  ctb
#'
#'  ### Example 3: Compute the LPPL trust indicator  with bitcoin prices
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  path.bit <- system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat <- xts(dat[,2], order.by=as.Date(dat[,1]))
#'  b1<-lppl.confidence.trust.B(x=dat,time.scale=750,max.window=600,step.shrinking.window=150,
#'  boot.rep=10,clust_number=8)
#'  b1
#'  }

lppl.confidence.trust.B<- function(x, time.scale, max.window, step.shrinking.window,  boot.rep, clust_number)
{

  lppl.roll <- rollapply(as.xts(x), width = time.scale,
                         FUN = function(z) {
                           boot.rep=boot.rep
                           stepf=seq(0, max.window, by=step.shrinking.window)
                           a1 <- matrix(0,12,length(stepf));
                           cl <- makeCluster(clust_number) #
                           registerDoParallel(cl) #
                           r <- foreach ( i = stepf, .packages=c('zoo', 'tseries', 'bubble')) %dopar%
                           {
                             lppl.par <- lppl_estimate_rob_3all(z[(1+i):length(z)], max.win.tc=0.2)
                             conditions1.pos.bub <- 0
                             conditions1.neg.bub <- 0
                             if (lppl.par$param.JLS[1]>0 & lppl.par$param.JLS[1]< 1 & lppl.par$param.JLS[5]<0 & lppl.par$crash.rate > 0 & lppl.par$KPSS.res <0.463) {
                               conditions1.pos.bub <- 1
                             }
                             if (lppl.par$param.JLS[1]>0 & lppl.par$param.JLS[1]< 1 & lppl.par$param.JLS[5]>0 & lppl.par$crash.rate < 0 & lppl.par$KPSS.res <0.463) {
                               conditions1.neg.bub <- 1
                             }

                             a.pos.boot <- matrix(0,1,boot.rep)
                             a.neg.boot <- matrix(0,1,boot.rep)
                             for(j in 1:boot.rep)
                             {
                               z.est<- lppl_simulate_boot(length(lppl.par$resids), lppl.par$param.JLS, lppl.par$resids)
                               z.boot.est <- lppl_estimate_rob_3all(z.est)
                               conditions2.boot.pos <- 0
                               conditions2.boot.neg <- 0
                               if (z.boot.est$param.JLS[1]>0 & z.boot.est$param.JLS[1]< 1 & z.boot.est$param.JLS[5]<0 & z.boot.est$crash.rate > 0 & z.boot.est$KPSS.res <0.463) {
                                 conditions2.boot.pos <- 1
                               }
                               if (z.boot.est$param.JLS[1]>0 & z.boot.est$param.JLS[1]< 1 & z.boot.est$param.JLS[5]>0 & z.boot.est$crash.rate < 0 & z.boot.est$KPSS.res <0.463) {
                                 conditions2.boot.neg <- 1
                               }
                               a.pos.boot[,j] <- conditions2.boot.pos
                               a.neg.boot[,j] <- conditions2.boot.neg
                             }
                             rownames(a.pos.boot)<-c("cond2")
                             rownames(a.neg.boot)<-c("cond2")
                             a2.pos.boot=apply(a.pos.boot,1,mean)
                             a2.neg.boot=apply(a.neg.boot,1,mean)

                             all.data=c(lppl.par$param.JLS, length(z[(1+i):length(z)]), conditions1.pos.bub, conditions1.neg.bub, a2.pos.boot, a2.neg.boot )

                           }
                           stopCluster(cl) #
                           for (i in 1 : length(stepf)) {a1[,i] <- r[[i]][]} ##
                           rownames(a1)<-c("bet","ome","phi","A","B","C", "tc", "avg.length", "confidence.pos.bub", "confidence.neg.bub", "Trust.pos.bub", "Trust.neg.bub")
                           a1[7,] <- (a1[7,]-a1[8,])+length(z)
                           a2a=apply(a1[1:10,],1,mean)
                           a2b=apply(a1[11:12,],1,median)
                           a2=c(a2a,a2b)
                           return(a2)
                         }, by.column = FALSE, align = "right")

  ### PLOT
  zx <- as.zoo(x)[, 1]
  zpos.bubble <- as.zoo(lppl.roll[,9])[, 1]
  zneg.bubble <- as.zoo(lppl.roll[,10])[, 1]
  ## this draws on top using semi-transparent colors.
  rgb <- hcl(c(0, 0, 260), c = c(100, 0, 100), l = c(50, 90, 50), alpha = 0.5)
  par(mar=c(4,3,3,2)) #for Rstudio users -otherwise the plot dimensions are too big
  par(mfcol = c(4,1))
  # 1) plot the price
  plot(zx, log="y", xlab="", main="Price series (in log scale)") #plot in log scale
  xblocks(zpos.bubble > 0.05, col = rgb[1])
  xblocks(zneg.bubble > 0.05, col = rgb[2])
  # 2) plot the DS LPPLS Confidence
  plot(coredata(lppl.roll[,9]), type="l", main="LPPL Confidence indicator for positive and negative bubbles")
  lines(coredata(lppl.roll[,10]))
  legend("bottomleft", c("Confidence indicator for negative bubbles", "Confidence indicator for positive bubbles"), bg="transparent", lty=1:2)
  # 3) plot the DS LPPLS Trust
  plot(coredata(lppl.roll[,11]), type="l", main="LPPL Trust indicator for positive and negative bubbles")
  lines(coredata(lppl.roll[,12]))
  legend("bottomleft", c("Trust indicator for negative bubbles", "Trust indicator for positive bubbles"), bg="transparent", lty=1:2)

  col.index<-seq(1,nrow(lppl.roll),1)
  #crash dates in xts time format
  crash.lockin<-coredata(lppl.roll[,7]-time.scale)+index(lppl.roll)
  #crash dates in rows index format
  lppl.roll[,7]<-coredata(lppl.roll[,7]-time.scale)+ col.index
  # 4) plot the crash lock-in plot
  plot(index(lppl.roll), crash.lockin, xaxt="n", yaxt="n", type="l", main="Crash lock-in plot (rolling estimates for the critical time tc)")
  axis.Date(1,at=index(lppl.roll), format="%m-%Y")
  axis.Date(2,at=crash.lockin, format="%m-%Y")
  par(mfcol = c(1,1))
  lppl.roll<-cbind.data.frame(lppl.roll, col.index, crash.lockin)

  return(lppl.roll)
}
