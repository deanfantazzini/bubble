#' Compute the LPPLS Confidence indicator which satisfies the filtering condition 1 in Table 1 by Sornette et al. (2015)
#'
#' This function computes the LPPLS Confidence indicator which satisfies the filtering condition 1 in Table 1 by Sornette et al. (2015)
#'
#' @param x is a T x 1 numeric data vector (or xts object)
#' @param time.scale is the (maximum) rolling window size
#' @param step.shrinking.window is the amount by which the time.scale is decreased at each step
#' @param max.window is the maximum amount to be substracted from time.scale to create shrinking estimation windows. The initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by the amount given by step.shrinking.window  and a new LPPL estimation is performed with window size (t2-t1+step.shrinking.window);
#' this process continues till the estimation window size is (t2-t1+max.window). As a consequence, there are [max.window/step.shrinking.window+1] estimations with endpoint given by t2,
#' whereas the starting dates are given by t1+i, where i=0,step.shrinking.window,2*step.shrinking.window...,max.window.
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
#' - num.osc: the average estimate of the LPPL indicator called 'number of oscillations' (=omega/(2*pi) x log(abs((tc-t1)/(tc-t2))) ) over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - damping: the average estimate of the LPPL indicator called 'Damping' (= (beta x abs(B))/(omega x abs(C)) ) over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T
#'
#' - rel.err: the average estimate of the LPPL indicator called 'relative error' (= abs((Y-Yfit)/Yfit) ) over all rolling estimation window sizes with endpoint t2, where t2=time.scale,time.scale+1, ..., T.
#'
#' - lppl.Confidence: the average fraction of fitting windows for which the LPPLS calibrations satisfy the filtering condition 1 in Table 1 by Sornette et al. (2015). This is called the LPPL Confidence indicator.
#'
#' - col.index: numerical sequential index from 1 to the total number of rows of the final xts object.
#'
#' - crash lock-in: same information as in tc, but in date format.
#'
#' The function also compute a plot with 3 components:
#'
#' - the original price series (with the bubble period highlighted in light red)
#'
#' - the fraction of fitted LPPL satisfying the filtering condition 1 (over time)
#'
#' - the crash lock-in plot, that is the  rolling estimates for the critical time tc
#'
#' @details
#' This function computes the LPPLS Confidence indicator which satisfies the filtering condition 1 in Table 1 by Sornette et al. (2015). More specifically, it estimates the LPPL model using the 3-step procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016) where, for each endpoint t2, the initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by the amount given by step.shrinking.window  and a new LPPL estimation is performed with window size (t2-t1+step.shrinking.window);
#' this process continues till the estimation window size is (t2-t1+max.window). As a consequence, there are [max.window/step.shrinking.window+1] estimations with endpoint given by t2,
#' whereas the starting dates are given by t1+i, where i=0,step.shrinking.window,2*step.shrinking.window...,max.window. All these estimations are then used to compute
#' the filtering conditions 1 in Table 1 by Sornette et al. (2015).
#' The estimated parameters and indicators with the rolling LPPL regression are saved in the xts object lppl.roll.
#' The original price series (with the bubble period highlighted in light red), the fraction of fitted LPPL satisfying the filtering condition 1, and the crash lock-in plot are all reported in a single plot with 3 rows.
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
#'  #Compute the LPPL confidence indicator  with bitcoin prices
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  path.bit <- system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat <- xts(dat[,2], order.by=as.Date(dat[,1]))
#'  a1<-lppl.confidence.A(x=dat,time.scale=750,max.window=600,step.shrinking.window=150,clust_number=8)
#'  a1
#'  }

lppl.confidence.A<- function(x, time.scale, max.window, step.shrinking.window, clust_number)
{
  lppl.roll <- rollapply(as.xts(x), width = time.scale,
                         FUN = function(z) {
                           stepf=seq(0, max.window, by=step.shrinking.window)
                           a1 <- matrix(0,12,length(stepf));
                           cl <- makeCluster(clust_number) #
                           registerDoParallel(cl) #
                           r <- foreach ( i = stepf, .packages=c('zoo', 'tseries', 'bubble')) %dopar%
                           {
                             lppl.par <- lppl_estimate_rob_3all(z[(1+i):length(z)], max.win.tc=0.2)
                             num.osc <-  (lppl.par$param.JLS[2]/(2*pi))*log(abs( lppl.par$param.JLS[7] / (lppl.par$param.JLS[7] - length(z[(1+i):length(z)]))  ))
                             damping <- (lppl.par$param.JLS[1]*abs(lppl.par$param.JLS[5]))/(lppl.par$param.JLS[2]*abs(lppl.par$param.JLS[6]))
                             conditions1 <- 0
                             if (lppl.par$param.JLS[1]>0.01 & lppl.par$param.JLS[1]< 1.2 & lppl.par$param.JLS[2]>2 & lppl.par$param.JLS[2]< 25 & lppl.par$param.JLS[7]>0.95*(length(z[(1+i):length(z)])) &  lppl.par$param.JLS[7]<1.11*(length(z[(1+i):length(z)])) & num.osc>2.5 & damping >0.8 & lppl.par$relative.error < 0.05) {
                               conditions1 <-1
                             }
                             all.data=c(lppl.par$param.JLS, length(z[(1+i):length(z)]), num.osc, damping, lppl.par$relative.error,  conditions1 )

                           }
                           stopCluster(cl) #
                           for (i in 1 : length(stepf)) {a1[,i] <- r[[i]][]} ##
                           rownames(a1)<-c("bet","ome","phi","A","B","C", "tc", "avg.length", "num.osc", "damping", "rel.err", "lppl.Confidence")
                           a1[7,] <- (a1[7,]-a1[8,])+length(z)
                           a2 <- apply(a1,1,mean)
                           return(a2)
                         }, by.column = FALSE, align = "right")

   ### PLOTS
  zx <- as.zoo(x)[, 1]
  zbubble <- as.zoo(lppl.roll[,12])[, 1]
  ## this draws on top using semi-transparent colors.
  rgb <- hcl(c(0, 0, 260), c = c(100, 0, 100), l = c(50, 90, 50), alpha = 0.5)
  par(mar=c(4,3,3,2)) #for Rstudio users -otherwise the plot dimensions are too big
  par(mfcol = c(3,1))
  # 1) plot the price
  plot(zx, log="y", xlab="", main="Price series") #plot in log scale
  xblocks(zbubble > 0.05, col = rgb[1])
  # 2) plot the DS LPPLS Confidence
  plot(coredata(lppl.roll[,12]), type="l", main="DS LPPLS Confidence")

  col.index<-seq(1,nrow(lppl.roll),1)
  #crash dates in xts time format
  crash.lockin<-coredata(lppl.roll[,7]-time.scale)+index(lppl.roll)
  #crash dates in rows index format
  lppl.roll[,7]<-coredata(lppl.roll[,7]-time.scale)+ col.index
  # 3) plot the crash lock-in plot
  plot(index(lppl.roll), crash.lockin, xaxt="n", yaxt="n", type="l", main="Crash lock-in plot (rolling estimates for the critical time tc)")
  axis.Date(1,at=index(lppl.roll), format="%m-%Y")
  axis.Date(2,at=crash.lockin, format="%m-%Y")
  par(mfcol = c(1,1))
  lppl.roll<-cbind.data.frame(lppl.roll, col.index, crash.lockin)
  colnames(lppl.roll)[13]<-"col.index"
  return(lppl.roll)
}



#' Compute the LPPLS Confidence indicator which satisfies the filtering conditions used by Fantazzini (2016)
#'
#' This function computes the LPPLS Confidence indicator which satisfies the filtering conditions used by Fantazzini (2016)
#'
#' @param x is a T x 1 numeric data vector (or xts object)
#' @param time.scale is the (maximum) rolling window size
#' @param step.shrinking.window is the amount by which the time.scale is decreased at each step
#' @param max.window is the maximum amount to be substracted from time.scale to create shrinking estimation windows. The initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by the amount given by step.shrinking.window  and a new LPPL estimation is performed with window size (t2-t1+step.shrinking.window);
#' this process continues till the estimation window size is (t2-t1+max.window). As a consequence, there are [max.window/step.shrinking.window+1] estimations with endpoint given by t2,
#' whereas the starting dates are given by t1+i, where i=0,step.shrinking.window,2*step.shrinking.window...,max.window.
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
#' - confidence.pos.bub: the average fraction of fitting windows for which the LPPLS calibrations satisfy the filtering conditions used by Fantazzini (2016) for a positive bubble.
#'
#' - confidence.neg.bub: the average fraction of fitting windows for which the LPPLS calibrations satisfy the filtering conditions used by Fantazzini (2016) for a negative bubble.
#'
#' - col.index: numerical sequential index from 1 to the total number of rows of the final xts object.
#'
#' - crash lock-in: same information as in tc, but in date format.
#'
#' The function also compute a plot with 3 components:
#'
#' - the original price series (with a positive bubble period highlighted in light red, while a negative bubble in light green)
#'
#' - the fraction of fitted LPPL satisfying the filtering conditions for a positive bubble and a negative bubble (over time)
#'
#' - the crash lock-in plot, that is the rolling estimates for the critical time tc
#'
#' @details
#' This function computes the LPPLS Confidence indicator which satisfies the filtering conditions used by Fantazzini (2016). More specifically, it estimates the LPPL model using the 3-step procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016) where, for each endpoint t2, the initial maximum rolling window size
#' for estimating the LPPL model is time.scale (=t2-t1); then, the window size is decreased by 1 step and a new LPPL estimation is performed with window size (t2-t1+1);
#' this process continues till the estimation window size is (t2-t1+max.window).  As a consequence, there are (max.window+1) estimations for each endpoint t2, which are then used to compute
#' the filtering conditions  used by Fantazzini (2016).
#' The estimated parameters and indicators with the rolling LPPL regression are saved in the xts object lppl.roll.
#' The original price series (with a positive bubble period highlighted in light red, while a negative bubble in light green), the fraction of fitted LPPL satisfying the filtering conditions for a positive bubble and a negative bubble, and the crash lock-in plot are all reported in a single plot with 3 rows.
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
#'
#' @references Fantazzini, D. (2016). The oil price crash in 2014/15: Was there a (negative) financial bubble?. Energy Policy, 96, 383-396.
#' @references Geraskin, P., & Fantazzini, D. (2013). Everything you always wanted to know about log-periodic power laws for bubble modeling but were afraid to ask. The European Journal of Finance, 19(5), 366-391.
#'
#' @examples
#'
#'  \dontrun{
#'  #Compute the LPPL confidence indicator  with bitcoin prices
#'  #load data on bitcoin downloaded from coindesk:  http://www.coindesk.com/price/
#'  path.bit <- system.file("extdata", "coindesk-bpi-USD-close.csv", package = "bubble")
#'  dat <- read.table(path.bit, dec = ".", sep =",", header = TRUE)
#'  dat <- xts(dat[,2], order.by=as.Date(dat[,1]))
#'  b1<-lppl.confidence.B(x=dat,time.scale=750,max.window=600,step.shrinking.window=150,clust_number=8)
#'  b1
#'  }


lppl.confidence.B<- function(x, time.scale, max.window, step.shrinking.window, clust_number)
{
  lppl.roll <- rollapply(as.xts(x), width = time.scale,
                         FUN = function(z) {
                           stepf=seq(0, max.window, by=step.shrinking.window)
                           a1 <- matrix(0,10,length(stepf));
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
                             all.data=c(lppl.par$param.JLS, length(z[(1+i):length(z)]), conditions1.pos.bub, conditions1.neg.bub )

                           }
                           stopCluster(cl) #
                           for (i in 1 : length(stepf)) {a1[,i] <- r[[i]][]} ##
                           rownames(a1)<-c("bet","ome","phi","A","B","C", "tc", "avg.length", "confidence.pos.bub", "confidence.neg.bub")
                           a1[7,] <- (a1[7,]-a1[8,])+length(z)
                           a2=apply(a1,1,mean)
                           return(a2)
                         }, by.column = FALSE, align = "right")

  ### PLOT
  zx <- as.zoo(x)[, 1]
  zpos.bubble <- as.zoo(lppl.roll[,9])[, 1]
  zneg.bubble <- as.zoo(lppl.roll[,10])[, 1]
  ## this draws on top using semi-transparent colors.
  rgb <- hcl(c(0, 0, 260), c = c(100, 0, 100), l = c(50, 90, 50), alpha = 0.5)
  par(mar=c(4,3,3,2)) #for Rstudio users -otherwise the plot dimensions are too big
  par(mfcol = c(3,1))
  # 1) plot the price
  plot(zx, log="y", xlab="", main="Price series") #plot in log scale
  xblocks(zpos.bubble > 0.05, col = rgb[1])
  xblocks(zneg.bubble > 0.05, col = rgb[2])

  # 2) plot the DS LPPLS Confidence
  plot(coredata(lppl.roll[,9]), type="l", main="LPPL Confidence indicator for positive and negative bubbles")
  lines(coredata(lppl.roll[,10]))
  legend("bottomleft", c("Confidence indicator for negative bubbles", "Confidence indicator for positive bubbles"), bg="transparent", lty=1:2)

  col.index<-seq(1,nrow(lppl.roll),1)
  #crash dates in xts time format
  crash.lockin<-coredata(lppl.roll[,7]-time.scale)+index(lppl.roll)
  #crash dates in rows index format
  lppl.roll[,7]<-coredata(lppl.roll[,7]-time.scale)+ col.index
  # 3) plot the crash lock-in plot
  plot(index(lppl.roll), crash.lockin, xaxt="n", yaxt="n", type="l", main="Crash lock-in plot (rolling estimates for the critical time tc)")
  axis.Date(1,at=index(lppl.roll), format="%m-%Y")
  axis.Date(2,at=crash.lockin, format="%m-%Y")
  par(mfcol = c(1,1))
  lppl.roll<-cbind.data.frame(lppl.roll, col.index, crash.lockin)
  return(lppl.roll)
}
