#' 1st step LPPL estimation procedure  by Geraskin and Fantazzini (2013) and Fantazzini(2016)
#'
#' This function performs the 1st step of the LPPL estimation procedure  by Geraskin and Fantazzini (2013) and Fantazzini (2016)
#'
#' @param x is a T x 1 data vector
#' @param par_start is a 6 x 1 vector of starting values for the parameters to be estimated
#' @param max.win.tc is a scalar setting the max window size (in percentage terms) used to fix the critical time tc in this estimation step
#' @return par_est1 is a 6 x 1 vector of estimated parameters
#' @details
#' This function performs the 1st step of the LPPL estimation procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016) using the LPPL formula
#' by Filimonov and Sornette (2013):
#' The critical time tc is fixed at tc = t2 + 0.1 x (t2-t1), where t2 and t1 are
#' the last and the first observation of the
#' estimation sample, respectively. The remaining LPPL parameters [A, B, C1, C2, beta, omega]
#' are estimated by using a quasi-Newton method algorithm.
#' The starting values for the parameter vector is given by:
#'
#' par_start[1] = beta
#'
#' par_start[2] = omega
#'
#' par_start[3] = A
#'
#' par_start[4] = B
#'
#' par_start[5] = C1
#'
#' par_start[6] = C2
#'
#' @export
#' @importFrom zoo coredata
#' @importFrom stats nlminb

lppl_estimate_rob_1s=function(x, par_start = c(bet=0.5, ome=6, A=log(x[NROW(x)]), B=0, C1=0, C2=0), max.win.tc=0.1){

	# Compose function to be MINIMIZED:
	LPPL_LLH1 = function(parm){
		bet = parm[1]; ome = parm[2]; A = parm[3]; B = parm[4]; C1=parm[5]; C2=parm[6];
		n=NROW(x)       # number of observations
		sdum=rep(1,n)
		tt_sim=seq(1, n, 1)
		tc=NROW(x)+max.win.tc*NROW(x)
		f_t=(tc - tt_sim)^bet;
		g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) )
		h_t=( (tc - tt_sim)^bet )*sin( ome*log(tc - tt_sim) )
		Y=log(x)
		et=Y-sdum*A- B*f_t - C1*g_t - C2*h_t
		et=coredata(et)
		ll=sum(et*et)
 	}

	# Estimate Parameters:
	 fit = nlminb(par_start, objective=LPPL_LLH1, control = list(trace=0), hessian=TRUE)

	par_est1=c( (fit$par) )
	return(par_est1)
}


#' 2nd step LPPL estimation procedure  by Geraskin and Fantazzini (2013) and Fantazzini(2016)
#'
#' This function performs the 2nd step of the LPPL estimation procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016)
#'
#' @param x is a T x 1 data vector
#' @param par is a 6 x 1 vector containing the parameters estimated in the 1st step [beta, omega, A, B, C1, C2] and which are kept fixed in the 2nd step
#' @param max.win.tc is a scalar (in percentage terms) used to set the starting value for the critical time tc
#' @return par_est2 is a 1 x 1 scalar containing the estimated parameter tc
#' @details
#' This function performs the 2nd step of the LPPL estimation procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016) using the LPPL formula
#'  by Filimonov and Sornette (2013):
#' Keeping fixed the LPPL parameters  [beta, omega, A, B, C1, C2] computed in the
#' first stage, the critical time tc is estimated in a second step by using a
#'  quasi-Newton method algorithm.
#'
#' @export
#' @importFrom zoo coredata
#' @importFrom stats nlminb

lppl_estimate_rob_2s=function(x,par, max.win.tc=0.1){

	#  Create the starting value for the critical time tc:
	par_start = c(tc=NROW(x)+max.win.tc*NROW(x))

	#  Compose function to be MINIMIZED:
	LPPL_LLH2 = function(parm){
		tc=parm[1]
		bet = par[1]; ome = par[2]; A = par[3]; B = par[4]; C1=par[5]; C2=par[6];
		n=NROW(x)       # number of observations
		sdum=rep(1,n)
		tt_sim=seq(1, n, 1)
		f_t=(tc - tt_sim)^bet
		g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) )
		h_t=( (tc - tt_sim)^bet )*sin( ome*log(tc - tt_sim) )
		Y=log(x)
		et=Y-sdum*A- B*f_t - C1*g_t - C2*h_t
		et=coredata(et)
		ll=sum(et*et)
 	}

	#Estimate Parameters:
  fit = nlminb(par_start, objective=LPPL_LLH2, control = list(trace=0), hessian=TRUE)

	par_est2=c( (fit$par) )
	return(par_est2)
}


#' 3rd step LPPL estimation procedure  by Geraskin and Fantazzini (2013) and Fantazzini (2016)
#'
#' This function performs the 3rd step of the LPPL estimation procedure  by Geraskin and Fantazzini (2013) and Fantazzini (2016)
#'
#' @param x is a T x 1 data vector
#' @param par1 is a 6 x 1 vector containing the parameters estimated in the 1st step [beta, omega, A, B, C1, C2]
#' @param par2 is a 1 x 1 scalar containing the parameter estimated in the 2nd step [tc]
#' @return par_est is a 8 x 1 vector containing the estimated 7 LPPL parameters [beta, omega, A, B, C1, C2, tc], togother with the KPSS test statistic computed with the LPPL residuals to check their stationarity (a new condition introduced by Lin et al. (2014))
#' @details
#' This function performs the 3rd step of the LPPL estimation procedure
#' by Geraskin and Fantazzini (2013) and Fantazzini (2016) using the LPPL
#' formula by Filimonov and Sornette (2013):
#' Using  the estimated parameters in  the first and second stages as starting
#' values, it estimates all the 7 LPPL parameters [beta, omega, A, B, C1, C2, tc]
#' by using a quasi-Newton method algorithm. Moreover, it computes the KPSS test
#' statistic with the LPPL residuals to check their stationarity (a new condition introduced by Lin et al. (2014)).
#'
#' @export
#' @importFrom zoo coredata
#' @importFrom stats nlminb
#' @importFrom tseries kpss.test
#'
#' @examples
#'
#'  \dontrun{
#'  tparm=c(0.35, 4.15, 2.07, 7.16,-0.43, 0.035, 0.00007, 530)
#'  aa=lppl_simulate(500,tparm)
#'  # 1st estimation step
#'  bb1=lppl_estimate_rob_1s(aa)
#'  bb1
#'  # 2nd estimation step
#'  bb2=lppl_estimate_rob_2s(aa,bb1)
#'  bb2
#'  # 3rd estimation step
#'  bb3=lppl_estimate_rob_3s(aa,bb1,bb2)
#'  bb3
#'  # The original C parameter can be retrieved using standard trigonometric
#'  # functions
#'  C.param=bb3[5]/cos(atan(bb3[6]/bb3[5]))
#'
#'  # The first major condition for a bubble to occur within the LPPL framework
#'  # is that 0 < beta < 1, which guarantees that the crash hazard rate accelerates.
#'  # The second major condition is that the crash rate should be non-negative,
#'  # as  highlighted by  Bothmer and Meister (2003), which imposes that
#'  # b = [- B x beta - |C| x sqrt(beta^2  + omega^2)] >=0.
#'  # For the original parameters we have:
#'  -tparm[5]*tparm[1]-abs(tparm[6])*sqrt(tparm[1]^2+tparm[2]^2)
#'  # while for the estimated paramaters we have:
#'  -bb3[4]*bb3[1]-abs(C.param)*sqrt(bb3[1]^2+bb3[2]^2)
#'    }

lppl_estimate_rob_3s=function(x,par1,par2){

	# Create the starting values for the parameter vector:
	par_start = c(par1,par2)

	# Compose function to be MINIMIZED:
	LPPL_LLH3 = function(parm){
		bet = parm[1]; ome = parm[2]; A = parm[3]; B = parm[4]; C1=parm[5]; C2=parm[6]; tc =parm[7];
		n=NROW(x)      # number of observations
		sdum=rep(1,n)
		tt_sim=seq(1, n, 1)
		f_t=(tc - tt_sim)^bet
		g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) )
		h_t=( (tc - tt_sim)^bet )*sin( ome*log(tc - tt_sim) )
		Y=log(x)
		et=Y-sdum*A- B*f_t - C1*g_t - C2*h_t
		et=coredata(et)
		ll=sum(et*et)
 	}

	# Estimate Parameters:
	  fit = nlminb(par_start, objective=LPPL_LLH3, control = list(trace=0), hessian=TRUE)

	#Computation of residuals and KPSS test statistic
	bet = fit$par[1]; ome = fit$par[2]; A = fit$par[3]; B = fit$par[4]; C1=fit$par[5]; C2=fit$par[6]; tc =fit$par[7];
		n=NROW(x)
		sdum=rep(1,n)
		tt_sim=seq(1, n, 1)
		f_t=(tc - tt_sim)^bet
		g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) )
		h_t=( (tc - tt_sim)^bet )*sin( ome*log(tc - tt_sim) )
		Y=log(x)
		et=Y-sdum*A- B*f_t - C1*g_t - C2*h_t
		kpsst=kpss.test(et)

	par_est=c( (fit$par), kpsst$statistic )
	return(par_est)
}

#' 3-step LPPL estimation procedure by Geraskin and Fantazzini (2013) and Fantazzini (2016)
#'
#' This function performs all the 3 steps of the LPPL estimation procedure  by Geraskin and Fantazzini (2013) and Fantazzini (2016)
#'
#' @param x is a T x 1 data vector
#' @param par_start is a 6 x 1 vector of starting values for the parameters to be estimated in the first step
#' @param max.win.tc is a scalar setting the max window size (in percentage terms) used to fix the critical time tc in the 1st step, and to set the starting value for tc in the 2nd step
#' @return results is a list containing the following objects:
#'
#' - a vector named param.JLS containing the estimated parameters according to the original LPPL structure by Johansen et al. (2001): [beta, omega, phi, A, B, C, tc]
#'
#' - a vector named param.FS containing the estimated parameters according to the LPPL structure by Filimonov and Sornette (2013): [beta, omega, A, B, C1, C2, tc]
#'
#' - the (scalar) KPSS test statistic computed with the LPPL residuals to check their stationarity
#'
#' - the crash hazard rate computed according to Bothmer and Meister (2003)
#'
#' - the Relative Error of the model fit: mean( (Y-Yfit)/Yfit )
#'
#' - the model residuals
#'
#' @details
#'
#' This function performs all the 3 steps of the LPPL estimation procedure  by Geraskin and Fantazzini (2013) and Fantazzini (2016) using the LPPL
#' formula by Filimonov and Sornette (2013). See the functions for each of these estimation steps for more details.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'  tparm=c(0.35, 4.15, 2.07, 7.16,-0.43, 0.035, 0.00007, 530)
#'  aa=lppl_simulate(500,tparm)
#'  est.all<-lppl_estimate_rob_3all(aa)
#'  est.all
#' }

lppl_estimate_rob_3all=function(x, par_start = c(bet=0.5, ome=6, A=log(x[NROW(x)]), B=0, C1=0, C2=0), max.win.tc=0.1){
  bb1=lppl_estimate_rob_1s(x, par_start, max.win.tc)
  bb2=lppl_estimate_rob_2s(x,bb1, max.win.tc)
  bb3=lppl_estimate_rob_3s(x,bb1,bb2)
  C.param=bb3[5]/cos(atan(bb3[6]/bb3[5]))
  crash.rate=-bb3[4]*bb3[1]-abs(C.param)*sqrt(bb3[1]^2+bb3[2]^2)
  names(crash.rate)<-"crash.rate"
  param.JLS=c(bb3[1:2],acos(bb3[5]/C.param), bb3[3:4],C.param, bb3[7])
  names(param.JLS)<-c("bet","ome","phi","A","B","C", "tc")
  param.FS=c(bb3[1:7])
  names(param.FS)<-c("bet","ome","A","B","C1","C2","tc")

  bet = bb3[1]; ome = bb3[2]; A = bb3[3]; B = bb3[4]; C1=bb3[5]; C2=bb3[6]; tc =bb3[7];
  n=NROW(x)
  sdum=rep(1,n)
  tt_sim=seq(1, n, 1)
  f_t=(tc - tt_sim)^bet
  g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) )
  h_t=( (tc - tt_sim)^bet )*sin( ome*log(tc - tt_sim) )
  Y=log(x)
  et=Y-sdum*A- B*f_t - C1*g_t - C2*h_t
  Yfit = sdum*A + B*f_t + C1*g_t + C2*h_t
  relative.error=abs(mean( (Y-Yfit)/Yfit ))

  results <-list(param.JLS=param.JLS, param.FS=param.FS, KPSS.res=bb3[8], crash.rate=crash.rate, relative.error=relative.error, resids=et)
  return(results)
}

#' Simulate a LPPL model using the formula in Filimonov and Sornette (2013)
#'
#' This function simulates a LPPL model using the formula in Filimonov and Sornette (2013)
#'
#' @param T is the number of simulated time steps
#' @param true_parm is a a 8 x 1 vector containing the parameters for the simulation
#' @return y a T  x 1 vector of simulated data
#' @details
#' This function simulates a LPPL model using the formula in Filimonov and Sornette (2013),
#'  with parameter vector given by true_param and make a plot of it. More specifically,
#'  the LPPL parameters in the true_param vector are given below:
#'
#' true_parm[1] = beta
#'
#' true_parm[2] = omega
#'
#' true_parm[3] = A
#'
#' true_parm[4] = B
#'
#' true_parm[5] = C1
#'
#' true_parm[6] = C2
#'
#' true_parm[7] = the critical time tc
#'
#' true_parm[8] = sigma (i.e. the error term variance)
#'
#' @export
#' @importFrom graphics plot
#' @importFrom stats rnorm
#'
#' @examples
#'
#'  \dontrun{
#'  tparm=c(0.35, 4.15, 7.16,-0.43, -0.016, 0.031, 530, 0.00007)
#'  a1=lppl_simulate_FS(500,tparm)
#'  est.all<-lppl_estimate_rob_3all(a1)
#'  est.all$param.FS
#'  }
#'

lppl_simulate_FS=function(T=500, true_parm){
  bet=true_parm[1]; ome=true_parm[2]; A= true_parm[3]; B =true_parm[4];
  C1= true_parm[5]; C2= true_parm[6];
  tc=true_parm[7]; ws=true_parm[8];
  tt_sim=seq(1, T, 1);
  sdum=rep(1,T);
  f_t=(tc - tt_sim)^bet;
  g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) )
  h_t=( (tc - tt_sim)^bet )*sin( ome*log(tc - tt_sim) )
  x=exp(A*sdum +B*f_t + C1*g_t + C2*h_t +sqrt(ws)*rnorm(T) );
  plot(x, type="l")
  return(x)
}

#' Simulate a LPPL model using bootstrap resampling of residuals and the formula in Filimonov and Sornette (2013)
#'
#' This function simulates a LPPL model using bootstrap resampling of residuals and the formula in Filimonov and Sornette (2013)
#'
#' @param T is the number of simulated time steps
#' @param true_parm is a a 7 x 1 vector containing the parameters for the simulation
#' @param resids is a T x 1 vector of residuals from a previous LPPL estimation
#' @return y a T x 1 vector of simulated data
#' @details
#' This function simulates a LPPL model using the formula in Filimonov and Sornette (2013),
#' error terms given by the bootstrp resampling of residuals from a previous LPPL estimation,
#'  and with a parameter vector given by true_param. The simulated data are then plotted. More specifically,
#'  the LPPL parameters in the true_param vector are given below:
#'
#' true_parm[1] = beta
#'
#' true_parm[2] = omega
#'
#' true_parm[3] = A
#'
#' true_parm[4] = B
#'
#' true_parm[5] = C1
#'
#' true_parm[6] = C2
#'
#' true_parm[7] = the critical time tc
#'
#' @export
#' @importFrom graphics plot
#' @importFrom zoo coredata
#'
#' @examples
#'
#'  \dontrun{
#'  tparm=c(0.35, 4.15, 7.16,-0.43, -0.016, 0.031, 530, 0.00007)
#'  a1=lppl_simulate_FS(500,tparm)
#'  est.all<-lppl_estimate_rob_3all(a1)
#'  resids=est.all$resids
#'  tparm.boot=c(0.35, 4.15, 7.16,-0.43, -0.016, 0.031, 530)
#'  aa=lppl_simulate_FS_boot(500,tparm.boot, resids)
#'  }
#'
lppl_simulate_FS_boot=function(T=500, true_parm, resids){
  bet=true_parm[1]; ome=true_parm[2]; A= true_parm[3]; B =true_parm[4];
  C1= true_parm[5]; C2= true_parm[6];
  tc=true_parm[7];
  tt_sim=seq(1, T, 1)
  sdum=rep(1,T)
  f_t=(tc - tt_sim)^bet
  g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) )
  h_t=( (tc - tt_sim)^bet )*sin( ome*log(tc - tt_sim) )
  x=exp(A*sdum +B*f_t + C1*g_t + C2*h_t +sample(coredata(resids),replace = TRUE) )
  plot(x, type="l")
  return(x)
}

