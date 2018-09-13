#' Simulate a LPPL model
#'
#' This function simulates a LPPL model
#'
#' @param T is the number of simulated time steps
#' @param true_parm is a a 8 x 1 vector containing the parameters for the simulation
#' @return y a T  x 1 vector of simulated data
#' @details
#' This function simulates a LPPL model with parameter vector given by true_param and make a plot of it. More specifically, the LPPL parameters in the true_param vector are given below:
#'
#' true_parm[1] = beta
#'
#' true_parm[2] = omega
#'
#' true_parm[3] = phi
#'
#' true_parm[4] = A
#'
#' true_parm[5] = B
#'
#' true_parm[6] = C
#'
#' true_parm[7] = sigma (i.e. the error term variance)
#'
#' true_parm[8] = the critical time tc
#'
#' @export
#' @importFrom graphics plot
#' @importFrom stats rnorm
#'
#' @examples
#'
#'  \dontrun{
#'  tparm=c(0.353689, 9.154368, 2.074608, 7.166421,-0.434324, 0.035405, 0.000071, 530)
#'  aa=lppl_simulate(500,tparm)
#'  }
#'

lppl_simulate=function(T=500, true_parm){
    bet=true_parm[1]; ome=true_parm[2]; phi=true_parm[3];
    A= true_parm[4]; B =true_parm[5]; C= true_parm[6]; ws=true_parm[7];
    tc=true_parm[8];
    tt_sim=seq(1, T, 1);
    sdum=rep(1,T);
    f_t=(tc - tt_sim)^bet;
    g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) + phi );
    x=exp(A*sdum +B*f_t + C*g_t +sqrt(ws)*rnorm(T) );
    plot(x, type="l", xlab = "Time index", ylab = "Price")
    return(x)
}


#' Estimation of the LPPL model using a nonlinear optimization
#'
#' This function estimates a LPPL model using a nonlinear optimization
#'
#' @param x is a T x 1 numeric data vector
#' @return par_est is a 7 x 1 vector of estimated parameters
#' @details
#' This function estimates the LPPL model by Johansen, Ledoit, and Sornette (2000) using the original (two-step) nonlinear optimization, see section 4.1 in Geraskin and Fantazzini (2013) for a compact review.
#' The returned parameter vector contained the following parameters:
#'
#' par_est[1] = beta
#'
#' par_est[2] = omega
#'
#' par_est[3] = phi
#'
#' par_est[4] = tc (i.e. the critical time)
#'
#' par_est[5] = A
#'
#' par_est[6] = B
#'
#' par_est[7] = C
#'
#' We remark that this estimation method is not recommended, due to the frequent presence of many
#' local minima of the cost function where the minimization algorithm can get trapped.
#' It was included in this package for historical (and teaching) reasons.
#'
#' @export
#' @importFrom zoo coredata
#' @importFrom stats nlminb
#'
#' @examples
#'
#'  \dontrun{
#'  tparm=c(0.353689, 9.154368, 2.074608, 7.166421,-0.434324, 0.035405, 0.000071, 530)
#'  aa=lppl_simulate(500,tparm)
#'
#'  bb=lppl_estimate(aa); bb;
#'  }
#'

lppl_estimate=function(x){

  # Initialize Model Parameters and Bounds:
  par_start = c(bet=0.5, ome=9, phi=pi, tc=NROW(x) +15)
  lowerBounds = c(bet=0.000001, ome=1, phi=0.00001, tc=NROW(x)+1)
  upperBounds = c(bet=2, ome=21, phi=2*pi, tc=NROW(x)+10000)

  # Compose function to be MINIMIZED:
  LPPL_LLH = function(parm){
    bet = parm[1]; ome = parm[2]; phi = parm[3]; tc =parm[4];
    n=NROW(x);       # number of observations
    sdum=rep(1,n);
    tt_sim=seq(1, n, 1);
    f_t=(tc - tt_sim)^bet;
    g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) + phi );
    X=cbind(sdum,f_t,g_t);
    Y=log(x)
    b_sol=solve(t(X)%*%X)%*%(t(X)%*%Y);
    et=Y-sdum*b_sol[1]- b_sol[2]*f_t - b_sol[3]*g_t;
    et=coredata(et)
    ll=sum(et*et);
    ll
  }

  print(LPPL_LLH(par_start))

  # Estimate Parameters and Compute Numerically Hessian:
  fit = nlminb(par_start, objective=LPPL_LLH, lower = lowerBounds, upper = upperBounds, control = list(trace=3), hessian=TRUE)

  bet = fit$par[1]; ome = fit$par[2]; phi = fit$par[3]; tc = fit$par[4];
  n=NROW(x);       # number of observations
  sdum=rep(1,n);
  tt_sim=seq(1, n, 1);
  f_t=(tc - tt_sim)^bet;
  g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) + phi );
  X=cbind(sdum,f_t,g_t);
  Y=log(x)
  b_sol=solve(t(X)%*%X)%*%(t(X)%*%Y)
  names(b_sol)= c("A","B", "C")
  par_est=c( (fit$par), (b_sol) )
  return(par_est)
}


#' Compute the residuals given the parameters estimated with the function lppl_estimate
#'
#' This function computes the residuals given the parameters estimated with the function lppl_estimate
#'
#' @param x is a T x 1 numeric data vector
#' @param par_est is a 7 x 1 vector of parameters estimated with the function lppl_estimate
#' in this order: ["bet" "ome" "phi" "tc"  "A"   "B"   "C"]
#' @return et is a vector of residuals
#' @details
#' This function computes the residuals given the parameters estimated with the function lppl_estimate
#'
#' @export
#' @importFrom zoo coredata
#'
#' @examples
#'
#'  \dontrun{
#'  tparm <- c(0.353689, 9.154368, 2.074608, 7.166421,-0.434324, 0.035405, 0.000071, 30)
#'  aa <- lppl_simulate(500,tparm)
#'  bb <- lppl_estimate(aa)
#'  bb
#'  resids <- lppl.resids.from.basic.estimate(aa,bb)
#'  }
#'

lppl.resids.from.basic.estimate=function(x,par_est){

  bet = par_est[1]; ome = par_est[2]; phi = par_est[3]; tc = par_est[4];
  A = par_est[5]; B = par_est[6]; C = par_est[7];
  n=NROW(x)
  sdum=rep(1,n)
  tt_sim=seq(1, n, 1)
  f_t=(tc - tt_sim)^bet
  g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) + phi )
  X=cbind(sdum,f_t,g_t)
  Y=log(x)
  et=Y - sdum*A - B*f_t - C*g_t
  et=coredata(et)

  return(et)
}


#' Simulate a LPPL model using bootstrap resampling of residuals
#'
#' This function simulates a LPPL model using bootstrap resampling of residuals
#'
#' @param T is the number of simulated time steps
#' @param true_parm is a a 8 x 1 vector containing the parameters for the simulation
#' @param resids is a T x 1 vector of residuals from a previous LPPL estimation
#' @return y a T x 1 vector of simulated data
#' @details
#' This function simulates a LPPL model with parameter vector given by true_param,
#' and error terms given by the bootstrp resampling of residuals from a previous LPPL estimation.
#' The simulated data are then plotted. More specifically,
#' the LPPL parameters in the true_param vector are given below:
#'
#' true_parm[1] = beta
#'
#' true_parm[2] = omega
#'
#' true_parm[3] = phi
#'
#' true_parm[4] = A
#'
#' true_parm[5] = B
#'
#' true_parm[6] = C
#'
#' true_parm[7] = the critical time tc
#'
#' @export
#' @importFrom graphics plot
#' @importFrom zoo coredata
#' @examples
#'
#'  \dontrun{
#'  tparm=c(0.353689, 9.154368, 2.074608, 7.166421,-0.434324, 0.035405, 0.000071, 530)
#'  aa=lppl_simulate(500,tparm)
#'  bb=lppl_estimate(aa); bb;
#'  resids <- lppl.resids.from.basic.estimate(aa,bb)
#'  tparm.boot=c(0.353689, 9.154368, 2.074608, 7.166421,-0.434324, 0.035405, 530)
#'  aa=lppl_simulate_boot(500,tparm.boot, resids)
#'  }
#'

lppl_simulate_boot=function(T=500, true_parm, resids){
  bet=true_parm[1]; ome=true_parm[2]; phi=true_parm[3];
  A= true_parm[4]; B =true_parm[5]; C= true_parm[6];
  tc=true_parm[7]
  tt_sim=seq(1, T, 1)
  sdum=rep(1,T)
  f_t=(tc - tt_sim)^bet
  g_t=( (tc - tt_sim)^bet )*cos( ome*log(tc - tt_sim) + phi )
  x=exp(A*sdum +B*f_t + C*g_t +sample(coredata(resids),replace = TRUE))
  plot(x, type="l", xlab = "Time index", ylab = "Price")
  return(x)
}

