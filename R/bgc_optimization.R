# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr



source("R/constants.R")
source("R/bgc_SS_constant.R")
source("R/bgc_SS_hotspot.R")


#' @title AIC for least square
#'
#' @description Function that gives the AIC of a least-square model
#' The expression if given with the correction for small sample size according to:
#' Banks, H. T., and M. L. Joyner. 2017. AIC under the framework of least squares estimation. Applied Mathematics Letters 74:33–45.
#'
#' @param n: number of observations
#' @param np: number of parameters
#' @param SSres: residual sum of squares
#'
#' @return AIC
#' @export
#'
AICls <- function(n,np,SSres){
  n*log(SSres/n) + 2*(np + 1)*(n + 2)/(n - np)
}


# Functions to minimize the sum of squares of the different models


#' @title Sum of squares minimization of model M
#'
#' @description Function that searches for the three parameters that minimize the sum_of_squares function
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: GC content
#' @param Mmin: minimum for the range of M, default value = -5
#' @param Mmax: maximum or the range of M, default value = 5
#' @param Maxit: maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm: number of updates in the method (option for optim, see manual)
#' @param Verbose: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr: to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @return A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_M(sfsWS,sfsSW,0.5)
#' LS$param$mutbias # mutation bias
#' LS$criteria$AIC # model AIC


least_square_M <- function(WS,SW,GC,
                           Mmin=MMIN,Mmax=MMAX,
                           Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  n <- length(WS)
  Minit <- mean(log(WS/SW)) + log(1 - GC) - log(GC)
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  e1max <- min(SW/(SW + rev(WS))) * ONE
  e2max <- min(WS/(WS + rev(SW))) * ONE
  # The factor ONE is to avoid 0 values in the SFS
  init <- c(Minit,e1max/2,e2max/2)
  # Boundaries for optimization
  inf <- c(Mmin,0,0)
  sup <- c(Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_M else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_M,
    gr = gradient,
    WS = WS, SW = SW, GC = GC,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(n,length(init),SSres)
  return( list(
    "param"=list("mutbias"=exp(minSSE$par[1]),
                 "e1"=minSSE$par[2],
                 "e2"=minSSE$par[3]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}



#' @title Sum of squares minimization of model B
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: GC content
#' @param Bmin: minimum for the range of B, default value = -100
#' @param Bmax: maximum for the range of B, default value = 100
#' @param Maxit: maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm: number of updates in the method (option for optim, see manual)
#' @param Verbose: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr: to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @return A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_B(sfsWS,sfsSW,0.5)
#' LS$param$B # population-scaled gBGC
#' LS$criteria$R2 # R2 of the model



least_square_B <- function(WS,SW,GC,
                           Bmin=BMIN,Bmax=BMAX,
                           Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  n <- length(WS)
  x <- c(1:n)/(n+1)
  reginit <- lm(log(WS/SW) ~ x - 1)
  Binit <- reginit$coef
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  e1max <- min(SW/(SW + rev(WS))) * ONE
  e2max <- min(WS/(WS + rev(SW))) * ONE
  # The factor ONE is to avoid 0 values in the SFS
  init <- c(Binit,e1max/2,e2max/2)
  # Boundaries for optimization
  inf <- c(Bmin,0,0)
  sup <- c(Bmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_B else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_B,
    gr = gradient,
    WS = WS, SW = SW, GC = GC,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(n,length(init),SSres)
  return( list(
    "param"=list("B"=minSSE$par[1],
                 "e1"=minSSE$par[2],
                 "e2"=minSSE$par[3]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}


#' @title Sum of squares minimization of model BM
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: GC content
#' @param Bmin: minimum for the range of B, default value = -100
#' @param Bmax: maximum for the range of B, default value = 100
#' @param Mmin: minimum for the range of M, default value = -5
#' @param Mmax: maximum or the range of M, default value = 5
#' @param Maxit: maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm: number of updates in the method (option for optim, see manual)
#' @param Verbose: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr: to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @return A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square(sfsWS,sfsSW,0.5)
#' LS$B # population-scaled gBGC
#' LS$mutbias # mutation bias


least_square_BM <- function(WS,SW,GC,
                            Bmin=BMIN,Bmax=BMAX,Mmin=MMIN,Mmax=MMAX,
                            Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  n <- length(WS)
  x <- c(1:n)/(n+1)
  reginit <- lm(log(WS/SW) ~ x)
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  Binit <- reginit$coef[2]
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  e1max <- min(SW/(SW + rev(WS))) * ONE
  e2max <- min(WS/(WS + rev(SW))) * ONE
  # The factor ONE is to avoid 0 values in the SFS
  init <- c(Binit,Minit,e1max/2,e2max/2)
  # Boundaries for optimization
  inf <- c(Bmin,Mmin,0,0)
  sup <- c(Bmax,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_BM else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_BM,
    gr = gradient,
    WS = WS, SW = SW, GC = GC,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(n,length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B"=minSSE$par[1],
                 "mutbias"=exp(minSSE$par[2]),
                 "e1"=minSSE$par[3],
                 "e2"=minSSE$par[4]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}


#' @title Sum of squares minimization of model Hotpsot 1
#'
#' @description Function that searches for the six parameters that minimize the sum_of_squares function
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: GC content
#' @param Bmin: minimum for the range of B0, default value = -100
#' @param Bmax: maximum for the range of B0, default value = 100
#' @param Mmin: minimum for the range of M, default value = -5
#' @param Mmax: maximum or the range of M, default value = 5
#' @param Maxit: maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm: number of updates in the method (option for optim, see manual)
#' @param Verbose: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr: to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @return A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_hotspot2(sfsWS,sfsSW,0.5)
#' LS$B # population-scaled hotspot gBGC
#' LS$f # proportion of hotspots


least_square_hotspot1 <- function(WS,SW,GC,
                                  Bmin=BMIN,Bmax=BMAX,Mmin=MMIN,Mmax=MMAX,
                                  Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  n <- length(WS)
  x <- c(1:n)/(n+1)
  # Arbitray starting value fo f. Note that we assume 0 ≤ f ≤ 1/2
  finit <- 0.1
  reginit <- lm(log(WS/SW) ~ x)
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  Binit <- reginit$coef[2]/finit
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  e1max <- min(SW/(SW + rev(WS))) * ONE
  e2max <- min(WS/(WS + rev(SW))) * ONE
  # The factor ONE is to avoid 0 values in the SFS
  init <- c(Binit,finit,Minit,e1max/2,e2max/2)
  # Boundaries for optimization
  inf <- c(Bmin,0,Mmin,0,0)
  sup <- c(Bmax,1/2,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_hotspot1 else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_hotspot1,
    gr = gradient,
    WS = WS, SW = SW, GC = GC,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(n,length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B"=minSSE$par[1],
                 "f"=minSSE$par[2],
                 "mutbias"=exp(minSSE$par[3]),
                 "e1"=minSSE$par[4],
                 "e2"=minSSE$par[5]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}



#' @title Sum of squares minimization of model Hotpsot 2
#'
#' @description Function that searches for the six parameters that minimize the sum_of_squares function
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: GC content
#' @param B0min: minimum for the range of B0, default value = -100
#' @param B0max: maximum for the range of B0, default value = 100
#' @param B1min: minimum for the range of B1, default value = -100
#' @param B1max: maximum for the range of B1, default value = 100
#' @param Mmin: minimum for the range of M, default value = -5
#' @param Mmax: maximum or the range of M, default value = 5
#' @param Maxit: maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm: number of updates in the method (option for optim, see manual)
#' @param Verbose: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr: to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @return A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_hotspot2(sfsWS,sfsSW,0.5)
#' LS$B0 # population-scaled background gBGC
#' LS$B1 # population-scaled hotspot gBGC
#' LS$f # proportion of hotspots


least_square_hotspot2 <- function(WS,SW,GC,
                            B0min=BMIN,B0max=BMAX,B1min=BMIN,B1max=BMAX,Mmin=MMIN,Mmax=MMAX,
                            Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  n <- length(WS)
  x <- c(1:n)/(n+1)
  # Arbitray starting value fo f. Note that we assume 0 ≤ f ≤ 1/2
  finit <- 0.1
  reginit <- lm(log(WS/SW) ~ x)
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  # Starting values such that B1 = 5*B0
  B0init <- reginit$coef[2]/(1 + 4*finit)
  B1init <- 5*reginit$coef[2]/(1 + 4*finit)
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  e1max <- min(SW/(SW + rev(WS))) * ONE
  e2max <- min(WS/(WS + rev(SW))) * ONE
  # The factor ONE is to avoid 0 values in the SFS
  init <- c(B0init,B1init,finit,Minit,e1max/2,e2max/2)
  # Boundaries for optimization
  inf <- c(B0min,B1min,0,Mmin,0,0)
  sup <- c(B0max,B1max,1/2,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_hotspot2 else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_hotspot2,
    gr = gradient,
    WS = WS, SW = SW, GC = GC,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(n,length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B0"=minSSE$par[1],
                 "B1"=minSSE$par[2],
                 "f"=minSSE$par[3],
                 "mutbias"=exp(minSSE$par[4]),
                 "e1"=minSSE$par[5],
                 "e2"=minSSE$par[6]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}



#' @title Sum of squares minimization of model Hotpsot 2bis
#'
#' @description Function that searches for the five parameters that minimize the sum_of_squares function
#' This model is the same as model 2 except that f is fixed by the user
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: GC content
#' @param f: proportion of hotspots (fixed by the user: 0 ≤ f ≤ 1/2)
#' @param B0min: minimum for the range of B0, default value = -100
#' @param B0max: maximum for the range of B0, default value = 100
#' @param B1min: minimum for the range of B1, default value = -100
#' @param B1max: maximum for the range of B1, default value = 100
#' @param Mmin: minimum for the range of M, default value = -5
#' @param Mmax: maximum or the range of M, default value = 5
#' @param Maxit: maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm: number of updates in the method (option for optim, see manual)
#' @param Verbose: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr: to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @return A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_hotspot2(sfsWS,sfsSW,0.5)
#' LS$B0 # population-scaled background gBGC
#' LS$B1 # population-scaled hotspot gBGC
#' LS$f # proportion of hotspots


least_square_hotspot2bis <- function(WS,SW,GC,f,
                                  B0min=BMIN,B0max=BMAX,B1min=BMIN,B1max=BMAX,Mmin=MMIN,Mmax=MMAX,
                                  Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  if(f<0 | f>1/2) {
    print("ERROR: GC content must be between 0 and 1/2")
    return(NA)
  }
  n <- length(WS)
  x <- c(1:n)/(n+1)
  # Arbitray starting value fo f. Note that we assume 0 ≤ f ≤ 1/2
  reginit <- lm(log(WS/SW) ~ x)
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  # Starting values such that B1 = 5*B0
  B0init <- reginit$coef[2]/(1 + 4*f)
  B1init <- 5*reginit$coef[2]/(1 + 4*f)
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  e1max <- min(SW/(SW + rev(WS))) * ONE
  e2max <- min(WS/(WS + rev(SW))) * ONE
  # The factor ONE is to avoid 0 values in the SFS
  init <- c(B0init,B1init,Minit,e1max/2,e2max/2)
  # Boundaries for optimization
  inf <- c(B0min,B1min,Mmin,0,0)
  sup <- c(B0max,B1max,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_hotspot2bis else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_hotspot2bis,
    gr = gradient,
    WS = WS, SW = SW, GC = GC, f = f,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(n,length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B0"=minSSE$par[1],
                 "B1"=minSSE$par[2],
                 "mutbias"=exp(minSSE$par[3]),
                 "e1"=minSSE$par[4],
                 "e2"=minSSE$par[5]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}
