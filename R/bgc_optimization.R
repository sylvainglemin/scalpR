# Estimation of GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

source("R/constants.R")
source("R/bgc_SS_constant.R")

# Functions to minimize the sum of squares of the different models

# The mean of the sum of square function gives directly the mean sum of square
# The two functions are set equal for homogeneity in function call
least_square_M <- sum_of_squares_M


#' @title Sum of squares minimization of the model with constant gBGC but without mutation bias
#'
#' @description Function that searches for the three parameters that minimize the sum_of_squares function
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param Bmin: minimum for the range of B, default value = -100
#' @param Bmax: maximum for the range of B, default value = 100
#' @param Maxit: maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm: number of updates in the method (option for optim, see manual)
#' @param Verbose: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr: to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @return The list of optimized parameters
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square(sfsWS,sfsSW,0.5)
#' LS$B # population-scaled gBGC


least_square_B <- function(WS,SW,GC,Bmin=BMIN,Bmax=BMAX,
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
  n <- length(WS) + 1
  x <- c(1:(n-1))/n
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
  SStot <- sum_of_squares_M(SW,WS,GC)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  return( list("B"=minSSE$par[1],"e1"=minSSE$par[2],"e2"=minSSE$par[3],"SStot"=SStot,"SSres"=SSres,"R2"=R2) )
}


#' @title Sum of squares minimization of the model with constant gBGC and mutation bias
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
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
#' @return The list of optimized parameters
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square(sfsWS,sfsSW,0.5)
#' LS$B # population-scaled gBGC
#' LS$mutbias # mutation bias


least_square_BM <- function(WS,SW,GC,Bmin=BMIN,Bmax=BMAX,Mmin=MMIN,Mmax=MMAX,
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
  n <- length(WS) + 1
  x <- c(1:(n-1))/n
  reginit <- lm(log(WS/SW) ~ x)
  Minit <- -reginit$coef[1]*(1-GC)/GC
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
  SStot <- sum_of_squares_M(SW,WS,GC)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  return(
    list("B"=minSSE$par[1],
         "mutbias"=exp(minSSE$par[2]),
         "e1"=minSSE$par[3],"e2"=minSSE$par[4],
         "SStot"=SStot,
         "SSres"=SSres,
         "R2"=R2)
  )
}
