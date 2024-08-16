# Estimation of GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

# Constants used when functions are not defined in 0
ZERO <- 10^(-6)
ONE <- 1 - ZERO


#' @title Sum of squares
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = B * j/n - log(mut_bias)
#' where Tws and Tsw are the "True" SFSs, B = 4Neb is the population-scaled gBGC, mut_bias the mutation bias towards AT, and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = B
#' par[2] = M (i.e. log(mut_bias))
#' par[3] = e1
#' par[4] = e2
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: the GC content
#'
#' @return The weighted sum of squares
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' param <- c(1,2,0.02,0.01)
#' sum_of_squares(param,sfsWS,sfsSW)

sum_of_squares <- function(par,WS,SW,GC) {
  if(length(par)!=4) {
    print("ERROR: a four values vector must be given as par")
    return(NA)
  }
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be between 0 and 1")
    return(NA)
  }
  B <- par[1]
  M <- par[2]
  e1 <- par[3]
  e2 <- par[4]
  n <- length(WS) + 1
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression withou 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- WSt2*SWt2/(WSt2 + SWt2)
  x <- c(1:(n-1))/n
  y <- log(WSt2/SWt2) # - 1/(2*WSt) + 1/(2*SWt)
  ypred <- B*x - M - log(GC) + log(1 - GC)
  removeNA <- !is.na(y)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}

# USE OF THE GRADIENT --> MUCH LONGER!
# TO IMPROVE OR TO REMOVE

#' @title Gradient of Sum of squares
#'
#' @description Gradient of the sum of squares function
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = B
#' par[2] = M (i.e. log(mut_bias))
#' par[3] = e1
#' par[4] = e2
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: the GC content
#'
#' @return The gradient function
#'
#' @examples

gr_sum_of_squares <- function(par,WS,SW,GC) {
  if(length(par)!=4) {
    print("ERROR: a four values vector must be given as par")
    return(NA)
  }
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be between 0 and 1")
    return(NA)
  }
  B <- par[1]
  M <- par[2]
  e1 <- par[3]
  e2 <- par[4]
  n <- length(WS) + 1
  # True SFS as a function of observed one.
  WSt <- ((1-e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1-e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression withou 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- WSt2*SWt2/(WSt2 + SWt2)
  x <- c(1:(n-1))/n
  y <- log(WSt2/SWt2) # - 1/(2*WSt) + 1/(2*SWt)
  ypred <- B*x - M - log(GC) + log(1 - GC)
## TO BE COMPLETED
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- (SW + rev(WS))/((1 - e1)*SW - e1*rev(WS))
  dy2 <- (rev(SW) + WS)/(e2*rev(SW) - (1 - e2)*WS)
  removeNA <- !is.na(y)
  w <- w[removeNA]
  x <- x[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  grB <- sum(w*(-2*x*(y - ypred))/sum(w))
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB,grM,gre1,gre2))
}




#' @title Sum of squares minimization
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = B
#' par[2] = M (i.e. log(mut_bias))
#' par[3] = e1
#' par[4] = e2
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param Bmin: minimum for the range of B, default value = -100
#' @param Bmax: maximum for the range of B, default value = 100
#' @param Mmin: minimum for the range of M, default value = -5
#' @param Mmax: maximum or the range of M, default value = 5
#' @param MAXIT: maximum number of iterations (option for optim), see manual, default value = 100
#' @param FACTR: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param LMM: number of updates in the method (option for optim, see manual)
#' @param VERBOSE: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#'
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


least_square <- function(WS,SW,GC,Bmin=-100,Bmax=100,Mmin=-5,Mmax=5,MAXIT=100,FACTR=10^7,LMM=20,VERBOSE=0,USEGR=F) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  n <- length(WS) + 1
  x <- c(1:(n-1))/n
  reginit <- lm(log(WS/SW) ~ x)
  Minit <- -reginit$coef[1]*(1-GC)/GC
  Binit <- reginit$coef[2]
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  # The factor ONE is to avoid 0 values in the SFS
  e1max <- min(SW/(SW + rev(WS))) * ONE
  e2max <- min(WS/(WS + rev(SW))) * ONE
  init <- c(Binit,Minit,e1max/2,e2max/2)
  # Boundaries for optimization
  inf <- c(Bmin,Mmin,0,0)
  sup <- c(Bmax,Mmax,e1max,e2max)
  if(USEGR) gradient <- gr_sum_of_squares else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares,
    gr = gradient,
    WS = WS, SW = SW, GC = GC,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=MAXIT,factr=FACTR,lmm=LMM,trace=VERBOSE))
  return( list("B"=minSSE$par[1],"mutbias"=exp(minSSE$par[2]),"e1"=minSSE$par[3],"e2"=minSSE$par[4],"val"=minSSE$value))
}
