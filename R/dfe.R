# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr



# TO DO: CHANGE FUNCTIONS NAMES WHICH ARE THE SAME AS FOR BGC



# Libraries
library(hypergeo)
library(VGAM)


# Constants used when functions are not defined in 0
ZERO <- 10^(-6)
ONE <- 1-ZERO

#' @title Expected number of deleterious SNPs in frequency i/n - approximation for low frequencies
#'
#' @description Function that gives the expected ratio of the non-synonymous over synonymous SFS
#' under a gamma distribution of deleterious mutations. Accurate approximations for low frequencies
#'
#' @param Sdel: the population-scale mean effect of deleterious mutations
#' @param shape: the shape of the distribution
#' @param n: the sample size
#' @param i: snp category
#'
#'
#' @return Expected number of deleterious SNPs in frequency i/n normalized by the neutral expectation
#' @export
#'
#' @examples
#' sfs_del_low(1000,0.5,20,2)
#'
sfs_del_low <- function(Sdel,shape,n,i) {
  sfsi <- n*hypergeo(A = i,B = shape,C = n, z = -Sdel/shape) / (n - i) -
    (n*shape + i*Sdel)*zeta(1+shape)*(shape/Sdel)^(shape+1) / (n - i)
  return(sfsi)
}

#' @title Expected number of deleterious SNPs in frequency i/n - approximation for high frequencies
#'
#' @description Function that gives the expected ratio of the non-synonymous over synonymous SFS
#' under a gamma distribution of deleterious mutations. Accurate approximations for low frequencies
#' We use the third approximation from the Mathematica notebook
#'
#' @param Sdel: the population-scale mean effect of deleterious mutations
#' @param shape: the shape of the distribution
#' @param n: the sample size
#' @param i: snp category
#'
#'
#' @return Expected number of deleterious SNPs in frequency i/n normalized by the neutral expectation
#' @export
#'
#' @examples
#' sfs_del_high(1000,0.5,20,15)
#'
sfs_del_high <- function(Sdel,shape,n,i) {
  sfsi <- shape*zeta(1+shape)*(shape/Sdel)^shape*
    (i/n)^(-(1+shape)*zeta(2+shape)/(2*zeta(1+shape)))
  return(sfsi)
}
# # Less approximate expression but not necessarily better
# sfs_del_high <- function(Sdel,shape,n,i) {
#   sfsi <- (1+shape)*log(shape) - shape*log(Sdel) +lgamma(n+1) +log(zeta(1+shape)) +
#     lgamma(i-(1+shape)*zeta(2+shape)/(2*zeta(1+shape))) -
#     lgamma(i) -lgamma(n+1-(1+shape)*zeta(2+shape)/(2*zeta(1+shape)))
#   return(exp(sfsi))
# }

#' @title Expected number of advantageous SNPs in frequency i/n - approximation for high frequencies
#'
#' @description Function that gives the expected ratio of the non-synonymous over synonymous SFS
#' under an exponential distribution of advantageous mutations. Accurate approximations for low frequencies
#' We use the third approximation from the Matematica notebook
#'
#' @param Sadv: the population-scale mean effect of deleterious mutations
#' @param n: the sample size
#' @param i: snp category
#'
#' @return Expected number of deleterious SNPs in frequency i/n normalized by the neutral expectation
#' @export
#'
#' @examples
#' sfs_adv(10,20,15)
#'
sfs_adv_i <- function(Sadv,n,i) {
  sfsi <- n*(digamma(1 - i/n + 1/Sdel)-digamma(1/Sadv))/((n-i)*Sadv)
  return(sfsi)
}

#' @title SFS of deleterious effects of mutations
#'
#' @description Function that gives the expected ratio of the non-synonymous over synonymous SFS
#' under a gamma distribution of deleterious mutations
#'
#' @param Sdel: the population-scale mean effect of deleterious mutations
#' @param shape: the shape of the distribution
#' @param n: the sample size
#'
#' @return The expected deleterious non-synonymous SFS normalized byt the neutral SFS
#' @export
#'
#' @examples
#' sfs_del(1000,0.5,20)
#'
# Function with a threshold at n/2
sfs_del <- function(Sdel,shape,n) {
  sfs <- sapply(c(1:(n-1)), function(i) ifelse(i<=n/2,sfs_del_low(Sdel,shape,n,i),sfs_del_high(Sdel,shape,n,i)))
  return(log(Re(sfs)))
}
# Continuous weighting
sfs_del <- function(Sdel,shape,n,sm=20) {
  w <- function(s,x) exp(-s*(x-1/2)) /(1+exp(-s*(x-1/2)))
  sfs <- sapply(c(1:(n-1)), function(i)
    sfs_del_low(Sdel,shape,n,i)*w(sm,i/n)+
      sfs_del_high(Sdel,shape,n,i)*w(sm,1-i/n) )
  return(log(Re(sfs)))
}

#' @title SFS of advantageous effects of mutations
#'
#' @description Function that gives the expected ratio of the non-synonymous over synonymous SFS
#' under a gamma distribution of deleterious mutations
#'
#' @param Sadv: the population-scale mean effect of advantageous mutations
#' @param n: the sample size
#'
#' @return The expected advantageous non-synonymous SFS normalized by the neutral SFS
#' @export
#'
#' @examples
#' sfs_del(1000,0.5,20)
#'
sfs_adv <- function(Sadv,n) {
  sfs <- sapply(c(1:(n-1)), function(i) sfs_adv_i(Sadv,n,i))
  return(log(Re(sfs)))
}


#' @title SFS total
#'
#' @description Function that gives the expected ratio of the non-synonymous over synonymous SFS
#' under a gamma distribution of deleterious mutations
#'
#' @param Sdel: the population-scale mean effect of deleterious mutations
#' @param shape: the shape of the distribution
#' @param Sadv: the population-scale mean effect of advantageous mutations
#' @param p: proportion of advantageous mutations
#' @param n: the sample size
#'
#' @return The expected non-synonymous SFS normalized by the neutral SFS
#' @export
#'
#' @examples
#'
sfs_tot <- function(Sdel,shape,Sadv,p,n) {
  sfs <- (1-p)*exp(sfs_del(Sdel,shape,n)) + p*exp(sfs_adv(Sadv,n))
  return(log(Re(sfs)))
}

#' @title Sum of squares - model with deleterious mutations only
#'
#' @description
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = Sdel
#' par[2] = shape
#' par[3] = e
#' @param qns: the ratio of -the length of non-synonymous over synonymous positions
#' @param syn: the synonymous observed SFS
#' @param nonsyn: the non-synonymous observed SFS
#'
#' @return The weighted sum of square
#' @export
#'
#' @examples
#'
sum_of_squares_del <- function(par,qns,syn,nonsyn) {
  if(length(par)!=3) {
    print("ERROR: a three values vector  must be given as par")
    return(NA)
  }
  if(length(syn)!=length(nonsyn)) {
    print("ERROR: the two SFSs must have the same length")
    return(NA)
  }
  Sdel <- 10^par[1]
  shape <- par[2]
  e <- par[3]
  n <- length(syn)+1
  # True SFS as a function of observed one.
  synt <- ((1-e)*syn-e*rev(syn))/(1-2*e)
  # To avoid negative values
  synt <- sapply(synt, function(x) max(x,0))
  nonsynt <- ((1-e)*nonsyn-e*rev(nonsyn))/(1-2*e)
  nonsynt <- sapply(nonsynt, function(x) max(x,0))
  # Same expression without 1-2e that simplifies in ratios
  synt2 <- ((1-e)*syn-e*rev(syn))
  synt2 <- sapply(synt2, function(x) max(x,0))
  nonsynt2 <- ((1-e)*nonsyn-e*rev(nonsyn))
  nonsynt2 <- sapply(nonsynt2, function(x) max(x,0))
  w <- synt2*nonsynt2/(synt2+nonsynt2)
  y <- log(nonsynt2/synt2) - 1/(2*nonsynt2) + 1/(2*synt2)
  ypred <- sfs_del(Sdel,shape,n) + log(qns)
  removeNA <- !is.na(y) & !is.infinite(y)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}

#' @title Sum of squares minimization
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = Sdel
#' par[2] = shape
#' par[3] = e
#' @param qns: the ratio of -the length of non-synonymous over synonymous positions
#' @param syn: the syn observed SFS
#' @param nonsyn: the nonsyn observed SFS
#' @param Sdmin: minimum for the range of Sdel, default value = 10
#' @param Sdmax: maximum for the range of Sdel, default value = 10^8
#' @param smin: minimum for the range of shape, default value = 0.01
#' @param smax: maximum or the range of shape, default value = 10
#' @param MAXIT: maximum number of iterations (option for optim), see manual, default value = 100
#' @param FACTR: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param LMM: number of updates in the method (option for optim, see manual)
#' @param VERBOSE: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#'
#' @return The list of optimized parameters
#' @export
#'
#' @examples
least_square_del <- function(qns,syn,nonsyn,Sdmin=1,Sdmax=8,smin=0.01,smax=10,MAXIT=100,FACTR=10^7,LMM=20,VERBOSE=0) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(syn)!=length(nonsyn)) {
    print("ERROR: the two SFSs must have the same length")
    return(NA)
  }
  n <- length(syn)+1
  #To obtain starting value: regression on the first half part of the sfs
  n2 <- round(n/2)
  x <- log(c(1:(n-1))/n)[c(1:n2)]
  reginit <- lm(log(nonsyn[c(1:n2)]/syn[c(1:n2)]) ~ x,na.action = "na.omit")
  shape_init <- max(-reginit$coef[2],0.1)
  names(shape_init)<-NULL
  Sdel_init <- log(shape_init) -reginit$coef[1]/shape_init # Sdel_init <- log(shape_init) - reginit$coef[1]/shape_init
  Sdel_init <- ifelse(Sdel_init>6,6,ifelse(Sdel_init<1,1,Sdel_init))
  names(Sdel_init) <- NULL
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  # The factor ONE is to avoid 0 values in the SFS
  e_s <- min(syn/(syn + rev(syn)),na.rm = T) * ONE
  e_ns <- min(nonsyn/(nonsyn + rev(nonsyn)),na.rm = T) * ONE
  e_max <- min(e_s,e_ns)
  init <- c(Sdel_init,shape_init,e_max/2)
  # Boundaries for optimization
  inf <- c(Sdmin,smin,0)
  sup <- c(Sdmax,smax,e_max)
  minSSE <- optim(init,sum_of_squares_del,qns=qns,syn=syn,nonsyn=nonsyn,lower = inf,upper = sup,method="L-BFGS-B",control=list(parscale=abs(init),maxit=MAXIT,factr=FACTR,lmm=LMM,trace=VERBOSE))
  return( list("SSE"=minSSE$value,"Sdel"=10^(minSSE$par[1]),"shape"=minSSE$par[2],"e"=minSSE$par[3]) )
}

#' @title Sum of squares - total
#'
#' @description
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = Sdel
#' par[2] = shape
#' par[3] = Sadv
#' par[4] = p
#' par[5] = e
#' @param syn: the syn observed SFS
#' @param nonsyn: the nonsyn observed SFS
#'
#' @return The weighted sum of square
#' @export
#'
#' @examples
#'
sum_of_squares_tot <- function(par,qns,syn,nonsyn) {
  if(length(par)!=5) {
    print("ERROR: a five values vector must be given as par")
    return(NA)
  }
  if(length(syn)!=length(nonsyn)) {
    print("ERROR: the two SFSs must have the same length")
    return(NA)
  }
  Sdel <- 10^par[1]
  shape <- par[2]
  Sadv <- 10^par[3]
  p <- par[4]
  e <- par[5]
  n <- length(syn)+1
  # True SFS as a function of observed one.
  synt <- ((1-e)*syn-e*rev(syn))/(1-2*e)
  synt <- sapply(synt, function(x) max(x,0))
  nonsynt <- ((1-e)*nonsyn-e*rev(nonsyn))/(1-2*e)
  nonsynt <- sapply(nonsynt, function(x) max(x,0))
  # Same expression without 1-2e that simplifies in ratios
  synt2 <- ((1-e)*syn-e*rev(syn))
  nonsynt2 <- ((1-e)*nonsyn-e*rev(nonsyn))
  synt2 <- ((1-e)*syn-e*rev(syn))
  synt2 <- sapply(synt2, function(x) max(x,0))
  nonsynt2 <- ((1-e)*nonsyn-e*rev(nonsyn))
  nonsynt2 <- sapply(nonsynt2, function(x) max(x,0))
  w <- synt2*nonsynt2/(synt2+nonsynt2)
  w <- synt2*nonsynt2/(synt2+nonsynt2)
  y <- log(nonsynt2/synt2) - 1/(2*nonsynt2) + 1/(2*synt2)
  ypred <- sfs_tot(Sdel,shape,Sadv,p,n) + log(qns)
  removeNA <- !is.na(y) & !is.infinite(y)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}

#' @title Sum of squares minimization
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param par: a vector with the for parameters of the model.
#' @param syn: the syn observed SFS
#' @param nonsyn: the nonsyn observed SFS
#' @param Sdmin: minimum for the range of Sdel, default value = 10
#' @param Sdmax: maximum for the range of Sdel, default value = 10^8
#' @param Samin: minimum for the range of Sadv, default value = 10^(-3)
#' @param Samax: maximum for the range of Sadv, default value = 10^5
#' @param smin: minimum for the range of shape, default value = 0.01
#' @param smax: maximum or the range of shape, default value = 10
#' @param MAXIT: maximum number of iterations (option for optim), see manual, default value = 100
#' @param FACTR: level of control the convergence (option for optim, see manual), default value = 10^7
#' @param LMM: number of updates in the method (option for optim, see manual)
#' @param VERBOSE: from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#'
#' @return The list of optimized parameters
#' @export
#'
#' @examples
least_square_tot <- function(qns,syn,nonsyn,Sdmin=1,Sdmax=8,Samin=-3,Samax=5,smin=0.01,smax=10,MAXIT=100,FACTR=10^7,LMM=20,VERBOSE=0) {
  # Determination of initial values for optimization: use of the simple regression
  if(length(syn)!=length(nonsyn)) {
    print("ERROR: the two SFSs must have the same length")
    return(NA)
  }
  #For testing
  # Sdmin=1;Sdmax=8;Samin=-3;Samax=5;smin=0.01;smax=10;MAXIT=100;FACTR=10^7;LMM=20;VERBOSE=5
  # nsample <- 20
  # syn <- sapply(10000/c(1:(n_sample-1)),function(x) rpois(1,x))
  # nonsyn <- sapply(10000*exp(sfs_tot(1000,0.3,1,0.01,n_sample))/c(1:(n_sample-1)),function(x) rpois(1,x))
  #To obtain starting value: regression on the first half part of the sfs
  n <- length(syn)+1
  n2 <- round(n/2)
  x <- log(c(1:(n-1))/n)[c(1:n2)]
  reginit <- lm(log(nonsyn[c(1:n2)]/syn[c(1:n2)]) ~ x)
  shape_init <- max(-reginit$coef[2],0.1)
  names(shape_init)<-NULL
  Sdel_init <- log(shape_init) -reginit$coef[1]/shape_init # Sdel_init <- log(shape_init) - reginit$coef[1]/shape_init
  Sdel_init <- min(max(Sdel_init,6),1)
  names(Sdel_init) <- NULL
  p_init <- 0.01 # Should not be 0 to avoid problem with the optim option: parscale=abs(init)
  Sadv_init <- log(max(0.5,(nonsyn[n-1]/syn[n-1])/p_init))
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  # The factor ONE is to avoid 0 values in the SFS
  e_s <- min(syn/(syn + rev(syn)),na.rm = T) * ONE
  e_ns <- min(nonsyn/(nonsyn + rev(nonsyn)),na.rm = T) * ONE
  e_max <- min(e_s,e_ns)
  init <- c(Sdel_init,shape_init,Sadv_init,p_init,e_max/2)
  # Boundaries for optimization
  inf <- c(Sdmin,smin,Samin,0,0)
  sup <- c(Sdmax,smax,Samax,1,e_max)
  minSSE <- optim(init,sum_of_squares_tot,qns=qns,syn=syn,nonsyn=nonsyn,lower = inf,upper = sup,method="L-BFGS-B",control=list(parscale=abs(init),maxit=MAXIT,factr=FACTR,lmm=LMM,trace=VERBOSE))
  return( list("SSE"=minSSE$value,"Sdel"=10^(minSSE$par[1]),"shape"=minSSE$par[2],"Sadv"=10^(minSSE$par[3]),"p"=minSSE$par[4],"e"=minSSE$par[5]) )
}


# TEST ####

theta <- 100000
n_sample <- 10
Sdel <- 1000
b <- .3
Sadv <- 10
p <- 0.01
sfs_s <- sapply(theta/c(1:(n_sample-1)),function(x) rpois(1,x))
sfs_ns <- sapply(2*theta*exp(sfs_tot(Sdel,b,Sadv,p,n_sample))/c(1:(n_sample-1)),function(x) rpois(1,x))
plot(sfs_ns/sfs_s)

least_square_del(qns = 2,syn = sfs_s,nonsyn = sfs_ns)
least_square_tot(qns = 2,syn = sfs_s,nonsyn = sfs_ns)

