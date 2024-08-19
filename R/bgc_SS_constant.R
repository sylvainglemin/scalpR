# Estimation of GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

source("R/constants.R")


# Series of SS functions and their corresponding gradient where parameters are constant


######################################### #
# NULL MODEL ##############################
######################################### #

# This model correspond to a fit of the intercept only
# It gives the mean and the total sum of square of the dataset


#' @title Sum of squares of the null model
#'
#' @description Function that return the weighted total sum of squares for the model with only the intercept
#'
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#'
#' @return The weighted sum of squares
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' sum_of_squares_NULL(sfsWS,sfsSW)


sum_of_squares_NULL <- function(WS,SW,GC) {
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  w <- WS*SW/(WS + SW)
  y <- log(WS/SW)  - 1/(2*WS) + 1/(2*SW)
  removeNA <- !is.na(y)
  w <- w[removeNA]
  SW <- SW[removeNA]
  WS <- WS[removeNA]
  y_mean <- mean(y)
  return(
    list(
      "mean" = y_mean,
      "SStot"=sum(w*(y - y_mean)^2)/sum(w)
    )
  )
}


######################################### #
# NO BGC, MUTATION BIAS ###################
######################################### #

#' @title Sum of squares of model M
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = -log(mutbias) + log(GC) - log(1 - GC)
#' where Tws and Tsw are the "True" SFSs and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#'@param par: a vector with the for parameters of the model.
#' par[1] = M (i.e. log(mut_bias))
#' par[2] = e1
#' par[3] = e2
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: GC content
#'
#' @return The weighted sum of squares
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' sum_of_squares_M(sfsWS,sfsSW,0.5)

sum_of_squares_M <- function(par,WS,SW,GC) {
  if(length(par)!=3) {
    print("ERROR: a three values vector must be given as par")
    return(NA)
  }
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  M <- par[1]
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  # w <- WSt2*SWt2/(WSt2 + SWt2)
  # Here the weight is written as a function of corrected SFSs
  # This complexifies the whole equation as e1 and e2 appear in w
  # The gradient function is also more complicated
  # Instead we use the weight as a function of observed SFSs
  w <- WS*SW/(WS + SW)
  y <- log(WSt2/SWt2)  - 1/(2*WSt) + 1/(2*SWt)
  ypred <- rep(- M - log(GC) + log(1 - GC),n)
  removeNA <- !is.na(y)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}


#' @title Gradient of Sum of squares model M
#'
#' @description Gradient of the sum of squares function
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = M (i.e. log(mut_bias))
#' par[2] = e1
#' par[3] = e2
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: the GC content
#'
#' @return The gradient function
#'
#' @examples

gr_sum_of_squares_M <- function(par,WS,SW,GC) {
  if(length(par)!=3) {
    print("ERROR: a three values vector must be given as par")
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
  M <- par[1]
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1-e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1-e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- WS*SW/(WS + SW)
  x <- c(1:(n-1))/n
  y <- log(WSt2/SWt2) - 1/(2*WSt) + 1/(2*SWt)
  ypred <- rep(- M - log(GC) + log(1 - GC),n)
  ## TO BE COMPLETED
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- (SW + rev(WS))/((1 - e1)*SW - e1*rev(WS)) + # Log term
    1/(2*((1 - e2)*WS - e2*rev(SW))) + # First ratio term
    ((1 - e2)*rev(WS) - e2*SW)/(2*((1 - e1)*SW - e1*rev(WS))^2) # Second ratio term
  dy2 <- (rev(SW) + WS)/(e2*rev(SW) - (1 - e2)*WS) + # Log term
    (e1*WS - (1 - e1)*rev(SW))/(2*((1 - e2)*WS - e2*rev(SW))^2) + # First ratio term
    1/(2*(e1*rev(WS) - (1 - e1)*SW))
  removeNA <- !is.na(y)
  w <- w[removeNA]
  x <- x[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grM,gre1,gre2))
}



######################################### #
# CONSTANT BGC, NO MUTATION BIAS ##########
######################################### #


#' @title Sum of squares of model B
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = B * j/n + log(GC) - log(1 - GC)
#' where Tws and Tsw are the "True" SFSs, B = 4Neb is the population-scaled gBGC, and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par: a vector with the three parameters of the model.
#' par[1] = B
#' par[2] = e1
#' par[3] = e2
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
#' sum_of_squares(param,sfsWS,sfsSW,0.5)

sum_of_squares_B <- function(par,WS,SW,GC) {
  if(length(par)!=3) {
    print("ERROR: a three values vector must be given as par")
    return(NA)
  }
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  B <- par[1]
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  # w <- WSt2*SWt2/(WSt2 + SWt2)
  # Here the weight is written as a function of corrected SFSs
  # This complexifies the whole equation as e1 and e2 appear in w
  # The gradient function is also more complicated
  # Instead we use the weight as a function of observed SFSs
  w <- WS*SW/(WS + SW)
  x <- c(1:n)/(n+1)
  y <- log(WSt2/SWt2)  - 1/(2*WSt) + 1/(2*SWt)
  ypred <- B*x - log(GC) + log(1 - GC)
  removeNA <- !is.na(y)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}


#' @title Gradient of Sum of squares of model B
#'
#' @description Gradient of the sum of squares function
#'
#' @param par: a vector with the for parameters of the model.
#' par[1] = B
#' par[2] = e1
#' par[3] = e2
#' @param WS: the WS observed SFS
#' @param SW: the SW observed SFS
#' @param GC: the GC content
#'
#' @return The gradient function
#'
#' @examples

gr_sum_of_squares_B <- function(par,WS,SW,GC) {
  if(length(par)!=3) {
    print("ERROR: a three values vector must be given as par")
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
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1-e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1-e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- WS*SW/(WS + SW)
  x <- c(1:n)/(n+1)
  y <- log(WSt2/SWt2) - 1/(2*WSt) + 1/(2*SWt)
  ypred <- B*x - log(GC) + log(1 - GC)
  ## TO BE COMPLETED
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- (SW + rev(WS))/((1 - e1)*SW - e1*rev(WS)) + # Log term
    1/(2*((1 - e2)*WS - e2*rev(SW))) + # First ratio term
    ((1 - e2)*rev(WS) - e2*SW)/(2*((1 - e1)*SW - e1*rev(WS))^2) # Second ratio term
  dy2 <- (rev(SW) + WS)/(e2*rev(SW) - (1 - e2)*WS) + # Log term
    (e1*WS - (1 - e1)*rev(SW))/(2*((1 - e2)*WS - e2*rev(SW))^2) + # First ratio term
    1/(2*(e1*rev(WS) - (1 - e1)*SW))
  removeNA <- !is.na(y)
  w <- w[removeNA]
  x <- x[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  grB <- sum(w*(-2*x*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB,gre1,gre2))
}


######################################### #
# CONSTANT BGC, MUTATION BIAS #############
######################################### #

#' @title Sum of squares of model BM
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = B * j/n - log(mut_bias) + log(GC) - log(1 - GC)
#' where Tws and Tsw are the "True" SFSs, B = 4Neb is the population-scaled gBGC, mut_bias the mutation bias towards AT, and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par: a vector with the four parameters of the model.
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
#' sum_of_squares(param,sfsWS,sfsSW,0.5)

sum_of_squares_BM <- function(par,WS,SW,GC) {
  if(length(par)!=4) {
    print("ERROR: a four values vector must be given as par")
    return(NA)
  }
  if(length(WS)!=length(SW)) {
    print("ERROR: the two SFSs, WS and SW, must have the same length")
    return(NA)
  }
  if(GC<=0 | GC>=1) {
    print("ERROR: GC content must be strictly between 0 and 1")
    return(NA)
  }
  B <- par[1]
  M <- par[2]
  e1 <- par[3]
  e2 <- par[4]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  # w <- WSt2*SWt2/(WSt2 + SWt2)
  # Here the weight is written as a function of corrected SFSs
  # This complexifies the whole equation as e1 and e2 appear in w
  # The gradient function is also more complicated
  # Instead we use the weight as a function of observed SFSs
  w <- WS*SW/(WS + SW)
  x <- c(1:n)/(n+1)
  y <- log(WSt2/SWt2)  - 1/(2*WSt) + 1/(2*SWt)
  ypred <- B*x - M - log(GC) + log(1 - GC)
  removeNA <- !is.na(y)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}


#' @title Gradient of Sum of squares of model BM
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

gr_sum_of_squares_BM <- function(par,WS,SW,GC) {
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
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1-e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1-e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- WS*SW/(WS + SW)
  x <- c(1:n)/(n+1)
  y <- log(WSt2/SWt2) - 1/(2*WSt) + 1/(2*SWt)
  ypred <- B*x - M - log(GC) + log(1 - GC)
## TO BE COMPLETED
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- (SW + rev(WS))/((1 - e1)*SW - e1*rev(WS)) + # Log term
    1/(2*((1 - e2)*WS - e2*rev(SW))) + # First ratio term
    ((1 - e2)*rev(WS) - e2*SW)/(2*((1 - e1)*SW - e1*rev(WS))^2) # Second ratio term
  dy2 <- (rev(SW) + WS)/(e2*rev(SW) - (1 - e2)*WS) + # Log term
    (e1*WS - (1 - e1)*rev(SW))/(2*((1 - e2)*WS - e2*rev(SW))^2) + # First ratio term
    1/(2*(e1*rev(WS) - (1 - e1)*SW))
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
