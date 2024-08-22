# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

# TO DO CLEAN AND COMPLETE

# Various functions to manipulate SFS or compute statistics


# Harmonic number
hn <- function(n) sum(1/c(1:n))



# Function to compute the skweness of a frequency spectrum and testing the asymmetry
# Typically for testing the asymmetry of the GC unfolded spectrum

skewness_sfs <- function(sfs) {
  Nclass <- length(sfs)+1
  Nobs <- sum(sfs)
  freq <- c(1:(Nclass-1))/Nclass
  mean <- sum(sfs*freq)/Nobs
  skew <- (sqrt(Nobs*(Nobs-1))/(Nobs-2))*sum(sfs*(freq-mean)^3/Nobs)/(sum(sfs*(freq-mean)^2)/Nobs)^(3/2)
  SES <- sqrt(6*Nobs*(Nobs-1)/((Nobs-2)*(Nobs+1)*(Nobs+3))) # Standard error of skewness
  pval <- 2*pnorm(abs(skew/SES),0,1,lower.tail=F)
  return(list(skewness=skew,p.value=pval))
}


# Fonction to project the SFS (as a vector) from n to m < n chromosomes
# Return the reduced SFS vector

project_sfs <- function(sfs,m){
  n <- length(sfs)
  output <- rep(0,m)
  index<-c(1:m)
  for(k in 1:n){
    output <- output + sfs[k]*choose(k,index)*choose(n-k,m-index)/choose(n,m)
  }
  return(output)
}

