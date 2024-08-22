# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

source("R/constants.R")

# COMMENTS TO ADD TO EACH FUNCTION

# Functions giving the expected number of SNP (ens) in frequency j/n under different models


# Neutrality
ens_neutral <- function(theta,j) {
	theta/j
}

# Constant gBGC with intensity B. B can be positive (W->S mutations) or negative (S->W) mutations
ens_constant <- function(theta,B,n,j) {
	if (B==0) return(ens_neutral(theta,j))
	((n*theta)/(j*(n-j))) * ((1-exp(-B*(1-j/n)))/(1-exp(-B)))
}

# Hotspot1: gBGC with intensity B in hotspot (fraction f) and 0 otherwise (fraction 1-f)
ens_hotspot1 <- function(theta,B,f,n,j) {
	if (B==0) return(ens_neutral(theta,j))
	f * ens_constant(theta,B,n,j) + (1-f) * ens_neutral(theta,j)
}

# Hotspot2: gBGC with intensity B1 in hotspot (fraction f) and B0 otherwise (fraction 1-f)
ens_hotspot2 <- function(theta,B0,B1,f,n,j) {
	if (B0==0) return(ens_hotspot1(theta,B1,f,n,j))
	if (B1==0) return(ens_hotspot1(theta,B0,1-f,n,j))
	if (B0==0 & B1==0) return(ens_neutral(theta,j))
	f * ens_constant(theta,B1,n,j) + (1-f) * ens_constant(theta,B0,n,j)
}



# EXPECTATIONS WITH ORIENTATION ERRORS ####

# Neutrality
ens_neutral_err <- function(theta1,theta2,e1,e2,n,j) {
  (1-e1)*ens_neutral(theta1,j) + e2*ens_neutral(theta2,n-j)
}

# Constant gBGC with intensity B. B can be positive (W->S mutations) or negative (S->W) mutations
ens_constant_err <- function(theta1,theta2,B,e1,e2,n,j) {
  (1-e1)*ens_constant(theta1,B,n,j) + e2*ens_constant(theta2,-B,n,n-j)
}

# Hotspot1: gBGC with intensity B in hotspot (fraction f) and 0 otherwise (fraction 1-f)
ens_hotspot1_err <- function(theta1,theta2,B,f,e1,e2,n,j) {
  (1-e1)*ens_hotspot1(theta1,B,f,n,j) + e2*ens_hotspot1(theta2,-B,f,n,n-j)
}

# Hotspot2: gBGC with intensity B1 in hotspot (fraction f) and B0 otherwise (fraction 1-f)
ens_hotspot2_err <- function(theta1,theta2,B0,B1,f,e1,e2,n,j) {
  (1-e1)*ens_hotspot2(theta1,B0,B1,f,n,j) + e2*ens_hotspot2(theta2,-B0,-B1,f,n,n-j)
}





