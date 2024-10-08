# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

# Files with constants and default parameters

# Constants used when functions are not defined in 0
ZERO <- 10^(-8)
ONE <- 1 - ZERO

# Default parameters for optimization
# gBGC coefficient
BMIN <- -100
BMAX <- 100
# Log of the mutation bias
MMIN <- -5
MMAX <- 5
# Option parameters of the optim function
MAXIT <- 100
FACTR <- 10^7
LMM <- 20
VERBOSE <- 0
USEGR <- T
