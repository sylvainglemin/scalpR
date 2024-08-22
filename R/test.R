# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

source("R/bgc_optimization.R")

sfsWS <- c(4000,2000,1000,500,300,200,80,50,20,20,30)
sfsSW <- c(10000,4500,2000,800,400,150,50,10,8,6,10)

least_square_M(sfsWS,sfsSW,0.65)
least_square_B(sfsWS,sfsSW,0.65)
least_square_BM(sfsWS,sfsSW,0.65)
least_square_hotspot1(sfsWS,sfsSW,0.65)
least_square_hotspot2(sfsWS,sfsSW,0.65)
least_square_hotspot2bis(sfsWS,sfsSW,0.65,0.1)

least_square_hotspot1(sfsWS,sfsSW,0.65,Usegr = F)
least_square_hotspot1(sfsWS,sfsSW,0.65)

least_square_hotspot2bis(sfsWS,sfsSW,0.65,0.1,Usegr = F)
least_square_hotspot2bis(sfsWS,sfsSW,0.65,0.1)

system.time(
  for(i in 1:20) least_square_hotspot1(sfsWS,sfsSW,0.2,Usegr = F)
)
system.time(
  for(i in 1:20) least_square_hotspot1(sfsWS,sfsSW,0.2,Usegr = T)
)

system.time(
  for(i in 1:20) least_square_hotspot2(sfsWS,sfsSW,0.2,Usegr = F)
)
system.time(
  for(i in 1:20) least_square_hotspot2(sfsWS,sfsSW,0.2,Usegr = T)
)


system.time(
  for(i in 1:20) least_square_hotspot2bis(sfsWS,sfsSW,0.2,0.1,Usegr = F)
  )
system.time(
  for(i in 1:20) least_square_hotspot2bis(sfsWS,sfsSW,0.2,0.1,Usegr = T)
)


par <- c(0.3,-1.0,1.2,0.0001,0.0002)
WS <- sfsWS
SW <- sfsSW

sum_of_squares_hotspot2bis(par,WS,SW,0.7,0.61)


