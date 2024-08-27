# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

source("R/bgc_optimization.R")

sfsWS <- c(4000,2000,1000,500,300,200,80,50,30,20,10)
sfsSW <- c(10000,4500,2000,800,400,200,50,20,8,6,7)
e1 <- 0.01
e2 <- 0.03
sfsWSe <- sfsWS*(1-e1) + e2*rev(sfsSW)
sfsSWe <- sfsSW*(1-e2) + e1*rev(sfsWS)

plot(sfsWS,log="y")
plot(sfsSW,log="y")

plot(sfsWSe,log="y")
plot(sfsSWe,log="y")

plot(log(sfsWS/sfsSW))
points(log(sfsWSe/sfsSWe),col="red")

sum(sfsWSe + sfsSWe)
sum(sfsWS + sfsSW)
sum(sfsWSc + sfsSWc)

sfsWSc <- ((1-e2)*sfsWSe-e2*rev(sfsSWe))/(1-e1-e2)
sfsSWc <- ((1-e1)*sfsSWe-e1*rev(sfsWSe))/(1-e1-e2)


plot(log(sfsWS/sfsSW))

least_square_M(sfsWS,sfsSW,0.65)
least_square_B(sfsWS,sfsSW,0.65)
least_square_BM(sfsWS,sfsSW,0.65)
least_square_hotspot1(sfsWS,sfsSW,0.65)
least_square_hotspot2(sfsWS,sfsSW,0.65)
least_square_hotspot2bis(sfsWS,sfsSW,0.65,0.2)

least_square_hotspot1(sfsWS,sfsSW,0.65,Usegr = F)$criteria$SSres
least_square_hotspot1(sfsWS,sfsSW,0.65)$criteria$SSres

least_square_hotspot2bis(sfsWS,sfsSW,0.65,0.2,Usegr = F)
least_square_hotspot2bis(sfsWS,sfsSW,0.65,0.28)$param$B1

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

d_expected_log_ratio(sfsWS,sfsSW,0.0001,0.0001)


