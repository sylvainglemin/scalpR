source("R/bgc_optimization.R")

sfsWS <- c(1000,500,300,200,80,50)
sfsSW <- c(2000,800,400,150,50,10)
LS <- least_square_BM(sfsWS,sfsSW,0.65,Usegr = F)
c(LS$B,LS$mutbias,LS$R2)
LS <- least_square_B(sfsWS,sfsSW,0.65,Usegr = T)
c(LS$B,LS$SSres)

least_square_M(sfsWS,sfsSW,0.65)
least_square_B(sfsWS,sfsSW,0.65)
least_square_BM(sfsWS,sfsSW,0.65)

total_ss(sfsWS,sfsSW)

system.time(
  for(i in 1:1000) least_square(sfsWS,sfsSW,0.6,USEGR = F)
  )
system.time(
  for(i in 1:1000) least_square(sfsWS,sfsSW,0.6,USEGR = T)
)
