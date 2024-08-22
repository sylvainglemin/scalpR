# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr


# Derivatives of several sub-functions used in the computation of gradients

# COMMENTS OF FUNCTIONS TO ADD

d_expected_log_ratio <- function(WS,SW,e1,e2){
  d1 <- (SW + rev(WS))/((1 - e1)*SW - e1*rev(WS)) + # Log term
    1/(2*((1 - e2)*WS - e2*rev(SW))) + # First ratio term
    ((1 - e2)*rev(WS) - e2*SW)/(2*((1 - e1)*SW - e1*rev(WS))^2) # Second ratio term
  d2 <- (rev(SW) + WS)/(e2*rev(SW) - (1 - e2)*WS) + # Log term
    (e1*WS - (1 - e1)*rev(SW))/(2*((1 - e2)*WS - e2*rev(SW))^2) + # First ratio term
    1/(2*(e1*rev(WS) - (1 - e1)*SW))
  return(list("d1"=d1,"d2"=d2))
}

ratio_B <- function(B,x){
  if(B==0) return(1 - x)
  return(
    (exp(B) - exp(B*x))/((exp(B) - 1))
    )
}
ratio_B <- Vectorize(ratio_B)

d_ratio_B <- function(B,x){
  if(B==0) return(
    x*(1 - x)/2
    )
  return(
    (-exp(B) + exp(B*x)*(exp(B)* (1 - x) + x))/(-1 + exp(B))^2
  )
}
d_ratio_B <- Vectorize(d_ratio_B)

d_hotspot1 <- function(B,f,x) {
  dB <- -f*(
    d_ratio_B(B,x)/((1-f)*(1-x) + f*ratio_B(B,x)) +
      d_ratio_B(-B,x)/((1-f)*(1-x) + f*ratio_B(-B,x))
  )
  df <- ((1-x)*ratio_B(-B,x)-(1-x)*ratio_B(B,x))/
    ( ((1-f)*(1-x) + f*ratio_B(B,x))*
        ((1-f)*(1-x) + f*ratio_B(-B,x)) )
  return(list("dB"=dB,"df"=df))
}

d_hotspot2 <- function(B0,B1,f,x) {
  dB0 <- -(1-f)*(
    d_ratio_B(B0,x)/((1-f)*ratio_B(B0,x) + f*ratio_B(B1,x)) +
    d_ratio_B(-B0,x)/((1-f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  )
  dB1 <- -f*(
    d_ratio_B(B1,x)/((1-f)*ratio_B(B0,x) + f*ratio_B(B1,x)) +
    d_ratio_B(-B1,x)/((1-f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  )
  df <- (ratio_B(B0,x)*ratio_B(-B1,x)-ratio_B(-B0,x)*ratio_B(B1,x))/
    ( ((1-f)*ratio_B(B0,x) + f*ratio_B(B1,x))*
    ((1-f)*ratio_B(-B0,x) + f*ratio_B(-B1,x)) )
  return(list("dB0"=dB0,"dB1"=dB1,"df"=df))
}
