Elipsa <- function(a, b, Xp, Yp, om){
  
  if(a < b){
    
    stop("Main half-axis smaller then minor half-axis")
  }

vec <- matrix(vector(mode = "numeric"), ncol = 2)
vec_local <- vector(mode = "numeric", 2)

i <- 0

for(n in 0:4000){
  
  if(n <= 1000){
    
    X_pomocne <- sqrt((a^2 * b^2) / (b^2 + a^2 * tan(i*pi/200)^2))
    Y_pomocne <- sqrt((a^2 * b^2) / (a^2 + b^2 * (1 / tan(i*pi/200)^2)))
    
    X <- Xp + (X_pomocne*cos(om*pi/200) - Y_pomocne*sin(om*pi/200))*15
    Y <- Yp + (X_pomocne*sin(om*pi/200) + Y_pomocne*cos(om*pi/200))*15
    
  }
  
  if( (1000 < n) && (n <= 2000) ){
    
    X_pomocne <- - sqrt((a^2 * b^2) / (b^2 + a^2 * tan(i*pi/200)^2))
    Y_pomocne <- sqrt((a^2 * b^2) / (a^2 + b^2 * (1 / tan(i*pi/200)^2)))
    
    X <- Xp + (X_pomocne*cos(om*pi/200) - Y_pomocne*sin(om*pi/200))*15
    Y <- Yp + (X_pomocne*sin(om*pi/200) + Y_pomocne*cos(om*pi/200))*15
    
  }
  
  if( (2000 < n) && (n <= 3000) ){
    
    X_pomocne <- - sqrt((a^2 * b^2) / (b^2 + a^2 * tan(i*pi/200)^2))
    Y_pomocne <- - sqrt((a^2 * b^2) / (a^2 + b^2 * (1 / tan(i*pi/200)^2)))
    
    X <- Xp + (X_pomocne*cos(om*pi/200) - Y_pomocne*sin(om*pi/200))*15
    Y <- Yp + (X_pomocne*sin(om*pi/200) + Y_pomocne*cos(om*pi/200))*15
    
  }
  
  if( (3000 < n) && (n<= 4000) ){
    
    X_pomocne <- sqrt((a^2 * b^2) / (b^2 + a^2 * tan(i*pi/200)^2))
    Y_pomocne <- - sqrt((a^2 * b^2) / (a^2 + b^2 * (1 / tan(i*pi/200)^2)))
    
    X <- Xp + (X_pomocne*cos(om*pi/200) - Y_pomocne*sin(om*pi/200))*15
    Y <- Yp + (X_pomocne*sin(om*pi/200) + Y_pomocne*cos(om*pi/200))*15
  
  }
  
  vec_local[1] <- X
  vec_local[2] <- Y
  
  vec <- rbind(vec, vec_local)
  
  i <- i + 0.1
  
}

return(vec)

} 

