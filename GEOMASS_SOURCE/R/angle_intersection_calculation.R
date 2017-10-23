angle_intersection_calculation <- function(values, gon, choice = FALSE){

  options(digits = 10)

  if(choice == TRUE){
    values = read.table(file = values, header = FALSE, sep = " ", dec = ".")
  }

  X_A = values[1,1]
  Y_A = values[1,2]

  X_B = values[2,1]
  Y_B = values[2,2]

  omega_A = values[3,1]
  omega_B = values[3,2]

  if(gon == TRUE){

    X_P = X_A + ((X_B - X_A) * (1 / tan(omega_A*pi/200)) - (Y_B - Y_A)) / ((1 / tan(omega_A*pi/200)) + (1 / tan(omega_B*pi/200)))
    Y_P = Y_A + ((Y_B - Y_A) * (1 / tan(omega_A*pi/200)) + (X_B - X_A)) / ((1 / tan(omega_A*pi/200)) + (1 / tan(omega_B*pi/200)))

  }else{

    X_P = X_A + ((X_B - X_A) * (1 / tan(omega_A*pi/180)) - (Y_B - Y_A)) / ((1 / tan(omega_A*pi/180)) + (1 / tan(omega_B*pi/180)))
    Y_P = Y_A + ((Y_B - Y_A) * (1 / tan(omega_A*pi/180)) + (X_B - X_A)) / ((1 / tan(omega_A*pi/180)) + (1 / tan(omega_B*pi/180)))

  }

  result = matrix(as.numeric(c(X_P, Y_P)), ncol = 2)

  return(result)
}
