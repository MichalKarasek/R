rajon_calculation <- function(values, gon, choice = FALSE){

  if(choice == TRUE){
    values = read.table(file = values, header = FALSE, sep = " ", dec = ".")
  }

  X_stanovisko = values[1, 1]
  Y_stanovisko = values[1, 2]

  X_orientace = values[2, 1]
  Y_orientace = values[2, 2]

  omega = values[3, 1]
  distance = values[3, 2]

  smernik = atan2(Y_orientace - Y_stanovisko, X_orientace - X_stanovisko)

  if(gon == TRUE){

    yp = Y_stanovisko + sin(smernik + omega*pi/200) * distance
    xp = X_stanovisko + cos(smernik + omega*pi/200) * distance
  }else{

    yp = Y_stanovisko + sin(smernik + omega*pi/180) * distance
    xp = X_stanovisko + cos(smernik + omega*pi/180) * distance
  }

  Bod = matrix(as.numeric(c(xp, yp)), ncol = 2)

  return(Bod)

}
