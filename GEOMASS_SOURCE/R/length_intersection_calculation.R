length_intersection_calculation <- function(values, choice = FALSE){

  options(digits = 10)

  if(choice == TRUE){
    values = read.table(file = values, header = FALSE, sep = " ", dec = ".")
  }

  y1 = values[1, 2]
  x1 = values[1, 1]
  y2 = values[2, 2]
  x2 = values[2, 1]
  S1_p = values[3, 1]
  S2_p = values[3, 2]

  #vypocet stanoviska pomoci transformace

  S12 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))
  u = (S12 * S12 + S1_p * S1_p - S2_p * S2_p) / (2 * S12)

  k = sqrt(S1_p * S1_p - u * u)

  yp = y1 + k * ((x2 - x1) / S12) + u * ((y2 - y1) / S12)
  xp = x1 + u * ((x2 - x1) / S12) - k * ((y2 - y1) / S12)

  Result = matrix(as.numeric(c(xp, yp)))

  return(Result)
}
