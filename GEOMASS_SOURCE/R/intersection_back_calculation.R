intersection_back_calculation <- function(values, gon, choice = FALSE){

  options(digits = 10)

  if(choice == TRUE){
    values = read.table(file = values, header = FALSE, sep = " ", dec = ".")
  }

  y1 = values[1, 2]
  x1 = values[1, 1]
  y2 = values[2, 2]
  x2 = values[2, 1]
  y3 = values[3, 2]
  x3 = values[3, 1]
  wa_ = values[4, 1]
  wb_ = values[4, 2]

  if(gon == TRUE){

    wa = wa_* (pi / 200)
    wb = wb_* (pi / 200)
  }else{

    wa = wa_* (pi / 180)
    wb = wb_* (pi / 180)
  }

  yp_y2 = - (y2 - y1) + (x2 - x1) * (1 / tan(wa))
  xp_x2 = - (x2 - x1) - (y2 - y1) * (1 / tan(wa))

  p = yp_y2 - (y3 - y2) - (x3 - x2) * (1 / tan(wb))
  q = xp_x2 - (x3 - x2) + (y3 - y2) * (1 / tan(wb))

  l = - p / q
  k = - q / p

  J = k + l
  m = (yp_y2 + l * xp_x2) / J

  yp = y2 + k * m
  xp = x2 + m

  result = matrix(as.numeric(c(xp, yp)), ncol = 2)

  return(result)
}
