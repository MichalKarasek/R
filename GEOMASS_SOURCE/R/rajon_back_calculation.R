rajon_back_calculation <- function(values, gon, choice = FALSE){

  options(digits = 10)

  if(choice == TRUE){
    values = read.table(file = values, header = FALSE, sep = " ", dec = ".")
  }

  y1 = values[1, 2]
  x1 = values[1, 1]
  y2 = values[2, 2]
  x2 = values[2, 1]
  w_rad = values[3, 1]
  S1_p = values[3, 2]

  if(gon == TRUE){

    w = w_rad * (pi / 200)
  }else{

    w= w_rad * (pi / 180)
  }

  smer = atan2((y2 - y1), (x2 - x1))

  if(smer < 0){
    alfa = (smer + 400)
  }

  S12=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))

  # poznamka wb, mozny problem sinusu 1 a 2 kvadrant

  asinwb = (S1_p/ S12) * sin(w)
  wb = asin(asinwb)

  wa = pi - (wb + w)

  smer_ap = smer + wa

  yp = y1 + S1_p * sin(smer_ap)
  xp = x1 + S1_p * cos(smer_ap)

  result = matrix(as.numeric(c(xp, yp)), ncol = 2)

  return(result)
}
