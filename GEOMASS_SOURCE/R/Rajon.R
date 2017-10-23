Rajon<-function(filepath, choice = FALSE, gon = TRUE){

options(digits = 10)
Ro = 200/pi

values = read.table(file = filepath, header = FALSE, sep = " ", dec = ".")

dimension = dim(values)

if(dimension[1] == 6){

  sigma_smer = values[4, 1]
  sigma_delka = values[4, 2]

  X_stanovisko = values[1, 1]
  Y_stanovisko = values[1, 2]

  Chyba_X_stan = values[5, 1]
  Chyba_Y_stan = values[5, 2]

  X_orientace = values[2, 1]
  Y_orientace = values[2, 2]

  Chyba_X_ori = values[6, 1]
  Chyba_Y_ori = values[6, 2]

  if(choice == TRUE){

    Bod <- rajon_calculation(values, gon)

    X_urcovany_bod = Bod[1]
    Y_urcovany_bod = Bod[2]

  }else{

    X_urcovany_bod = values[3, 1]
    Y_urcovany_bod = values[3, 2]
  }

  Vzdalenost = sqrt((X_urcovany_bod - X_stanovisko)^2 + (Y_urcovany_bod - Y_stanovisko)^2)
  Vzdalenost_Stan_Or = sqrt((X_orientace - X_stanovisko)^2 + (Y_orientace - Y_stanovisko)^2)

  D = matrix(c(-1, 0, 1, 0, 0, 1),nrow = 2,ncol = 3)

  # Matice vah

  M_vah = diag(c(sigma_smer^2, sigma_smer^2, sigma_delka^2))
  M_vah_podklad = diag(c(Chyba_X_stan^2 , Chyba_X_stan^2 , Chyba_X_ori^2 , Chyba_Y_ori^2))

  # Matice urcovane konfigurace

  A1_1 = -Ro*((Y_urcovany_bod - Y_stanovisko)/Vzdalenost^2)
  A1_2 = Ro*((X_urcovany_bod - X_stanovisko)/Vzdalenost^2)
  A2_1 = (X_urcovany_bod - X_stanovisko)/Vzdalenost
  A2_2 = (Y_urcovany_bod - Y_stanovisko)/Vzdalenost

  Matice_urcovane_konfigurace = matrix(c(A1_1,A2_1,A1_2,A2_2), nrow = 2)

  # Matice dane konfigurae
  B1_1 = Ro * (((Y_urcovany_bod - Y_stanovisko) / Vzdalenost^2) - ((Y_orientace - Y_stanovisko) / Vzdalenost_Stan_Or^2))
  B1_2 = Ro * (- ((X_urcovany_bod - X_stanovisko) / Vzdalenost^2) + ((X_orientace - X_stanovisko) / Vzdalenost_Stan_Or^2))
  B1_3 = Ro * ((Y_orientace - Y_stanovisko) / Vzdalenost_Stan_Or^2)
  B1_4 = - Ro * ((X_orientace - X_stanovisko) / Vzdalenost_Stan_Or^2)
  B2_1 = - (X_urcovany_bod - X_stanovisko) / Vzdalenost
  B2_2 = - (Y_urcovany_bod - Y_stanovisko) / Vzdalenost
  B2_3 = 0
  B2_4 = 0

  Matice_dane_konfigurace = matrix(c(B1_1, B2_1, B1_2, B2_2, B1_3, B2_3, B1_4, B2_4),nrow = 2, ncol = 4)

  # Kovariancni matice mereni

  K = solve(Matice_urcovane_konfigurace) %*% D

  # Kovariancni matice podkladu

  L = solve(Matice_urcovane_konfigurace) %*% Matice_dane_konfigurace

  # Vliv mereni

  Vliv_mer = K %*% M_vah %*% t(K)

  #chyby
  sigma_xm = sqrt(Vliv_mer[1, 1])
  sigma_ym = sqrt(Vliv_mer[2, 2])
  sigma_xym = sqrt(((sigma_xm^2 + sigma_ym^2)) / 2)

  cm = sqrt((sigma_xm^2 - sigma_ym^2)^2 + 4 * Vliv_mer[1, 2]^2)
  am = sqrt((sigma_xm^2 + sigma_ym^2 + cm) / 2)
  bm = sqrt((sigma_xm^2 + sigma_ym^2 - cm) / 2)
  alfa2m = (atan2(2 * Vliv_mer[1, 2], (sigma_xm^2 - sigma_ym^2))) * Ro
  alfam = alfa2m / 2

  if(alfa2m<0){
    alfa2m = (alfa2m + 400) / 2
  }

  # Vliv mereni a dane konfigurace

  Vliv_celkovy = Vliv_mer + L %*% M_vah_podklad %*% t(L)

  #chyby
  sigma_x = sqrt(Vliv_celkovy[1, 1])
  sigma_y = sqrt(Vliv_celkovy[2, 2])
  sigma_xy = sqrt(((sigma_x^2 + sigma_y^2)) / 2)

  c = sqrt((sigma_x^2 - sigma_y^2)^2 + 4 * Vliv_celkovy[1, 2]^2)
  a = sqrt((sigma_x^2 + sigma_y^2 + c) / 2)
  b = sqrt((sigma_x^2 + sigma_y^2 - c) / 2)
  alfa2 = (atan2(2 * Vliv_celkovy[1, 2], (sigma_x^2 - sigma_y^2))) * Ro
  alfa = alfa2 / 2

  if(alfa < 0){
    alfa = (alfa + 400) / 2
  }

  vector_X = c(X_stanovisko, X_orientace, X_urcovany_bod)
  vector_Y = c(Y_stanovisko, Y_orientace, Y_urcovany_bod)
  X_max = vector_X[1]
  X_min = vector_X[1]
  Y_max = vector_Y[1]
  Y_min = vector_Y[1]

  for(i in 1:length(vector_X)){

    if(X_max < vector_X[i]){
      X_max = vector_X[i]
    }

    if(X_min > vector_X[i]){
      X_min = vector_X[i]
    }

    if(Y_max < vector_Y[i]){
      Y_max = vector_Y[i]
    }

    if(Y_min > vector_Y[i]){
      Y_min = vector_Y[i]
    }
  }

  Elipsa_mereni = Elipsa(am, bm, X_urcovany_bod, Y_urcovany_bod, alfam)
  Elipsa_podkladu = Elipsa(a, b, X_urcovany_bod, Y_urcovany_bod, alfa)

  plot(Elipsa_mereni, type = "l", main = "Elipsa chyb Rajon", sub = "Zvetseno 15 000", xlab = "X [m]", ylab = "Y [m]",
       xlim = range(X_min-500, X_max+500), ylim = range(Y_min-500, Y_max+500), col = "green", lty = 1, asp = 1)
  points(vector_X, vector_Y, pch=21, bg="red")
  lines(Elipsa_podkladu, type = "l", col = "blue", lty = 1)
  legend("topright", title = "Error ellipsis", c("Ellipse of measurements", "Ellipse of foundation"), col = c("green", "blue"), lty = c(1,1), bg = "gray90")

  ellipsis <- data.frame(
    Determined_point = c("Impact of measurement","Impact of measurement + background"),
    Sigma_x_mm = round(c(sigma_xm, sigma_x), digits = 2),
    Sigma_y_mm = round(c(sigma_ym, sigma_y), digits = 2),
    Sigma_xy_mm = round(c(sigma_xym, sigma_xy), digits = 2),
    Main_half_axis_a_mm = round(c(am, a), digits = 2),
    Minor_half_axis_b_mm = round(c(bm, b), digits = 2),
    Rotation_angle_gon = round(c(alfam, alfa), digits = 2),
    stringsAsFactors = FALSE
  )

  write.csv(ellipsis, file = "Output_Rajon.csv")
  return(ellipsis)

}else if(dimension[1] == 4){

  sigma_smer = values[4, 1]
  sigma_delka = values[4, 2]

  X_stanovisko = values[1, 1]
  Y_stanovisko = values[1, 2]

  X_orientace = values[2, 1]
  Y_orientace = values[2, 2]

  if(choice == TRUE){

    Bod <- rajon_calculation(values)

    X_urcovany_bod = Bod[1]
    Y_urcovany_bod = Bod[2]

  }else{

    X_urcovany_bod = values[3, 1]
    Y_urcovany_bod = values[3, 2]
  }

  Vzdalenost = sqrt((X_urcovany_bod - X_stanovisko)^2 + (Y_urcovany_bod - Y_stanovisko)^2)
  Vzdalenost_Stan_Or = sqrt((X_orientace - X_stanovisko)^2 + (Y_orientace - Y_stanovisko)^2)

  D = matrix(c(-1, 0, 1, 0, 0, 1),nrow = 2,ncol = 3)

  # Matice vah

  M_vah = diag(c(sigma_smer^2, sigma_smer^2, sigma_delka^2))

  # Matice urcovane konfigurace

  A1_1 = -Ro*((Y_urcovany_bod - Y_stanovisko)/Vzdalenost^2)
  A1_2 = Ro*((X_urcovany_bod - X_stanovisko)/Vzdalenost^2)
  A2_1 = (X_urcovany_bod - X_stanovisko)/Vzdalenost
  A2_2 = (Y_urcovany_bod - Y_stanovisko)/Vzdalenost

  Matice_urcovane_konfigurace = matrix(c(A1_1, A2_1, A1_2, A2_2), nrow = 2)

  # Kovariancni matice mereni

  K = solve(Matice_urcovane_konfigurace) %*% D

  # Vliv mereni

  Vliv_mer = K %*% M_vah %*% t(K)

  #chyby
  sigma_xm = sqrt(Vliv_mer[1, 1])
  sigma_ym = sqrt(Vliv_mer[2, 2])
  sigma_xym = sqrt(((sigma_xm^2 + sigma_ym^2)) / 2)

  cm = sqrt((sigma_xm^2 - sigma_ym^2)^2 + 4 * Vliv_mer[1, 2]^2)
  am = sqrt((sigma_xm^2 + sigma_ym^2 + cm) / 2)
  bm = sqrt((sigma_xm^2 + sigma_ym^2 - cm) / 2)
  alfa2m = (atan2(2 * Vliv_mer[1, 2], (sigma_xm^2 - sigma_ym^2))) * Ro
  alfam = alfa2m / 2

  if(alfa2m < 0){
    alfa2m = (alfa2m + 400) / 2
  }

  vector_X = c(X_stanovisko, X_orientace, X_urcovany_bod)
  vector_Y = c(Y_stanovisko, Y_orientace, Y_urcovany_bod)
  X_max = vector_X[1]
  X_min = vector_X[1]
  Y_max = vector_Y[1]
  Y_min = vector_Y[1]

  for(i in 1:length(vector_X)){

    if(X_max < vector_X[i]){
      X_max = vector_X[i]
    }

    if(X_min > vector_X[i]){
      X_min = vector_X[i]
    }

    if(Y_max < vector_Y[i]){
      Y_max = vector_Y[i]
    }

    if(Y_min > vector_Y[i]){
      Y_min = vector_Y[i]
    }
  }

  Elipsa_mereni = Elipsa(am, bm, X_urcovany_bod, Y_urcovany_bod, alfam)

  plot(Elipsa_mereni, type = "l", main = "Elipsa chyb Rajon", sub = "Zvetseno 15 000", xlab = "X [m]", ylab = "Y [m]",
       xlim = range(X_min-500, X_max+500), ylim = range(Y_min-500, Y_max+500), col = "green", lty = 1, asp = 1)
  points(vector_X, vector_Y, pch=21, bg="red")
  legend("topright", title = "Error ellipsis", c("Ellipse of measurements"), col = c("green"), lty = c(1), bg = "gray90")

  ellipsis <- data.frame(
    Determined_point = c("Impact of measurement"),
    Sigma_x_mm = round(c(sigma_xm), digits = 2),
    Sigma_y_mm = round(c(sigma_ym), digits = 2),
    Sigma_xy_mm = round(c(sigma_xym), digits = 2),
    Main_half_axis_a_mm = round(c(am), digits = 2),
    Minor_half_axis_b_mm = round(c(bm), digits = 2),
    Rotation_angle_gon = round(c(alfam), digits = 2),
    stringsAsFactors = FALSE
  )

  write.csv(ellipsis, file = "Output_Rajon.csv")
  return(ellipsis)

}else{
  stop("Invalid text file entry! Check help box for details.")
}
}
