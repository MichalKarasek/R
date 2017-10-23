Rajon_back <- function(filepath, choice = FALSE, gon = TRUE){

  options(digits = 10)
  values = read.table(file = filepath, header = FALSE, sep=" ", dec = ".")
  dimension = dim(values)

  if(dimension[1] == 5){

    ro = 200/ pi

    y1 = values[1, 2]
    x1 = values[1, 1]
    y2 = values[2, 2]
    x2 = values[2, 1]
    ss = values[5, 1]
    sm = values[4, 1]
    sd = values[4, 2]

    if(choice == TRUE){

      Bod = rajon_back_calculation(values, gon)
      print(Bod)

      xp = Bod[1]
      yp = Bod[2]

    }else{

      xp = values[3, 1]
      yp = values[3, 2]
    }

    S1p=sqrt((x1-xp)*(x1-xp)+(y1-yp)*(y1-yp))
    S2p=sqrt((x2-xp)*(x2-xp)+(y2-yp)*(y2-yp))

    D=matrix( c(-1,0,1,0,0,1),nrow=2,ncol=3)

    A1_1=ro*((y2-yp)/(S2p^2)-(y1-yp)/S1p^2)
    A1_2=ro*(-(x2-xp)/(S2p^2)+(x1-xp)/S1p^2)
    A1_3=-(x1-xp)/S1p
    A1_4=-(y1-yp)/S1p

    A1=matrix( c(A1_1,A1_3,A1_2,A1_4),nrow=2,ncol=2)

    A2_1=ro*((y1-yp)/(S1p^2))
    A2_2=-ro*((x1-xp)/(S1p^2))
    A2_3=-ro*((y2-yp)/(S2p^2))
    A2_4=ro*((x2-xp)/(S2p^2))
    A2_5=(x1-xp)/S1p
    A2_6=(y1-yp)/S1p

    A2=matrix( c(A2_1,A2_5,A2_2,A2_6,A2_3,0,A2_4,0),nrow=2,ncol=4)

    #kovariancni matice souradnic

    K=solve(A1)%*%D
    L=solve(A1)%*%A2

    El_1=c(sm^2,sm^2,sd^2)
    El=diag(El_1,3,3)

    Ex_1=c(ss^2,ss^2,ss^2,ss^2)
    Ex=diag(Ex_1,4,4)

    # kovarianci matice vlivu mereni

    EXm=K%*%El%*%t(K)

    #chyby

    sigma_xm=sqrt(EXm[1,1])
    sigma_ym=sqrt(EXm[2,2])
    sigma_xym=sqrt(((sigma_xm^2+sigma_ym^2))/2)

    cm=sqrt((sigma_xm^2-sigma_ym^2)^2+4*EXm[1,2]^2)
    am=sqrt((sigma_xm^2+sigma_ym^2+cm)/2)
    bm=sqrt((sigma_xm^2+sigma_ym^2-cm)/2)
    alfa2m=(atan2(2*EXm[1,2],(sigma_xm^2-sigma_ym^2)))*ro
    alfam=alfa2m/2
    if(alfa2m<0){alfam=(alfa2m+400)/2
    }

    #kovariancni matice souradnic

    EX=K%*%El%*%t(K)+L%*%Ex%*%t(L)

    #chyby

    sigma_x=sqrt(EX[1,1])
    sigma_y=sqrt(EX[2,2])
    sigma_xy=sqrt(((sigma_x^2+sigma_y^2))/2)

    c=sqrt((sigma_x^2-sigma_y^2)^2+4*EX[1,2]^2)
    a=sqrt((sigma_x^2+sigma_y^2+c)/2)
    b=sqrt((sigma_x^2+sigma_y^2-c)/2)
    alfa2=(atan2(2*EX[1,2],(sigma_x^2-sigma_y^2)))*ro
    alfa=alfa2/2
    if(alfa2<0){
      alfa=(alfa2+400)/2
    }

    vector_X <- c(x1,x2,xp)
    vector_Y <- c(y1,y2,yp)
    X_max <- vector_X[1]
    X_min <- vector_X[1]
    Y_max <- vector_Y[1]
    Y_min <- vector_Y[1]

    for(i in 1:length(vector_X)){

      if(X_max < vector_X[i]){
        X_max <- vector_X[i]
      }

      if(X_min > vector_X[i]){
        X_min <- vector_X[i]
      }

      if(Y_max < vector_Y[i]){
        Y_max <- vector_Y[i]
      }

      if(Y_min > vector_Y[i]){
        Y_min <- vector_Y[i]
      }
    }

    Elipsa_mereni <- Elipsa(am, bm, xp, yp, alfam)
    Elipsa_podkladu <- Elipsa(a, b, xp, yp, alfa)

    plot(Elipsa_mereni, type = "l", main = "Elipsa chyb Rajon zpet", sub = "Zvetseno 15 000", xlab = "X [m]", ylab = "Souradnice Y [m]",
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

    write.csv(ellipsis, file = "Output_Rajon_back.csv")
    return(ellipsis)

  }else if(dimension[1] == 4){

    ro = 200/ pi

    y1 = values[1, 2]
    x1 = values[1, 1]
    y2 = values[2, 2]
    x2 = values[2, 1]
    sm = values[4, 1]
    sd = values[4, 2]

    if(choice == TRUE){

      Bod = rajon_back_calculation(values, gon)
      print(Bod)

      xp = Bod[1]
      yp = Bod[2]

    }else{

      xp = values[3, 1]
      yp = values[3, 2]
    }

    S1p=sqrt((x1-xp)*(x1-xp)+(y1-yp)*(y1-yp))
    S2p=sqrt((x2-xp)*(x2-xp)+(y2-yp)*(y2-yp))

    D=matrix( c(-1,0,1,0,0,1),nrow=2,ncol=3)

    A1_1=ro*((y2-yp)/(S2p^2)-(y1-yp)/S1p^2)
    A1_2=ro*(-(x2-xp)/(S2p^2)+(x1-xp)/S1p^2)
    A1_3=-(x1-xp)/S1p
    A1_4=-(y1-yp)/S1p

    A1=matrix( c(A1_1,A1_3,A1_2,A1_4),nrow=2,ncol=2)

    #kovariancni matice souradnic

    K=solve(A1)%*%D

    El_1=c(sm^2,sm^2,sd^2)
    El=diag(El_1,3,3)

    # kovarianci matice vlivu mereni

    EXm=K%*%El%*%t(K)

    #chyby

    sigma_xm=sqrt(EXm[1,1])
    sigma_ym=sqrt(EXm[2,2])
    sigma_xym=sqrt(((sigma_xm^2+sigma_ym^2))/2)

    cm=sqrt((sigma_xm^2-sigma_ym^2)^2+4*EXm[1,2]^2)
    am=sqrt((sigma_xm^2+sigma_ym^2+cm)/2)
    bm=sqrt((sigma_xm^2+sigma_ym^2-cm)/2)
    alfa2m=(atan2(2*EXm[1,2],(sigma_xm^2-sigma_ym^2)))*ro
    alfam=alfa2m/2
    if(alfa2m<0){alfam=(alfa2m+400)/2
    }

    vector_X <- c(x1,x2,xp)
    vector_Y <- c(y1,y2,yp)
    X_max <- vector_X[1]
    X_min <- vector_X[1]
    Y_max <- vector_Y[1]
    Y_min <- vector_Y[1]

    for(i in 1:length(vector_X)){

      if(X_max < vector_X[i]){
        X_max <- vector_X[i]
      }

      if(X_min > vector_X[i]){
        X_min <- vector_X[i]
      }

      if(Y_max < vector_Y[i]){
        Y_max <- vector_Y[i]
      }

      if(Y_min > vector_Y[i]){
        Y_min <- vector_Y[i]
      }
    }

    Elipsa_mereni <- Elipsa(am, bm, xp, yp, alfam)

    plot(Elipsa_mereni, type = "l", main = "Elipsa chyb Rajon zpet", sub = "Zvetseno 15 000", xlab = "X [m]", ylab = "Y [m]",
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

    write.csv(ellipsis, file = "Output_Rajon_back.csv")
    return(ellipsis)

  }else{

    stop("Invalid text file entry! Check help box for details.")
  }
}
