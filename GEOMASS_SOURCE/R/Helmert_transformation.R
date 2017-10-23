Helmert_transformation<-function(filepath){
options(digits=10)
cor = read.table(file = filepath, header=TRUE, sep=" ", dec = ".")

rz=ncol(cor);
for (i in 1:rz){
if(cor[1,i]=="stop"){
k=i-2}
}
b1Y=matrix(nrow=1,ncol=k/2);
b1X=matrix(nrow=1,ncol=k/2);
b1y=matrix(nrow=1,ncol=k/2);
b1x=matrix(nrow=1,ncol=k/2);


for(i in 1:k/2){
b1Y[i]=cor[1,i+1]
b1X[i]=cor[1,(k/2+i+1)]
b1y[i]=cor[1,(k+2+i)]
b1x[i]=cor[1,(3*k/2+2+i)]
}
Yt=mean(b1Y)
Xt=mean(b1X)
yt=mean(b1y)
xt=mean(b1x)

yri=matrix(nrow=1,ncol=k/2);
xri=matrix(nrow=1,ncol=k/2);
Yri=matrix(nrow=1,ncol=k/2);
Xri=matrix(nrow=1,ncol=k/2);


for(i in 1:k/2){
yri[i]=b1y[i]-yt;
xri[i]=b1x[i]-xt;
Yri[i]=b1Y[i]-Yt;
Xri[i]=b1X[i]-Xt;
}

#melo by se blizit nule
syri=sum(yri);
sxri=sum(xri);
sYri=sum(Yri);
sYri=sum(Xri);

syY=sum(yri*Yri);
sxX=sum(xri*Xri);
sxY=sum(xri*Yri);
syX=sum(yri*Xri);


rx=matrix(nrow=1,ncol=k/2);
ry=matrix(nrow=1,ncol=k/2);
r=matrix(nrow=1,ncol=k/2);

for(i in 1:(k/2)){
rx[i]=xri[i]*xri[i];
ry[i]=yri[i]*yri[i];
r[i]=xri[i]*xri[i]+yri[i]*yri[i];
}

sx2y2=sum(xri*xri+yri*yri);
sr=sum(r);
#kontrola sx2y2-sum(r) is zero

lam1=(sxX+syY)/sx2y2;
lam2=(sxY-syX)/sx2y2;
q=sqrt(lam1*lam1+lam2*lam2);
om_rad=atan2(lam2,lam1);

if(om_rad<0){
om_rad=(om_rad+2*pi);
}

om_gon=om_rad*200/pi;

Yp=Yt-lam1*yt-lam2*xt
Xp=Xt-lam1*xt+lam2*yt

Y_id=matrix(nrow=1,ncol=k/2);
X_id=matrix(nrow=1,ncol=k/2);


for(i in 1:k/2){
Y_id[i]=Yp+lam1*b1y[i]+lam2*b1x[i];
X_id[i]=Xp+lam1*b1x[i]-lam2*b1y[i];
}

#POSUNY V MILIMETRECH!!!!

vy=(Y_id-b1Y)*1000;
vx=(X_id-b1X)*1000;

svy=sum(vy);
svx=sum(vx);

vv=vy*vy+vx*vx;
svv=sum(vv);
row=sqrt(vv);

#mira identity

s_v=sqrt(sum(row*row)/(k/2));

#stredni rozdil souradnicovych chyb
s_o=sqrt(sum(row*row)/(2*(k/2)-4));

#pro minimalni pocet (4) bodu je s_v=s_o

#stredni rozdil v poloze danych bodu
s_p=s_o*sqrt(2);

#stredni chyby (smer. odchylky) neznamych, s_x=s_y, proto s_y nikde neuvadim!

s_x=matrix(nrow=1,ncol=k/2);
for(i in 1:(k/2)){
s_x[i]=sqrt((s_o/1000)^2*(1/(k/2)+r[i]/sr))
}		#ACHTUNG v Metrech!!!!

s_x_mm=s_x*1000;

#stredni rozdil v poloze transformacniho bodu
sd=sqrt(s_x^2+s_x^2);
sd_mm=sd*1000;

#PODROBNE BODY

j=rz-2*k-3;
b2y=matrix(nrow=1,ncol=j/2);
b2x=matrix(nrow=1,ncol=j/2);

for(i in 1:j/2){
b2y[i]=cor[1,(2*k+i+3)]
b2x[i]=cor[1,(j/2+2*k+i+2)]
}

Y_p=matrix(nrow=1,ncol=j/2);
X_p=matrix(nrow=1,ncol=j/2);


for(i in 1:j/2){
Y_p[i]=Yp+lam1*b2y[i]+lam2*b2x[i];
X_p[i]=Xp+lam1*b2x[i]-lam2*b2y[i];
}


identic <- data.frame(
      Transformation_points = ("identical point"),
      original_Y = round(c(b1Y), digits = 3),
      original_X = round(c(b1X), digits = 3),
 	original_y = round(c(b1y), digits = 3),
      original_x = round(c(b1x), digits = 3),
	shift_vy_mm = round(c(vy), digits = 1),
      shift_vx_mm = round(c(vx), digits = 1),
	rho_mm = round(c(row), digits = 1),
	standard_deviation_mm = round(c(sd_mm), digits = 1),
	 stringsAsFactors = TRUE )



podrob <- data.frame(
      Determined_points = c("measured point"),
          new_Y = round(c(b2y),digits = 3),
          new_X = round(c(b2y), digits = 3),
	original_y = round(c(Y_p), digits = 3),
	original_x = round(c(X_p), digits = 3),
      stringsAsFactors = FALSE
    )

infos <- data.frame(
key_lambda1 = round(c(lam1),digits = 3),
key_lambda2 = round(c(lam2), digits = 3),
	ratio = round(c(q), digits = 8),
	omega_gon = round(c(om_gon), digits = 4),
rate_of_identity_mm = round(c(s_v), digits = 2),
mean_absolute_difference_of_coordinate_deviation_mm = round(c(s_p), digits = 3),
stringsAsFactors = FALSE)



write.csv(identic,file= "Output_Helmert_transformation.csv", row.names = FALSE, sep = "\t",
 append = TRUE)
write.table(podrob,file= "Output_Helmert_transformation.csv", row.names = FALSE, sep = "\t",
col.names = TRUE,
 append = TRUE)

write.table(infos,file= "Output_Helmert_transformation.csv", row.names = FALSE, sep = "\t",
col.names = TRUE,
 append = TRUE)

}

