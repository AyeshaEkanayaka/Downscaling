library(ggplot2)
library(viridis)
library(ggmap)
library(sf)
library(gridExtra)
library(ncdf4)
library(graphics)
library(raster)
library(rgdal)
library(data.table)
library(grid)
library(RColorBrewer)
library(fields)

#Big cluster
lon.A.1=140
lat.A.1=-25
lon.A.2=155
lat.A.2=-10
#Exact area 
#GBRMPA polygon. See the map


GBRMPA.lon.1=142
GBRMPA.lon.2=145
GBRMPA.lon.3=145
GBRMPA.lon.4=146
GBRMPA.lon.5=147
GBRMPA.lon.6=152.55
GBRMPA.lon.7=154
GBRMPA.lon.8=142

GBRMPA.lat.1=-10.40
GBRMPA.lat.2=-10.40
GBRMPA.lat.3=-12.59
GBRMPA.lat.4=-14.59
GBRMPA.lat.5=-17.29
GBRMPA.lat.6=-20.59
GBRMPA.lat.7=-24.29
GBRMPA.lat.8=-24.29


MUR.attribs=as.matrix(fread("MUR_attribs.csv"))

Predicted.temp.fusedguess=as.matrix(fread("Predicted_temp.fusedguess.csv"))
Predicted.sd.fusedguess=as.matrix(fread("Predicted_sd.fusedguess.csv"))

Trend.future=as.matrix(fread("Trend_weighted.csv",select=c(106:185)))
Interp.GCM.future=as.matrix(fread("Interp_GCMs_flat.csv",select=c(106:185)))



dim(MUR.attribs)
dim(Predicted.temp.fusedguess)
dim(Trend.future)
dim(Interp.GCM.future)

#Temperature at random locations


r=sample(1:nrow(MUR.attribs),4)
MUR.attribs[r,3:4]
dev.off()
par(mfrow=c(2,2))
par(mar = c(2,2,2,2),oma=c(0.3,0.5,0,0.3),mai=c(0.5,0.5,0.5,0.3))



min=min(c(Trend.future[r[1],],Interp.GCM.future[r[1],],Predicted.temp.fusedguess[r[1],-(1:3)]-Predicted.sd.fusedguess[r[1],-(1:3)],Predicted.temp.fusedguess[r[1],-(1:3)]+Predicted.sd.fusedguess[r[1],-(1:3)]))
max=max(c(Trend.future[r[1],],Interp.GCM.future[r[1],],Predicted.temp.fusedguess[r[1],-(1:3)]-Predicted.sd.fusedguess[r[1],-(1:3)],Predicted.temp.fusedguess[r[1],-(1:3)]+Predicted.sd.fusedguess[r[1],-(1:3)]))
plot(2020:2099,Interp.GCM.future[r[1],],type="l",col="red",ylab="MSE",xlab="years",main="location 1",ylim=c(min,max))

lines(2020:2099,Trend.future[r[1],],type="l",col="green")

lines(2020:2099,Predicted.temp.fusedguess[r[1],-(1:3)],type="l",col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[1],-(1:3)]-Predicted.sd.fusedguess[r[1],-(1:3)],lty=2,col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[1],-(1:3)]+Predicted.sd.fusedguess[r[1],-(1:3)],lty=2,col="blue")

min=min(c(Trend.future[r[2],],Interp.GCM.future[r[2],],Predicted.temp.fusedguess[r[2],-(1:3)]-Predicted.sd.fusedguess[r[2],-(1:3)],Predicted.temp.fusedguess[r[2],-(1:3)]+Predicted.sd.fusedguess[r[2],-(1:3)]))
max=max(c(Trend.future[r[2],],Interp.GCM.future[r[2],],Predicted.temp.fusedguess[r[2],-(1:3)]-Predicted.sd.fusedguess[r[2],-(1:3)],Predicted.temp.fusedguess[r[2],-(1:3)]+Predicted.sd.fusedguess[r[2],-(1:3)]))
plot(2020:2099,Interp.GCM.future[r[2],],type="l",col="red",ylab="MSE",xlab="years",main="location 2",ylim=c(min,max))

lines(2020:2099,Trend.future[r[2],],type="l",col="green")

lines(2020:2099,Predicted.temp.fusedguess[r[2],-(1:3)],type="l",col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[2],-(1:3)]-Predicted.sd.fusedguess[r[2],-(1:3)],lty=2,col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[2],-(1:3)]+Predicted.sd.fusedguess[r[2],-(1:3)],lty=2,col="blue")

min=min(c(Trend.future[r[3],],Interp.GCM.future[r[3],],Predicted.temp.fusedguess[r[3],-(1:3)]-Predicted.sd.fusedguess[r[3],-(1:3)],Predicted.temp.fusedguess[r[3],-(1:3)]+Predicted.sd.fusedguess[r[3],-(1:3)]))
max=max(c(Trend.future[r[3],],Interp.GCM.future[r[3],],Predicted.temp.fusedguess[r[3],-(1:3)]-Predicted.sd.fusedguess[r[3],-(1:3)],Predicted.temp.fusedguess[r[3],-(1:3)]+Predicted.sd.fusedguess[r[3],-(1:3)]))
plot(2020:2099,Interp.GCM.future[r[3],],type="l",col="red",ylab="MSE",xlab="years",main="location 3",ylim=c(min,max))

lines(2020:2099,Trend.future[r[3],],type="l",col="green")

lines(2020:2099,Predicted.temp.fusedguess[r[3],-(1:3)],type="l",col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[3],-(1:3)]-Predicted.sd.fusedguess[r[3],-(1:3)],lty=2,col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[3],-(1:3)]+Predicted.sd.fusedguess[r[3],-(1:3)],lty=2,col="blue")

min=min(c(Trend.future[r[4],],Interp.GCM.future[r[4],],Predicted.temp.fusedguess[r[4],-(1:3)]-Predicted.sd.fusedguess[r[4],-(1:3)],Predicted.temp.fusedguess[r[4],-(1:3)]+Predicted.sd.fusedguess[r[4],-(1:3)]))
max=max(c(Trend.future[r[4],],Interp.GCM.future[r[4],],Predicted.temp.fusedguess[r[4],-(1:3)]-Predicted.sd.fusedguess[r[4],-(1:3)],Predicted.temp.fusedguess[r[4],-(1:3)]+Predicted.sd.fusedguess[r[4],-(1:3)]))
plot(2020:2099,Interp.GCM.future[r[4],],type="l",col="red",ylab="MSE",xlab="years",main="location 4",ylim=c(min,max))

lines(2020:2099,Trend.future[r[4],],type="l",col="green")

lines(2020:2099,Predicted.temp.fusedguess[r[4],-(1:3)],type="l",col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[4],-(1:3)]-Predicted.sd.fusedguess[r[4],-(1:3)],lty=2,col="blue")
lines(2020:2099,Predicted.temp.fusedguess[r[4],-(1:3)]+Predicted.sd.fusedguess[r[4],-(1:3)],lty=2,col="blue")



#MSE

MSE.predicted.fusedguess=rep(NA,ncol(Interp.GCM.future))

MSE.Trend=rep(NA,ncol(Interp.GCM.future))


for(i in 1:ncol(Interp.GCM.future)){
  MSE.predicted.fusedguess[i]= mean(na.omit((Interp.GCM.future[,i]-Predicted.temp.fusedguess[,3+i])^2))
  MSE.Trend[i]= mean(na.omit((Interp.GCM.future[,i]-Trend.future[,i])^2))
}

dev.off()
min=min(MSE.Trend,MSE.predicted.fusedguess)
max=max(MSE.Trend,MSE.predicted.fusedguess)

plot(2020:2099,MSE.Trend,type="l",col="green",ylim=c(min,max),ylab="MSE",xlab="years",main="MSE values
(February)")
lines(2020:2099,MSE.predicted.fusedguess,type="l",col="blue")

legend(2020,0.78,legend=c("Trend","BGL"),fill=c("green","blue"))

