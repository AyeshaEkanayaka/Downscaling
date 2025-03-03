library(ggplot2)
library(gridExtra)
library(raster)
library(rgdal)
library(data.table)
library(RColorBrewer)

path_matlabRuntime="/Applications/MATLAB/MATLAB_Runtime/R2024a" #Where you matlab runtime is installed
path="/Volumes/Disk2/laptop-Backup/JRSSC/Git/"#Path to your working directory
Freq=12#Frequency of your data. For example if its is monthly data,frequency is 12. If it is daily data, frequency is 365.

source(paste0(path,"/src/Downscaling.R"))
Downscaling(path_matlabRuntime,path,Freq)

#Make a simple visual comparison between model SSTs and downscaled SSTs using ggplot
Model=as.matrix(fread(paste0(path,"/Model_data.csv")))
Obs=as.matrix(fread(paste0(path,"/Obs_data.csv")))
Downscaled=as.matrix(fread(paste0(path,"/BGL_Downscaled.csv")))


Model_fututure=Model[which(Model[,1]>=min(Downscaled[,1])-2 & max(Downscaled[,1])+2>=Model[,1] & Model[,2]>=min(Downscaled[,2])-2 & max(Downscaled[,2])+2>=Model[,2]),c(1:2,183:218)]
min=min(c(Downscaled[,3],Model_fututure[,3]))
max=max(c(Downscaled[,3],Model_fututure[,3]))

Downscaled=as.data.frame(Downscaled[,c(1:3)])
colnames(Downscaled)=c("lon","lat","SST")
Model_fututure=as.data.frame(Model_fututure[,1:3])
colnames(Model_fututure)=c("lon","lat","SST")


Temp.col=colorRampPalette(c("darkblue","yellow","red"))
col=Temp.col(1000)

# Big cluster
lon.A.1 <- min(Obs[,1])
lat.A.1 <- min(Obs[,2])
lon.A.2 <- max(Obs[,1])
lat.A.2 <- max(Obs[,2])

# Coastline .shp
coastlines <- readOGR(paste0(path, "ne-coastlines-10m/ne_10m_coastline.shp"))

crop_area.A <- extent(lon.A.1, lon.A.2, lat.A.1, lat.A.2)
coastlines_A <- crop(x = coastlines, y = crop_area.A)

p1=ggplot()+
  geom_path(data=coastlines_A,color="black",inherit.aes = FALSE,aes(x = long,y = lat,group=group))+
  geom_raster(data = Model_fututure , aes(x = lon, y = lat, fill = SST))+
  scale_fill_gradientn(colours = col)+
  labs(title="Model SSTs")+
  theme(text = element_text(size = 8))+
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(hjust = 0.5))

  
p2=ggplot()+
  geom_path(data=coastlines_A,color="black",inherit.aes = FALSE,aes(x = long,y = lat,group=group))+
  geom_raster(data = Downscaled , aes(x = lon, y = lat, fill = SST))+
  scale_fill_gradientn(colours = col)+
  labs(title="Downscaled SSTs")+
  theme(text = element_text(size = 8))+
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(hjust = 0.5))

grid.arrange(p1,p2, ncol = 2)
