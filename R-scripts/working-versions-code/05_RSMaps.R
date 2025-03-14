#This code is intended to produce maps Rate Summation project of stephensi models
#Sadie J. Ryan
#Sept 2024 

#SET YOUR DIRECTORY
setwd("R:/Ryan_Lab/Sadie/RSPaperMapping")

#Many R packages for spatial data deprecated during the course of this project, so code was adapted as best I could along the way.

library(sf)
library(raster)
library(terra)

#SET A TEMP FOLDER THAT CAN HANDLE CAST OFF RASTERS
rasterOptions(tmpdir='R:/XXXXXXXX/DUMPSTER/') 


##########################################################################################
##CURRENT temp data (WorldClim)
#Using version 1.4 - started in 2023, updates will appear

#rgeos package was deprecated during summer 2024, so you will need to download your worldclim data to a local folder to use it 
setwd("R:/Ryan_Lab/Sadie/RSPaperMapping/wc5")
y = stack(list.files(pattern = "bil", full.names = TRUE))
yy<-y*0.1

setwd("R:/Ryan_Lab/Sadie/RSPaperMapping")

##############################################################
#CLIMATE LAYER CRUNCHER MACHINE HERE
##############################################################
xx<-yy

#Simply excluding the limits on thermal curves

#Main figure panels A-D
#####Figure Set 1 - S(T) > 0 using lower CIs - panels A-D

#A - constant
# 15.8 - 35.8

#B empDTR12
# 15.8 - 33.1

#C rsAllTraitsDTR12
# 10.9 - 40.1

#D rsSuitDTR12
# 9.3 - 41.1

#####Figure Set 2 - S(T) > 0.25 using lower CIs - panels A-D

#A - constant
# 19.2 - 33.6

#B empDTR12
# 18.5 - 31.6

#C rsAllTraitsDTR12
# 19.2 - 33.3

#D rsSuitDTR12
# 15.7 - 37.0

#####Figure Set 3 - S(T) > 0.5 using lower CIs - panels A-D

#A - constant
# 21.1 - 31.9

#B empDTR12
# 20.3 - 30.2

#C rsAllTraitsDTR12
#  21.2 - 31.1

#D rsSuitDTR12
# 18.7 - 34.7



#A - constant
# 21.1 - 31.9
a<-xx
a[a<21.1]<-NA
a[a>31.9]<-NA


#B empDTR12
# 20.3 - 30.2
b<-xx
b[b<20.3]<-NA
b[b>30.2]<-NA

#C rsAllTraitsDTR12 
## 21.2 - 31.1
c<-xx
c[c<21.2]<-NA
c[c>31.1]<-NA

#D rsSuitDTR12 
## 18.7 - 34.7
d<-xx
d[d<18.7]<-NA
d[d>34.7]<-NA


#Turning it into 0,1s 
aa<-a
aa[aa>0]<-1

bb<-b
bb[bb>0]<-1

cc<-c
cc[cc>0]<-1

dd<-d
dd[dd>0]<-1

#Save the bricks
writeRaster(aa, filename="constant_May50.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(bb, filename="empDTR12_May50.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(cc, filename="rsAllTraitsDTR12_May50.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(dd, filename="rsSuitDTR12_May50.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

#To read these back in as stacks, use XX<-stack("TMean_LT_100.tif")
#The layers will have the filename in them, but not the original layer names. It still gets you to 12 months in order.

#Adding up months in the year for persistence
sum_aa<-sum(aa, na.rm=TRUE)
sum_bb<-sum(bb, na.rm=TRUE)
sum_cc<-sum(cc, na.rm=TRUE)
sum_dd<-sum(dd, na.rm=TRUE)

#Save the sums
writeRaster(sum_aa, filename="A_constant_sum_may50.tif", format="GTiff", overwrite=TRUE)
writeRaster(sum_bb, filename="B_empDTR12_sum_may50.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(sum_cc, filename="C_rsAllTraitsDTR12_sum_may50.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(sum_dd, filename="D_rsSuitDTR12_sum_may50.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

################################################################
#Create masks of Endemic Range Countries of Interest (COI), and for Africa
EndCOI<-read_sf("RSMAPPING/AS_COI_EndRange_OCT23a.shp")
Africa<-read_sf("RSMAPPING/Africa_OCT23.shp")

A_end<-mask(sum_aa, EndCOI)
A_Africa<-mask(sum_aa, Africa)

B_end<-mask(sum_bb, EndCOI)
B_Africa<-mask(sum_bb, Africa)

C_end<-mask(sum_cc, EndCOI)
C_Africa<-mask(sum_cc, Africa)

D_end<-mask(sum_dd, EndCOI)
D_Africa<-mask(sum_dd, Africa)

writeRaster(A_end, "A_end_may50.tif", format="GTiff", overwrite=TRUE)
writeRaster(A_Africa, "A_Africa_may50.tif", format="GTiff", overwrite=TRUE)

writeRaster(B_end, "B_end_may50.tif", format="GTiff", overwrite=TRUE)
writeRaster(B_Africa, "B_Africa_may50.tif", format="GTiff", overwrite=TRUE)

writeRaster(C_end, "C_end_may50.tif", format="GTiff", overwrite=TRUE)
writeRaster(C_Africa, "C_Africa_may50.tif", format="GTiff", overwrite=TRUE)

writeRaster(D_end, "D_end_may50.tif", format="GTiff", overwrite=TRUE)
writeRaster(D_Africa, "D_Africa_may50.tif", format="GTiff", overwrite=TRUE)

################################################################
#Messing around with some plotting in R

#################################################################
library(ggplot2)
#read in the 8 tifs
A_end<-rast("A_end_may50.tif")
B_end<-rast("B_end_may50.tif")
C_end<-rast("C_end_may50.tif")
D_end<-rast("D_end_may50.tif")

A_Af<-rast("A_Africa_may50.tif")
B_Af<-rast("B_Africa_may50.tif")
C_Af<-rast("C_Africa_may50.tif")
D_Af<-rast("D_Africa_may50.tif")

### Gotta change them all to df to plot
##A##
A_end_df<-as.data.frame(A_end, xy=TRUE)
A_end_df<-na.omit(A_end_df)
colnames(A_end_df) <- c("Longitude", "Latitude", "Months")

A_Af_df<-as.data.frame(A_Af, xy=TRUE)
A_Af_df<-na.omit(A_Af_df)
colnames(A_Af_df) <- c("Longitude", "Latitude", "Months")

##B##
B_end_df<-as.data.frame(B_end, xy=TRUE)
B_end_df<-na.omit(B_end_df)
colnames(B_end_df) <- c("Longitude", "Latitude", "Months")

B_Af_df<-as.data.frame(B_Af, xy=TRUE)
B_Af_df<-na.omit(B_Af_df)
colnames(B_Af_df) <- c("Longitude", "Latitude", "Months")


##C##
C_end_df<-as.data.frame(C_end, xy=TRUE)
C_end_df<-na.omit(C_end_df)
colnames(C_end_df) <- c("Longitude", "Latitude", "Months")

C_Af_df<-as.data.frame(C_Af, xy=TRUE)
C_Af_df<-na.omit(C_Af_df)
colnames(C_Af_df) <- c("Longitude", "Latitude", "Months")

##D##
D_end_df<-as.data.frame(D_end, xy=TRUE)
D_end_df<-na.omit(D_end_df)
colnames(D_end_df) <- c("Longitude", "Latitude", "Months")

D_Af_df<-as.data.frame(D_Af, xy=TRUE)
D_Af_df<-na.omit(D_Af_df)
colnames(D_Af_df) <- c("Longitude", "Latitude", "Months")

####
#Doing some viz from https://www.geeksforgeeks.org/how-to-make-world-map-with-ggplot2-in-r/

# load library tidyverse 
library(tidyverse) 
library(maps)
library(mapproj)
library(mapdata)
library(ggthemes)

# create data for world coordinates using  
# map_data() function 
#https://datavizpyr.com/how-to-make-world-map-with-ggplot2-in-r/
world_coordinates <- map_data("world") 
world_map = map_data("world")

################################################################################
#Example Plot to write to tif file nicely for publication

tiff('D_E_may50.tif',  units='in', width=8, height=8, res=300, compression = 'lzw')
ggplot() +
 geom_polygon(data = world_map, aes(x = long, y = lat, group = group), color=NA, 
              fill = "darkgrey") +
 coord_fixed(1.1) +
 coord_fixed(xlim = c(45, 95), ylim=c(5,40)) +

geom_tile(data=D_end_df, aes(x = Longitude, y = Latitude, fill = Months)) +
  scale_fill_viridis_c(option  = "F", direction = -1, name="Months") +

geom_polygon(data = world_map, aes(x = long, y = lat, group = group), color="grey", 
             fill = "NA") +
  coord_fixed(1.1) +
  coord_fixed(xlim = c(45, 95), ylim=c(5,40)) +
  theme_bw() +
  labs(x = "Lon", y= "Lat")
dev.off()
###############################################################################
  