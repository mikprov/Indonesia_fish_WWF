##  Indonesia fish data initial spatial dependence analysis

##  Date:     4-16-2015
##  Authors:  Andrew Tredennick
##            Mikaela Provost

####
####  Load libraries --------------------------------------------
####
library(ggmap)
library(cluster)
library(ggplot2)


####
####  Bring in data ---------------------------------------------
####
fish_data <- read.csv("../data/Fish_and_contextual_var_for_SAMSI.csv")


####
####  Look at spatial data --------------------------------------
####
fish_pts <- fish_data[,c('Lat', 'Lon')]
plot(fish_pts$Lon, fish_pts$Lat)
removes <- which(is.na(fish_pts))
ks <- kmeans(fish_pts[-removes,], centers = 20)
points(ks$centers[,"Lon"], ks$centers[,"Lat"], col="blue", pch=19)
str(ks$cluster)


####
####  Make a ggmap ---------------------------------------------
####
# Now get the base map using 'get_map', the location is also a Google search
mapALL <- get_map(location="selat dampier", source ="stamen", maptype="toner",
                  zoom=7)
# Put it together using 'ggmap' and the standard ggplot call to geom_point
ggmap(mapALL, extent="device")+
  geom_point(data=fish_data, aes(x=Lon, y=Lat), 
             size=4, alpha=0.6, color="dodgerblue")
