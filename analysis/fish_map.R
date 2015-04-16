##  Indonesia fish data

library(ggmap)
fish_data <- read.csv("Fish_and_contextual_var_for_SAMSI.csv")
fish_pts <- fish_data[,c('Lat', 'Lon')]

# Now get the base map using 'get_map', the location is also a Google search
mapALL <- get_map(location="indonesia", source ="stamen", maptype="toner",
                  zoom=6)
# Put it together using 'ggmap' and the standard ggplot call to geom_point
ggmap(mapALL, extent="device")+
  geom_point(data=fish_data, aes(x=Lon, y=Lat, color=Regional.Location), 
             size=4, alpha=0.6)
#   geom_point(data=allD, aes(x=Long, y=Lat),  shape=1, size=4)+
#   scale_color_discrete(name="Region")
#   geom_text(data=allD, aes(x=Long, y=Lat, label=Exclosure.Name),hjust=0, vjust=0)