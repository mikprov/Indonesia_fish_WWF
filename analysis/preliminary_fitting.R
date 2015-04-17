## Date: April 17, 2015

## Title: Preliminary fit to look at variograms

##  Authors:  Mikaela Provost
##            Andrew Tredennick

####
####  Load libraries --------------------------------------------
####

library(gstat)
library(sp)

####
####  Bring in data ---------------------------------------------
####
fish_data <- read.csv("../data/Fish_and_contextual_var_for_SAMSI.csv",
                      na.strings=c("NA","NaN","#N/A"))

####
####  Build dataframe for model function-------------------------
####

# identify columns to not include in linear model
toremove <- which(colnames(fish_data) %in% c("Include.","Site.ID","MPA","Treatment_original","Site.Name","Lat","Lon",
                                             "X4_Exposure_num","X4a_Exposure_ExpSemi","X4b_Exposure_SemiShel","X5_reef_slope_num",
                                             "X5a_reef_slope_fs","X5b_reef_slope_ws", "X6_reef_type_num","X6a_reef_type_pf",
                                             "X6b_reef_type_fb","X6c_reef_type_ba","Reef.Facing","Comment","X"))
# removing the columns
fish_data1 <- fish_data[,-toremove] 

####
####  Create linear model --------------------------------------
####

# output variable = Biomass
# input variables = all covariates
m1 <- lm(log(Biomass) ~ ., data = fish_data1)
summary(m1)

####
####  Plot veriogram -------------------------------------------
####

# create dataframe for variogram (colunms: residuals, Lat, Lon)
residuals <- resid(m1) #pull out residuals
fish_data1$Lat <- fish_data$Lat #re-add Lat column to fish data
fish_data1$Lon <- fish_data$Lon #re-add Lon column to fish data
fish_data_no_na <- fish_data[complete.cases(fish_data1),] #remove NAs from fish data

# The residuals with locations has to be a data frame for coordinates fxn to work
resids_df <- data.frame(Lon = fish_data_no_na$Lon, 
                        Lat = fish_data_no_na$Lat, 
                        residuals = residuals)
coordinates(resids_df) <- c("Lon","Lat")


####
####  Run variogram and plot result -------------------------------
####
varMod <- variogram(residuals~1, data=resids_df)
plot(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), 
     col="dodgerblue", type="l")
points(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), 
       col="dodgerblue", pch=21, bg="white")

