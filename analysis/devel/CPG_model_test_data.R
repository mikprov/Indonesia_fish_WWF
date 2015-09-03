##  Script to test Compound Poisson Gamma model for
##  zero-inflated continuous data (like pres/abs + biomass)

# Clear workspace 
rm(list=ls(all=TRUE))

####
####  Load libraries
####
library(rjags); library(ggplot2)
library(plyr);  library(reshape2)
library(xlsx)



####
####  Read in data
####
biomass_raw <- read.xlsx("/Users/atredenn/Dropbox/Fish Drivers Paper/DATA/fish_biomass_all_years_fishdrivers_6.17.2015.xlsx",
                         sheetIndex = 1)
toremove <- c("Contextual_variables", "sample.event", "num.transects.per.site.fdata",
              "MPA")
rmids <- which(colnames(biomass_raw) %in% toremove)
biomass_sub <- biomass_raw[,-rmids]
biomass <- melt(biomass_sub, id.vars = c("zone", "site_id", "year", "depth"))
biomass[which(biomass$zone=="USE"), "zone"] <- "Use"

variables <- read.xlsx("/Users/atredenn/Dropbox/Fish Drivers Paper/DATA/contextual_variables_fishdrivers_6.17.2015.xlsx",
                       sheetIndex = 1)
tmpids <- which(names(variables) %in% c("Site.ID", "Lat", "Lon"))
tokeep <- c(tmpids, grep("X", names(variables)))
variables_sub <- variables[,tokeep]
torm <- which(names(variables_sub) %in% c("X4_Exposure_numeric_code", 
                                          "X5_reef_slope_numeric_code", 
                                          "X6_reef_type_num"))
variables_sub <- variables_sub[,-torm]
variables_sub$X10_watershed_pollution_risk <- as.factor(variables_sub$X10_watershed_pollution_risk)
variables_sub$X11_monsoon_direction <- as.factor(variables_sub$X11_monsoon_direction)

##  Merge data frames
all_data <- merge(biomass, variables_sub, by.x = "site_id", by.y = "Site.ID")
sub_data <- subset(all_data, variable=="Carangidae")
sub_data <- sub_data[, c("site_id", "value", "X8_distance_to_settlement")]


####
####  Histrograms of biomass by species
####
ggplot(all_data, aes(x=value))+
  geom_bar()+
  facet_wrap("variable", scale="free")

# Plot proportion of sample with zeros for each species
count=1
out <- numeric(length(unique(all_data$variable)))
for(i in unique(all_data$variable)){
  tmp <- subset(all_data, variable==i)
  tmpzs <- length(which(tmp[,"value"]==0))/length(tmp[,"value"])
  out[count] <- tmpzs
  count <- count+1
  print(i)
}
zeros.data <- data.frame(family=unique(all_data$variable),
                         propzs=out)
ggplot(data=zeros.data, aes(x=family, y=propzs))+
  geom_bar(stat="identity")+
  coord_flip()+
  ylab("Proportion of observations that are 0")+
  xlab("Fish family")

subset(all_data, variable=="Nemipteridae")[,"value"]

####
####  Write JAGS model
####
model.string <- "
model{
  ### Latent process
  for(i in 1:nsite){
    log(mupred[i]) <- alpha 
  }

  ### Observation model
  ## Positive biomass data
  for(k in 1:npres){
    ## Number of patches (latent)
    isCensored[k] ~ dinterval(npatch[k], 1)
    npatch[k] ~ dpois(mupred[pres[k]])
    ## Observed positive biomass
    Y_a[k] <- a*npatch[k]
    Y[pres[k]] ~ dgamma(Y_a[k], b)
  }
  
  ## Evaluation of probabilities of zeros
  for(j in 1:nabs){
    # Probability of presence at site j
    proba[j] <- 1-exp(-mupred[abse[j]])
    Y[abse[j]] ~ dbern(proba[j])
  }

  ### Priors
  a <- pow(m,2)/pow(sd,2)
  b <- m/pow(sd,2)
  m ~ dunif(0,100)
  sd ~ dunif(0,100)
  alpha ~ dnorm(0, 0.01)

} #end model" 



####
####  Send model to JAGS and fit
####
Y.data <- sub_data
hist(Y.data[,"value"])
abse <- which(Y.data[,"value"]==0)
pres <- which(Y.data[,"value"]>0)
datalist <- list(Y=Y.data[,"value"], nsite=nrow(Y.data),
                 abse=abse, pres=pres, nabs=length(abse), npres=length(pres))
pars <- c("alpha", "a", "b")

nAdapt <- 100
nIter <- 1000
nSamp <- 2000
jm1 <- jags.model(textConnection(model.string), data=datalist, n.chains=1, n.adapt = nAdapt)
update(jm1, n.iter=nIter)
zm <- coda.samples(jm1, variable.names=pars, n.iter=nSamp, n.thin=1)
zmd <- as.data.frame(zm[[1]])

par(mfrow=c(3,2))
plot(density(zmd[,"a"]))
plot(zmd[,"a"], type="l")
plot(density(zmd[,"b"]))
plot(zmd[,"b"], type="l")
plot(density(zmd[,"alpha"]))
plot(zmd[,"alpha"], type="l")

