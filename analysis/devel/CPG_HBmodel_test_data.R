##  Script to test Compound Poisson Gamma model for
##  zero-inflated continuous data (like pres/abs + biomass)
##
##  This is the hierarchical version for fitting betas by species

# Clear workspace 
rm(list=ls(all=TRUE))



####
####  Load libraries -----------------------------------------------------------
####
library(rjags); library(ggplot2)
library(plyr);  library(reshape2)
library(xlsx)



####
####  Read in data -------------------------------------------------------------
####
##  Read in biomass data
datapath <- "/Users/atredenn/Dropbox/Fish Drivers Paper/DATA/"
file <- "fish_biomass_all_years_fishdrivers_9.1.2015.xlsx"
biomass_raw <- read.xlsx(paste0(datapath,file), sheetIndex = 1)

##  Remove some columns
toremove <- c("Contextual_variables", "sample.event", 
              "num.transects.per.site.fdata", "MPA", "hard_coral")
rmids <- which(colnames(biomass_raw) %in% toremove)
biomass_sub <- biomass_raw[,-rmids]

##  Convert data to long format
biomass <- melt(biomass_sub, id.vars = c("zone", "site_id", "year", "depth"))
biomass[which(biomass$zone=="USE"), "zone"] <- "Use"

##  Read in covariates
datapath <- "/Users/atredenn/Dropbox/Fish Drivers Paper/DATA/"
file <- "contextual_variables_fishdrivers_9.1.2015.xlsx"
variables <- read.xlsx(paste0(datapath,file), sheetIndex = 1)

##  Remove some columns
tmpids <- which(names(variables) %in% c("Site.ID", "Lat", "Lon"))
tokeep <- c(tmpids, grep("X", names(variables)))
variables_sub <- variables[,tokeep]
torm <- which(names(variables_sub) %in% c("X4_Exposure_numeric_code", 
                                          "X5_reef_slope_numeric_code", 
                                          "X6_reef_type_num"))
variables_sub <- variables_sub[,-torm]

##  Change some covariates to factors
variables_sub$X10_watershed_pollution_risk <- as.factor(variables_sub$X10_watershed_pollution_risk)
variables_sub$X11_monsoon_direction <- as.factor(variables_sub$X11_monsoon_direction)

##  Merge data frames
all_data <- merge(biomass, variables_sub, by.x = "site_id", by.y = "Site.ID")



####
####  Histrograms of biomass by species ----------------------------------------
####
ggplot(subset(all_data, variable!="average_biomass"), aes(x=value))+
  geom_bar(color="white", fill="purple", alpha=0.5)+
  facet_wrap("variable", scale="free")+
  theme_bw()

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

##  Remove Nemipteridae and average_biomass from data
model_data <- subset(all_data, variable!="Nemipteridae")
model_data <- subset(model_data, variable!="average_biomass")

####
####  Write JAGS model
####
model.string <- "
model{
  ### Latent process
  for(i in 1:nsite){
    log(mupred[i]) <- alpha + beta[fid[i]]*x[i]
  }

  ### Observation model
  ## Positive biomass data
  for(k in 1:npres){
    ## Number of patches (latent)
    isCensored[k] ~ dinterval(npatch[k], 1)
    npatch[k] ~ dpois(mupred[pres[k]])
    ## Observed positive biomass
    Y_a[k] <- a[fid[pres[k]]]*npatch[k]
    Y[pres[k]] ~ dgamma(Y_a[k], b[fid[pres[k]]])
  }
  
  ## Evaluation of probabilities of zeros
  for(j in 1:nabs){
    # Probability of presence at site j
    proba[j] <- 1-exp(-mupred[abse[j]])
    Y[abse[j]] ~ dbern(proba[j])
  }

  ### Priors
  for(f in 1:nfam){
    beta[f] ~ dnorm(beta_mu, sig_b)
    m[f] ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior
    sd[f] ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior
  }
  a <- pow(m,2)/pow(sd,2)
  b <- m/pow(sd,2)
  alpha ~ dnorm(0,0.001)
  beta_mu ~ dnorm(0,0.001)
  sig_b ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior

} #end model" 



####
####  Send model to JAGS and fit
####
Y.data <- model_data
hist(Y.data[,"value"])
abse <- which(Y.data[,"value"]==0)
pres <- which(Y.data[,"value"]>0)
fid <- as.numeric(as.factor(Y.data$variable))-1
nfam <- length(unique(fid))
x <- as.numeric(scale(model_data[,"X8_distance_to_settlement"], center = TRUE, scale = TRUE))
datalist <- list(Y=Y.data[,"value"], x=x, nsite=nrow(Y.data), fid=fid, nfam=nfam,
                 abse=abse, pres=pres, nabs=length(abse), npres=length(pres))
pars <- c("alpha", "a", "b", "beta")

nAdapt <- 100
nIter <- 1000
nSamp <- 1000
jm1 <- jags.model(textConnection(model.string), data=datalist, n.chains=1, n.adapt = nAdapt)
# update(jm1, n.iter=nIter)
zm <- coda.samples(jm1, variable.names=pars, n.iter=nSamp, n.thin=1)
zmd <- as.data.frame(zm[[1]])

par(mfrow=c(4,2))
plot(density(zmd[,"a[14]"]))
plot(zmd[,"a[14]"], type="l")
plot(density(zmd[,"b[14]"]))
plot(zmd[,"b[14]"], type="l")
plot(density(zmd[,"alpha"]))
plot(zmd[,"alpha"], type="l")
plot(density(zmd[,"beta[14]"]))
plot(zmd[,"beta[14]"], type="l")


