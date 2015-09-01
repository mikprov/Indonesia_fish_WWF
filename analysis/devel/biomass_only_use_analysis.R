##  Script to fit a hierarchical model for fish biomass
##
##  1. Fits a biomass model for USE zones only
##  2. Holds all covariates constant except distance to settlement to look
##     at effect on ecah family's biomass.


##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  9-1-205



# Clear workspace 
rm(list=ls(all=TRUE))

####
####  Load libraries -----------------------------------------------------------
####
library(rjags); library(ggplot2)
library(plyr);  library(reshape2)
library(xlsx);  library(ggmcmc)



####
####  Read in data -------------------------------------------------------------
####
##  Read in biomass data
datapath <- "/Users/atredenn/Dropbox/Fish Drivers Paper/DATA/"
file <- "fish_biomass_all_years_fishdrivers_6.17.2015.xlsx"
biomass_raw <- read.xlsx(paste0(datapath,file), sheetIndex = 1)

##  Remove some columns
toremove <- c("Contextual_variables", "sample.event", 
              "num.transects.per.site.fdata", "MPA")
rmids <- which(colnames(biomass_raw) %in% toremove)
biomass_sub <- biomass_raw[,-rmids]

##  Convert data to long format
biomass <- melt(biomass_sub, id.vars = c("zone", "site_id", "year", "depth"))
biomass[which(biomass$zone=="USE"), "zone"] <- "Use"

##  Read in covariates
datapath <- "/Users/atredenn/Dropbox/Fish Drivers Paper/DATA/"
file <- "contextual_variables_fishdrivers_6.17.2015.xlsx"
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
####  Subset data for USE and non-zero biomass ---------------------------------
####
# Subset for Use zones and family-level biomass
model_data <- subset(all_data, zone=="Use" & variable!="average_biomass")
# Only keep non-zero biomass values
model_data <- model_data[which(model_data[,"value"]>0),]



####
####  Set up predictor matrix
####
# Subset out the predictors
x_matrix <- all_data[,grep("X", names(all_data))] # all predictors
toscale <- sapply(x_matrix, is.numeric) # id numeric predictors for scaling
x_matrix[toscale] <- lapply(x_matrix[toscale], scale) # scale numerics
all_data[,grep("X", names(all_data))] <- x_matrix # replace with scaled version

# Make design matrix for predictors
preds <- paste(names(all_data[9:19]), sep="", collapse="+") 
model <- paste("value", preds, sep="~")
X <- model.matrix(as.formula(model) , data = all_data)



####
####  Write JAGS model ---------------------------------------------------------
####
model.string <- "
model{

  for(i in 1:nobs){
    ### Process model (linear regression)
    mu[i] <- a[yid[i]] + fint[fid[i]] + sum(b[fid[i],] * X[i, ])
    z[i] ~ dnorm(mu[i], prec.varsigma)
    
    ### Observation model (likelihood)
    ## Positive biomass data
    y[i] ~ dnorm(mu[i], prec.sigma)
  } #end observation loop
    
  ### Priors
  a_mu ~ dnorm(0,0.001)
  for(y in 1:nyrs){
    a[y] ~ dnorm(a_mu, sig_a)
  } #end years loop for random intercepts
  
  for(k in 1:ncovs){
    b_mu[k] ~ dnorm(0,0.001)
    sig_b[k] ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior
      for(f in 1:nfam){
        b[f,k] ~ dnorm(b_mu[k], sig_b[k])
      } #end family level loop
  } #end covariate loops

  sig_a ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior
  sig_f ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior
  prec.sigma ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior
  prec.varsigma ~ dt(0, 0.1, 1) T(0,) #half-cauchy prior
  for(f in 1:nfam){
    fint[f] ~ dnorm(0, sig_f)
  }

} #end model" 



####
####  Send model to JAGS and fit
####
y <- log(model_data[,"value"])
X <- X[,2:ncol(X)] #take out the added intercept column, we do that ourselves
nobs <- length(y)
ncovs <- ncol(X)
fid <- as.numeric(as.factor(all_data$variable))
nfam <- length(unique(fid))
nyrs <- length(unique(all_data$year))
yid <- as.numeric(as.factor(all_data$year))

datalist <- list(y=y, X=X, nobs=nobs, ncovs=ncovs, nfam=nfam, fid=fid,
                 nyrs=nyrs, yid=yid)
pars=c("a_mu", "a", "b_mu",  "b", "sig_b", "sig_a",
       "sig_m", "fint", "sig_f")

nAdapt <- 5000
nIter <- 10000
nSamp <- 20000
jm1 <- jags.model(textConnection(model.string), data=datalist, n.chains=1, n.adapt = nAdapt)
update(jm1, n.iter=nIter)
zm <- coda.samples(jm1, variable.names=pars, n.iter=nSamp, n.thin=1)
zmgg <- ggs(zm)
ggs_traceplot(zmgg, family="b_mu")

