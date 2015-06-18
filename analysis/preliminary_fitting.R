## Date: April 17, 2015

## Title: Preliminary fit to look at variograms

##  Authors:  Mikaela Provost
##            Andrew Tredennick

####
####  Load libraries --------------------------------------------
####

library(gstat)
library(sp)
library(xlsx)
library(ggplot2)
library(reshape2)
library(lme4)
library(rstan)
library(parallel)
library(ggmcmc)

####
##  Bring in data ---------------------------------------------
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

##  Merge data frames
all_data <- merge(biomass, variables_sub, by.x = "site_id", by.y = "Site.ID")


####
##  Run linear model on average biomass as test
####
avg_data <- subset(all_data, variable=="average_biomass")
x_matrix <- avg_data[,grep("X", names(avg_data))]
toscale <- sapply(x_matrix, is.numeric)
x_matrix[toscale] <- lapply(x_matrix[toscale], scale)
factors <- sapply(x_matrix, is.factor)
x_matrix[factors] <- lapply(x_matrix[factors], as.numeric)
x_matrix <- as.matrix(x_matrix)

lmmod <- lm(log(avg_data$value)~x_matrix)
summary(lmmod)
mixmod <- lmer(log(avg_data$value)~x_matrix+(1|avg_data$year))
summary(mixmod)


####
##  Set up Stan model
####
##  STAN model
model_string <- "
data{
  int<lower=0> obs_n; // observation length
  int<lower=0> yrs_n; // number of years
  int<lower=0> yid[obs_n]; // year id
  int<lower=0> covs_n; // number of covariates
  int<lower=0> fam_n; // number of families
  int<lower=0> fid[obs_n]; // family id
  vector[obs_n] Y; // observation vector
  vector[obs_n] X1; // design matrix (covariates)
}
parameters{
  real a_mu;
  vector[yrs_n] a;
  real b1_mu;
  vector[fam_n] b1;
  vector[fam_n] fint;
  real<lower=0.00001> sig_a;
  real<lower=0.00001> sig_b1;
  real<lower=0.00001> sig_f;
  real<lower=0.00001> sig_m;
}
transformed parameters{
  real mu[obs_n];
  for(n in 1:obs_n)
    mu[n] <- a[yid[n]] + fint[fid[n]] + b1[fid[n]]*X1[n];
}
model{
  // Priors
  a_mu ~ uniform(-100,100);
  a ~ normal(a_mu, sig_a);
  b1_mu ~ uniform(-100,100);
    for(f in 1:fam_n){
      b1[f] ~ normal(b1_mu, sig_b1);
    }
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  sig_f ~ cauchy(0,5);
  fint ~ normal(0, sig_f);
  sig_m ~ cauchy(0,5);

  //Likelihood
  Y ~ normal(mu, sig_m);
}
"

all_data <- subset(all_data, variable!="average_biomass")
# all_data <- all_data[which(all_data$value>0),]
x_matrix <- all_data[,grep("X", names(all_data))]
toscale <- sapply(x_matrix, is.numeric)
x_matrix[toscale] <- lapply(x_matrix[toscale], scale)
factors <- sapply(x_matrix, is.factor)
x_matrix[factors] <- lapply(x_matrix[factors], as.numeric)
x_matrix <- as.matrix(x_matrix)

obs_n <- nrow(all_data)
yrs_n <- length(unique(all_data$year))
yid <- as.numeric(as.factor(all_data$year))
covs_n <- ncol(x_matrix)
fam_n <- length(unique(all_data$variable))
fid <- as.numeric(as.factor(all_data$variable))-1
Y <- log(all_data$value)
X <- x_matrix

## Set reasonable initial values for three chains
inits <- list()
inits[[1]] <- list(a_mu = 0, a=rep(0,yrs_n), b1_mu=0, b1=rep(0,fam_n),
              fint = rep(0,fam_n), sig_b1=0.5, sig_a=0.5, sig_m=0.5, sig_f=0.5)

datalist <- list(obs_n=obs_n, yrs_n=yrs_n, yid=yid, covs_n=covs_n,
                 fam_n=fam_n, fid=fid, Y=Y, X1=X[,1])
pars=c("a_mu", "a", "b1_mu",  "b1", "sig_b1", "sig_a", "sig_m", "fint", "sig_f")

stanmodel <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)

mcmc_samples <- stan(fit=stanmodel, data=datalist, init=inits,
                     pars=pars, chains=1, iter=2000, warmup=1000)
print(mcmc_samples)
traceplot(mcmc_samples)





# 
# ####
# ##  Some quick plots -------------------
# ####
# ggplot(biomass, aes(x=zone, y=log(value)))+
#   geom_boxplot()+
#   facet_grid(variable~depth)
# 
# 
# ####
# ####  Build dataframe for model function-------------------------
# ####
# 
# # identify columns to not include in linear model
# toremove <- which(colnames(fish_data) %in% c("Include.","Site.ID","MPA","Treatment_original","Site.Name","Lat","Lon",
#                                              "X4_Exposure_num","X4a_Exposure_ExpSemi","X4b_Exposure_SemiShel","X5_reef_slope_num",
#                                              "X5a_reef_slope_fs","X5b_reef_slope_ws", "X6_reef_type_num","X6a_reef_type_pf",
#                                              "X6b_reef_type_fb","X6c_reef_type_ba","Reef.Facing","Comment","X"))
# # removing the columns
# fish_data1 <- fish_data[,-toremove] 
# 
# ####
# ####  Create linear model --------------------------------------
# ####
# 
# # output variable = Biomass
# # input variables = all covariates
# m1 <- lm(log(Biomass) ~ ., data = fish_data1)
# summary(m1)
# 
# ####
# ####  Plot veriogram -------------------------------------------
# ####
# 
# # create dataframe for variogram (colunms: residuals, Lat, Lon)
# residuals <- resid(m1) #pull out residuals
# fish_data1$Lat <- fish_data$Lat #re-add Lat column to fish data
# fish_data1$Lon <- fish_data$Lon #re-add Lon column to fish data
# fish_data_no_na <- fish_data[complete.cases(fish_data1),] #remove NAs from fish data
# 
# # The residuals with locations has to be a data frame for coordinates fxn to work
# resids_df <- data.frame(Lon = fish_data_no_na$Lon, 
#                         Lat = fish_data_no_na$Lat, 
#                         residuals = residuals)
# coordinates(resids_df) <- c("Lon","Lat")
# 
# 
# ####
# ####  Run variogram and plot result -------------------------------
# ####
# varMod <- variogram(residuals~1, data=resids_df)
# plot(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), 
#      col="dodgerblue", type="l")
# points(varMod$dist, varMod$gamma, ylim=c(0,max(varMod$gamma)), 
#        col="dodgerblue", pch=21, bg="white")
# 
