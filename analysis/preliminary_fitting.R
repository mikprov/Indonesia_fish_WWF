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
torm <- which(names(variables_sub) %in% c("X4_Exposure_numeric_code", 
                                         "X5_reef_slope_numeric_code", 
                                         "X6_reef_type_num"))
variables_sub <- variables_sub[,-torm]
variables_sub$X10_watershed_pollution_risk <- as.factor(variables_sub$X10_watershed_pollution_risk)
variables_sub$X11_monsoon_direction <- as.factor(variables_sub$X11_monsoon_direction)

##  Merge data frames
all_data <- merge(biomass, variables_sub, by.x = "site_id", by.y = "Site.ID")



###########################
##### BIOMASS MODEL #######
###########################
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
  //vector[obs_n] X1; // design matrix (covariates)
  matrix[obs_n, covs_n] X;
}
parameters{
  real a_mu;
  vector[yrs_n] a;
  vector[covs_n] b_mu;
  matrix[fam_n, covs_n] b;
  vector[fam_n] fint;
  real<lower=0.00001> sig_a;
  vector<lower=0.00001>[covs_n] sig_b;
  real<lower=0.00001> sig_f;
  real<lower=0.00001> sig_m;
}
transformed parameters{
  real mu[obs_n];
  real xbeta[obs_n];
  real xtmp[covs_n];
  for(n in 1:obs_n){
    for(k in 1:covs_n){
      xtmp[k] <- b[fid[n],k] * X[n, k];
    }
    xbeta[n] <- sum(xtmp);
    mu[n] <- a[yid[n]] + fint[fid[n]] + xbeta[n];
  }
}
model{
  // Priors
  a_mu ~ normal(0,1000);
  a ~ normal(a_mu, sig_a);
  for(k in 1:covs_n){
    b_mu[k] ~ normal(0,10);
    sig_b[k] ~ cauchy(0,5);
    for(f in 1:fam_n)
      b[f,k] ~ normal(b_mu[k], sig_b[k]);
  }
  sig_a ~ cauchy(0,5);
  sig_f ~ cauchy(0,5);
  fint ~ normal(0, sig_f);
  sig_m ~ cauchy(0,5);

  //Likelihood
  Y ~ normal(mu, sig_m);
}
"

all_data <- subset(all_data, variable!="average_biomass")
all_data <- all_data[which(all_data$value>0),]
x_matrix <- all_data[,grep("X", names(all_data))]
toscale <- sapply(x_matrix, is.numeric)
x_matrix[toscale] <- lapply(x_matrix[toscale], scale)
all_data[,grep("X", names(all_data))] <- x_matrix

## Set up design matrix for predictors
preds <- paste(names(all_data[9:19]), sep="", collapse="+") 
model <- paste("value", preds, sep="~")
X <- model.matrix(as.formula(model) , data = all_data)

## Declare variables for STAN
obs_n <- nrow(all_data)
yrs_n <- length(unique(all_data$year))
yid <- as.numeric(as.factor(all_data$year))
fam_n <- length(unique(all_data$variable))
fid <- as.numeric(as.factor(all_data$variable))-1
Y <- log(all_data$value)
X <- X[,2:ncol(X)] #take out the added intercept column, we do that ourselves
covs_n <- ncol(X)


## Set reasonable initial values for three chains
inits <- list()
inits[[1]] <- list(a_mu = 0, a=rep(0,yrs_n), b_mu=rep(0, covs_n), 
                   b=matrix(0,nrow=fam_n,ncol=covs_n),
                   fint = rep(0,fam_n), sig_b=rep(0.5,covs_n), 
                   sig_a=0.5, sig_m=0.5, sig_f=0.5)
inits[[2]] <- list(a_mu = 1, a=rep(1,yrs_n), b_mu=rep(1, covs_n), 
                   b=matrix(1,nrow=fam_n,ncol=covs_n),
                   fint = rep(1,fam_n), sig_b=rep(0.05,covs_n), 
                   sig_a=0.05, sig_m=0.05, sig_f=0.05)
inits[[3]] <- list(a_mu = 0.5, a=rep(0.5,yrs_n), b_mu=rep(0.5, covs_n), 
                   b=matrix(0.5,nrow=fam_n,ncol=covs_n),
                   fint = rep(0.5,fam_n), sig_b=rep(0.25,covs_n), 
                   sig_a=0.25, sig_m=0.25, sig_f=0.25)

datalist <- list(obs_n=obs_n, yrs_n=yrs_n, yid=yid, covs_n=covs_n,
                 fam_n=fam_n, fid=fid, Y=Y, X=X)
pars=c("a_mu", "a", "b_mu",  "b", "sig_b", "sig_a",
       "sig_m", "fint", "sig_f")

stanmodel <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)

# mcmc_samples <- stan(fit=stanmodel, data=datalist, init=inits,
#                      pars=pars, chains=1, iter=1000, warmup=500)
# print(mcmc_samples)
# plot(mcmc_samples)
# traceplot(mcmc_samples)

##  Run MCMC in parallel
rng_seed <- 123
sflist <-
  mclapply(1:3, mc.cores=3,
           function(i) stan(fit=stanmodel, data=datalist, pars=pars,
                            seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                            iter=2000, warmup=1000, init=list(inits[[i]])))
fit <- sflist2stanfit(sflist)
print(fit)
long <- ggs(fit)
stansumm <- as.data.frame(summary(fit)["summary"])
rhats <- stansumm["summary.Rhat"]

saveRDS(long, "fishbiomass_stanmcmc.RDS")
write.csv(rhats, "fishbiomass_stanmcmc_rhats.csv")



####################################
##### PRESENCE-ABSENCE MODEL #######
####################################
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
  int<lower=0,upper=1> Y[obs_n]; // observation vector
  matrix[obs_n, covs_n] X;
}
parameters{
  real a_mu;
  vector[yrs_n] a;
  vector[covs_n] b_mu;
  matrix[fam_n, covs_n] b;
  vector[fam_n] fint;
  real<lower=0.00001> sig_a;
  vector<lower=0.00001>[covs_n] sig_b;
  real<lower=0.00001> sig_f;
  real<lower=0.00001> sig_m;
  vector[obs_n] mu2;
}
transformed parameters{
  real mu[obs_n];
  real xbeta[obs_n];
  real xtmp[covs_n];
  for(n in 1:obs_n){
    for(k in 1:covs_n){
      xtmp[k] <- b[fid[n],k] * X[n, k];
    }
    xbeta[n] <- sum(xtmp);
    mu[n] <- a[yid[n]] + fint[fid[n]] + xbeta[n];
  }
}
model{
// Priors
  a_mu ~ normal(0,1000);
  a ~ normal(a_mu, sig_a);
  for(k in 1:covs_n){
    b_mu[k] ~ normal(0,10);
    sig_b[k] ~ cauchy(0,5);
    for(f in 1:fam_n){
      b[f,k] ~ normal(b_mu[k], sig_b[k]);
    }
  }
  sig_a ~ cauchy(0,5);
  sig_f ~ cauchy(0,5);
  fint ~ normal(0, sig_f);
  sig_m ~ cauchy(0,5);

  //Model error (dispersion)
  mu2 ~ normal(mu, sig_m);

  //Likelihood
  for(i in 1:obs_n)
    Y[i] ~ binomial(1, inv_logit(mu2[i]));
}
"

all_data <- subset(all_data, variable!="average_biomass")
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
Y <- all_data$value
Y[which(Y>0)] <- 1
Y[which(Y==0)] <- 0
X <- x_matrix

## Set reasonable initial values for three chains
inits <- list()
inits[[1]] <- list(a_mu = 0, a=rep(0,yrs_n), b_mu=rep(0, covs_n), 
                   b=matrix(0,nrow=fam_n,ncol=covs_n),
                   fint = rep(0,fam_n), sig_b=rep(0.5,covs_n), 
                   sig_a=0.5, sig_m=0.5, sig_f=0.5, mu2=rep(0.1, obs_n))
inits[[2]] <- list(a_mu = 1, a=rep(1,yrs_n), b_mu=rep(1, covs_n), 
                   b=matrix(1,nrow=fam_n,ncol=covs_n),
                   fint = rep(1,fam_n), sig_b=rep(0.05,covs_n), 
                   sig_a=0.05, sig_m=0.05, sig_f=0.05, mu2=rep(0.2, obs_n))
inits[[3]] <- list(a_mu = 0.5, a=rep(0.5,yrs_n), b_mu=rep(0.5, covs_n), 
                   b=matrix(0.5,nrow=fam_n,ncol=covs_n),
                   fint = rep(0.5,fam_n), sig_b=rep(0.25,covs_n), 
                   sig_a=0.25, sig_m=0.25, sig_f=0.25, mu2=rep(0.01, obs_n))

datalist <- list(obs_n=obs_n, yrs_n=yrs_n, yid=yid, covs_n=covs_n,
                 fam_n=fam_n, fid=fid, Y=Y, X=X)
pars=c("a_mu", "a", "b_mu",  "b", "sig_b", "sig_a",
       "sig_m", "fint", "sig_f")

stanmodel <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)

mcmc_samples <- stan(fit=stanmodel, data=datalist, init=list(inits[[1]]),
                     pars=pars, chains=1, iter=200, warmup=100)
print(mcmc_samples)
traceplot(mcmc_samples)
