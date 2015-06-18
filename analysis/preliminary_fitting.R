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
all_data <- all_data[which(all_data$value>0),]
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
}
transformed parameters{
real mu[obs_n];
for(n in 1:obs_n)
mu[n] <- inv_logit(a[yid[n]] + fint[fid[n]] + b1[fid[n]]*X1[n]);
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

//Likelihood
  Y ~ binomial(1,mu);
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
inits[[1]] <- list(a_mu = 0, a=rep(0,yrs_n), b1_mu=0, b1=rep(0,fam_n),
                   fint = rep(0,fam_n), sig_b1=0.5, sig_a=0.5, sig_m=0.5, sig_f=0.5)

datalist <- list(obs_n=obs_n, yrs_n=yrs_n, yid=yid, covs_n=covs_n,
                 fam_n=fam_n, fid=fid, Y=Y, X1=X[,1])
pars=c("a_mu", "a", "b1_mu",  "b1", "sig_b1", "sig_a", "fint", "sig_f")

stanmodel <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)

mcmc_samples <- stan(fit=stanmodel, data=datalist, init=inits,
                     pars=pars, chains=1, iter=2000, warmup=1000)
print(mcmc_samples)
traceplot(mcmc_samples)
