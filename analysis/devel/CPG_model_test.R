##  Script to test Compound Poisson Gamma model for
##  zero-inflated continuous data (like pres/abs + biomass)

# Clear workspace 
rm(list=ls(all=TRUE))



####
####  Load libraries
####
library(rjags); library(ggplot2)
library(plyr);  library(reshape2)



####
####  Simulate data
####
set.seed(123)
n.site <- 200
alpha0 <- -0.1
sigma <- 1
eps0 <- rnorm(n.site,0,sigma)
alpha <- 200
beta <- 2
Mu <- exp(alpha0)# + eps0)
N.s <- rpois(n.site, Mu)
Y.s <- N.s
for(i in 1:length(Y.s)){
  if(Y.s[i]>0){
    Y.s[i] <- sum(rgamma(Y.s[i], alpha, beta))
  }
}
hist(Y.s)
length(which(Y.s==0))
Y.data <- data.frame(site=c(1:n.site), biomass=Y.s)



####
####  Quantities of interest, derived
####
E.pres.prob <- 1-exp(-Mu)
mean(1-E.pres.prob)
E.pos.biomass <- ((Mu*alpha)/beta)*(1/(1-exp(-Mu)))
E.biomass <- (Mu*alpha)/beta
mean(E.biomass)
mean(E.pos.biomass)



####
####  Write JAGS model
####
model.string <- "
model{
  ### Latent process
  for(i in 1:nsite){
    log(mupred[i]) <- alpha + eps[i]
  }

  ### Observation model
  ## Positive biomass data
  for(k in 1:npres){
    ## Number of patches (latent)
    npatch[k] ~ dpois(mupred[pres[k]]) T(1,)
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
  a ~ dgamma(2,5)
  b ~ dgamma(2,5)
  alpha ~ dnorm(0, 0.01)
  sigma ~ dunif(0,100)
  for(i in 1:nsite){
    eps[i] ~ dnorm(0, sigma)
  }

  ### Derived quantities

} #end model" 



####
####  Send model to JAGS and fit
####
abse <- which(Y.data[,"biomass"]==0)
pres <- which(Y.data[,"biomass"]>0)
datalist <- list(Y=Y.data[,"biomass"], nsite=nrow(Y.data),
                 abse=abse, pres=pres, nabs=length(abse), npres=length(pres))
pars <- c("alpha", "a", "b", "npatch")

nAdapt <- 100
nIter <- 1000
nSamp <- 2000
jm1 <- jags.model(textConnection(model.string), data=datalist, n.chains=1, n.adapt = nAdapt)
update(jm1, n.iter=nIter)
zm <- coda.samples(jm1, variable.names=pars, n.iter=nSamp, n.thin=1)
head(zm)
zmd <- as.data.frame(zm[[1]])
