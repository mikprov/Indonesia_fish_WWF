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
n.site <- 100
alpha0 <- log(1.1)
sigma.s <- 0.25
alpha <- 200
beta <- 2
Mu <- rnorm(n.site, exp(alpha0), sigma.s)
N.s <- rpois(n.site, Mu)
N.s
Y.s <- N.s
for(i in 1:length(Y.s)){
  if(Y.s[i]>0){
    Y.s[i] <- sum(rgamma(Y.s[i], alpha, beta))
  }
}
hist(Y.s)
Y.data <- data.frame(site=c(1:n.site), biomass=Y.s)

####
####  Quantities of interest, derived
####
E.pres.prob <- 1-exp(-Mu)
mean(1-E.pres.prob)
E.pos.biomass <- ((Mu*alpha)/beta)*(1/(1-exp(-Mu)))
E.biomass <- (Mu*alpha)/beta


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
    proba[j] <- 1 - exp(-mupred[abse[j]])
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
  for(i in 1:nsite){
    mu.new[i] <- exp(alpha + eps[i])
    ngis.new[i] ~ dpois(mu.new[i])
    y.new ~ ifelse()    
if(ngis.new[i]>0){y.new[i] ~ dgamma(a*ngis.new[i], b)}
    if(ngis.new[1]==0){y.new[i] <- 0}
  }

  
}" #end model



####
####  Send model to JAGS and fit
####
abse <- which(Y.data[,"biomass"]==0)
pres <- which(Y.data[,"biomass"]>0)
datalist <- list(Y=Y.data[,"biomass"], nsite=nrow(Y.data),
                 abse=abse, pres=pres, nabs=length(abse), npres=length(pres))
pars <- c("alpha", "a", "b", "y.new")

nAdapt <- 1000
nIter <- 20000
nSamp <- 20000
jm1 <- jags.model(textConnection(model.string), data=datalist, n.chains=1, n.adapt = nAdapt)
update(jm1, n.iter=nIter)
zm <- coda.samples(jm1, variable.names=pars, n.iter=nSamp, n.thin=100)




