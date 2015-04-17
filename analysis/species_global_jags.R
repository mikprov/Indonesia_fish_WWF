##  Generic JAGS code for multi-level Bayesian model

##  This will be used to fit a predictive model to 
##    the WWF fish data. This code will be called by
##    a separate script that actually runs the JAGS model
##    based on defined data. The model currently does not
##    include a spatial effect.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4-16-2015

model{
  #Process model
  for(obs in 1:n_obs){
    mu[obs] <- intercept[species[obs]] + inprod(inprod(beta[1:n_covs],X[i,1:n_covs]))
    Y[obs] ~ dlnorm(mu[obs], tau[species[obs]])
  }
  
  #Priors
  for(spp in 1:n_species){
    intercept[spp] ~ dnorm(intercept_mu, intercept_tau)
  }
  for(covariate in 1:n_covariates) {
    beta_mu[covariate] ~ dnorm(0,1e-6)
    beta_tau[covariate] ~ dgamma(0.001, 0.001)
    for(spp in 1:n_species){
      beta[covariate, spp] ~ dnorm(beta_mu[covariate], beta_tau[covariate])
    }
  }
  intercept_tau ~ dgamma(0.001, 0.001)
}
