## Script to plot results from model fits

library(ggmcmc)
library(ggplot2)
library(plyr)
library(reshape2)

long <- readRDS("fishbiomass_stanmcmc.RDS")

short <- ddply(long, .(Parameter), summarise,
               mean_value = mean(value),
               sd_value = sd(value),
               lo95 = quantile(value, 0.025),
               up95 = quantile(value, 0.975),
               lo75 = quantile(value, 0.125),
               up75 = quantile(value, 0.875))

betas <- short[grep("b", short$Parameter),]
beta_mus <- betas[grep("b_mu", betas$Parameter),]
ggplot(beta_mus)+
  geom_hline(aes(yintercept=0), linetype=2)+
  geom_errorbar(aes(x=Parameter, ymin=lo95, ymax=up95), width=0.0001)+
  geom_errorbar(aes(x=Parameter, ymin=lo75, ymax=up75), 
                width=0.0001, size=1.5)+
  geom_point(aes(x=Parameter, y=mean_value), size=4)+
  coord_flip()