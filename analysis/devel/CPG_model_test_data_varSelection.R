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
sub_data <- subset(all_data, variable=="Serranidae")
# sub_data <- sub_data[, c("site_id", "value", "X8_distance_to_settlement")]


####
####  Set up predictor matrix
####
# Subset out the predictors
all.preds <- sub_data[,grep("X", names(sub_data))] # all predictors
x_matrix <- sub_data[,grep("X", names(sub_data))] # all predictors
toscale <- sapply(x_matrix, is.numeric) # id numeric predictors for scaling
x_matrix[toscale] <- lapply(x_matrix[toscale], scale) # scale numerics
sub_data[,grep("X", names(sub_data))] <- x_matrix # replace with scaled version

# Make design matrix for predictors
preds <- paste(names(sub_data[9:19]), sep="", collapse="+") 
model <- paste("value", preds, sep="~")
X <- model.matrix(as.formula(model) , data = sub_data)



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
    log(mupred[i]) <- alpha + inprod(X[i,],beta[])
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
  alpha~dnorm(0,0.0001)
  taub~dgamma(0.1,0.1)
  for(j in 1:ncovs){
    pind[j]~dbeta(2,8)
    ind[j]~dbern(pind[j])
    betaT[j]~dnorm(0,taub)
    beta[j]<-ind[j]*betaT[j]
  }
} #end model" 



####
####  Send model to JAGS and fit
####
Y.data <- sub_data
hist(Y.data[,"value"])
abse <- which(Y.data[,"value"]==0)
pres <- which(Y.data[,"value"]>0)
Xmod <- X[,2:ncol(X)]
datalist <- list(Y=Y.data[,"value"], nsite=nrow(Y.data), X=Xmod, ncovs=ncol(Xmod),
                 abse=abse, pres=pres, nabs=length(abse), npres=length(pres))
pars <- c("alpha", "a", "b", "betaT", "ind")

nAdapt <- 1000
nIter <- 5000
nSamp <- 10000
jm1 <- jags.model(textConnection(model.string), data=datalist, n.chains=1, n.adapt = nAdapt)
update(jm1, n.iter=nIter)
zm <- coda.samples(jm1, variable.names=pars, n.iter=nSamp, n.thin=20)
zmd <- as.data.frame(zm[[1]])
summary(zm)[[1]]

par(mfrow=c(3,2))
plot(density(zmd[,"a"]))
plot(zmd[,"a"], type="l")
plot(density(zmd[,"b"]))
plot(zmd[,"b"], type="l")
plot(density(zmd[,"alpha"]))
plot(zmd[,"alpha"], type="l")

par(mfrow=c(4,4))
beta.cols <- grep("beta", colnames(zmd))
for(i in 1:length(beta.cols)){
  plot(zmd[,beta.cols[i]], type="l")
}



####
####  Plot variable importance -------------------------------------------------
####
out.stats <- summary(zm)[[1]]
var.imp <- out.stats[grep("ind", rownames(out.stats)), "Mean"]
var.imp <- data.frame(variable = colnames(Xmod), importance=var.imp)
var.imp <- var.imp[with(var.imp, order(-importance)), ]
var.imp$rank <- rev(c(1:nrow(var.imp)))
var.imp$keep <- "no"
var.imp[which(var.imp$importance>0.2),"keep"] <- "yes"

library(ggthemes)
ggplot(var.imp, aes(x=as.factor(rank), y=importance, fill=keep))+
  geom_bar(stat="identity", alpha=0.5)+
  scale_fill_manual(values=c("grey", "blue"))+
  scale_x_discrete(labels = rev(var.imp$variable))+
  guides(fill=FALSE)+
  ylab("Mean posterior variable importance")+
  xlab("Predictor variable")+
  coord_flip()+
  theme_bw()


####
####  Traceplots for important predictors --------------------------------------
####
post.preds <- zmd[,grep("beta", colnames(zmd))]
var.imp2 <- out.stats[grep("ind", rownames(out.stats)), "Mean"]
var.imp2 <- data.frame(variable = colnames(Xmod), importance=var.imp2)
var.imp2$keep <- "no"
var.imp2[which(var.imp2$importance>0.2),"keep"] <- "yes"
cols_i_want <- which(var.imp2$keep=="yes")
imp.preds <- post.preds[,cols_i_want]
imp.preds$iteration <- c(1:nrow(imp.preds))
imp.preds <- melt(imp.preds, id.vars = "iteration")

ggplot(imp.preds, aes(x=iteration, y=value))+
  geom_line()+
  facet_wrap("variable", scales="free")

ggplot(imp.preds, aes(x=value))+
  geom_density(color=NA, fill="blue", alpha=0.5)+
  facet_wrap("variable", scales="free")
