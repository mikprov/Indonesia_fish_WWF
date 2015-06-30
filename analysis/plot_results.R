## Script to plot results from model fits

library(ggmcmc)
library(ggplot2)
library(plyr)
library(reshape2)


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
variable_names <- colnames(X[,2:ncol(X)])


####
##  Bring in model fit results
####
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
beta_mus$id <- substr(beta_mus$Parameter, 6, length(beta_mus$Parameter))
beta_mus$id <- as.numeric(unlist(strsplit(beta_mus$id, split=']')))
beta_mus <- beta_mus[with(beta_mus, order(id)), ]

png("global_covariates_caterpillar.png", width = 6, height = 4.5,
    units = "in", res=200)
g <- ggplot(beta_mus)+
  geom_hline(aes(yintercept=0), linetype=2)+
  geom_errorbar(aes(x=Parameter, ymin=lo95, ymax=up95), width=0.0001)+
  geom_errorbar(aes(x=Parameter, ymin=lo75, ymax=up75), 
                width=0.0001, size=1.5)+
  geom_point(aes(x=Parameter, y=mean_value), size=4)+
  scale_x_discrete(labels=variable_names)+
#   theme(axis.text.x = element_text(angle = 40, hjust = 1))+
  ylab("Standardized Coefficient Value")+
  xlab("Covariate")+
  ggtitle("Covariate Effects on Fish Biomass \n(Across Families)")+
  coord_flip()
print(g)
dev.off()