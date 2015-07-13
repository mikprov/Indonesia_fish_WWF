##  Script to plot residual variogram from biomass model to look
##    for spatial autocorrelation.

library(ggmcmc)
library(ggplot2)
library(plyr)
library(reshape2)
library(xlsx)


####
####  Bring in data
####
biomass_raw <- read.xlsx("/Users/atredenn/Dropbox/Fish Drivers Paper/DATA/fish_biomass_all_years_fishdrivers_6.17.2015.xlsx",
                         sheetIndex = 1)
toremove <- c("Contextual_variables", "sample.event", "num.transects.per.site.fdata",
              "MPA")
rmids <- which(colnames(biomass_raw) %in% toremove)
biomass_sub <- biomass_raw[,-rmids]
biomass <- melt(biomass_sub, id.vars = c("zone", "site_id", "year", "depth"))
biomass[which(biomass$zone=="USE"), "zone"] <- "Use"

rmbiom <- grep("average_biomass", biomass$variable)
fam_d <- biomass[-rmbiom, ]
families <- unique(fam_d$variable)

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


####
####  Make model predictions from coefficient point estimates
####
intercept <- short[grep("a_mu", short$Parameter),]
betas <- short[grep("b", short$Parameter),]
murms <- grep("b_mu", betas$Parameter)
beta_fams <- betas[-murms,]
sigrms <- grep("sig", beta_fams$Parameter)
beta_fams <- beta_fams[-sigrms,]
beta_mus <- betas[murms,]
beta_mus$id <- substr(beta_mus$Parameter, 6, length(beta_mus$Parameter))
beta_mus$id <- as.numeric(unlist(strsplit(beta_mus$id, split=']')))
beta_mus <- beta_mus[with(beta_mus, order(id)), ]

## For beta_fams -- family ID is the row, coefficient ID is column
beta_fams$coef_id <- NA
beta_fams$fam_id <- NA
for(i in 1:nrow(beta_fams)){
  tmp1 <- na.omit(as.numeric(unlist(strsplit(unlist(as.character(beta_fams$Parameter[i])), "[^0-9]+"))))
  beta_fams$fam_id[i] <- tmp1[1]
  beta_fams$coef_id[i] <- tmp1[2]
}
beta_fams <- beta_fams[with(beta_fams, order(fam_id, coef_id)), ]
families <- as.data.frame(families)
families$fam_id <- c(1:nrow(families))
variables <- as.data.frame(variable_names)
variables$coef_id <- c(1:nrow(variables))

beta_fams_coded <- merge(beta_fams, families, by="fam_id")
beta_fams_coded <- merge(beta_fams_coded, variables, by="coef_id")

Xdf <- as.data.frame(X)
Xdf$family <- all_data$variable

family_list <- unique(all_data$variable)
pred_list <- list()
for(do_fam in family_list){
  Xtmp <- subset(Xdf, family==do_fam)
  beta <- subset(beta_fams_coded, families==do_fam)[,"mean_value"]
  betas <- c(intercept[,"mean_value"], beta)
  xtmp <- as.matrix(Xtmp[,1:length(betas)])
  tmp_preds <-  as.data.frame(xtmp%*%betas)
  colnames(tmp_preds) <- "log_prediction"
  tmp_preds$id <- row.names(tmp_preds)
  tmp_preds$family <- do_fam
  pred_list[[do_fam]] <- tmp_preds
}

pred_df <- do.call(rbind, lapply(pred_list, data.frame, stringsAsFactors=FALSE))
obs_df <- all_data[, c("value", "Lat", "Lon")]
obs_df$id <- row.names(obs_df)
plot_df <- merge(obs_df, pred_df, by="id")
plot_df$prediction <- exp(plot_df$log_prediction)
plot(prediction~value, plot_df, xlim=c(0,200), ylim=c(0,200))
abline(0,1)
