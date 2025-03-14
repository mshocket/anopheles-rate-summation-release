## Rate Summation Project
## Script for performing rate summation and suitability calculations
## Written by Marta Shocket in 2022-24

##	Table of contents:
##
##	1. Set-up workspace
##  2. Generate hourly temperatures with the Parton-Logan model
##  3. Perform rate summation calculations on traits
##  4. Define S(T)/R0 function
##  5. Calculate versions of Suitability S(T)
##  6. Process rate summation calculation output for plotting
##  7. Join RS output with corresponding empirical TPCs and plot each panel
##  8. Calculate percent change compared to constant temperature model
##  9. Manuscript Figure 3
##  10. Process Suitability output
##  11. Manuscript Figure 4
##  12. Processing suitability output for mapping



##########
###### 1. Set-up workspace
##########

##### Load Packages
library(progress)
library(gridExtra)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot())

##### Load functions
source("R-scripts/working-versions-code/00_RSProjectFunctions.R")

##### Load fitted TPC models
load("saved-posteriors/constant_lifespan.Rdata")
load("saved-posteriors/constant_gamma.Rdata")
load("saved-posteriors/constant_eggs.Rdata")
load("saved-posteriors/constant_biterate.Rdata")

load("saved-posteriors/dtr9_lifespan.Rdata")
load("saved-posteriors/dtr9_gamma.Rdata")
load("saved-posteriors/dtr9_eggs.Rdata")
load("saved-posteriors/dtr9_biterate.Rdata")

load("saved-posteriors/dtr12_lifespan.Rdata")
load("saved-posteriors/dtr12_gamma.Rdata")
load("saved-posteriors/dtr12_eggs.Rdata")
load("saved-posteriors/dtr12_biterate.Rdata")

load("saved-posteriors/pea.Rdata")
load("saved-posteriors/MDR.Rdata")
load("saved-posteriors/EIP50.Rdata")
load("saved-posteriors/bc.Rdata")

Temp.xs <- seq(0, 45, 0.1) # temperature gradient


##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames here
# for the sensitivity analysis the predictions need to be matrices *NOT* data frames
predictions_bite_rate_constant <- as.data.frame(model_bite_rate_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_bite_rate_dtr9 <- as.data.frame(model_bite_rate_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_bite_rate_dtr12 <- as.data.frame(model_bite_rate_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_lifespan_constant <- as.data.frame(model_lifespan_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_lifespan_dtr9 <- as.data.frame(model_lifespan_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_lifespan_dtr12 <- as.data.frame(model_lifespan_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_eggs_constant <- as.data.frame(model_eggs_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_eggs_dtr9 <- as.data.frame(model_eggs_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_eggs_dtr12 <- as.data.frame(model_eggs_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_gamma_constant <- as.data.frame(model_gamma_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_gamma_dtr9 <- as.data.frame(model_gamma_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_gamma_dtr12 <- as.data.frame(model_gamma_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_bc <- as.data.frame(model_bc$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_eip50 <- as.data.frame(model_eip50$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_pea <- as.data.frame(model_pea$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_mdr <- as.data.frame(model_mdr$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)


##########
###### 2. Generate hourly temperatures with the Parton-Logan model
##########

# Temperature gradient that matches TPC predictions
Temp.xs <- seq(0, 45, 0.1)

# Generate hourly temperature sequences across the temperature gradient
LPtemps_dtr9 <- LoganPartonCalc(dtr = 9, Temp.xs)
LPtemps_dtr12 <- LoganPartonCalc(dtr = 12, Temp.xs)

# Set negative temperature values to 0, since TPCs stop at 0 on the low end
LPtemps_dtr9[LPtemps_dtr9 < 0 ] <- 0
LPtemps_dtr12[LPtemps_dtr12 < 0 ] <- 0

# Set temperature values > 45 to 45, since TPCs stop at 45 on the hight end
LPtemps_dtr9[LPtemps_dtr9 > 45 ] <- 45
LPtemps_dtr12[LPtemps_dtr12 > 45 ] <- 45


##########
###### 3. Perform rate summation calculations on traits
##########

# These take about 1.5-2 minutes each to run

predictions_bite_rate_rs9 <- RSCalcTempGrad(predictions_bite_rate_constant, LPtemps_dtr9, Temp.xs)
predictions_bite_rate_rs12 <- RSCalcTempGrad(predictions_bite_rate_constant, LPtemps_dtr12, Temp.xs)

predictions_lifespan_rs9 <- RSCalcTempGrad(predictions_lifespan_constant, LPtemps_dtr9, Temp.xs)
predictions_lifespan_rs12 <- RSCalcTempGrad(predictions_lifespan_constant, LPtemps_dtr12, Temp.xs)

predictions_eggs_rs9 <- RSCalcTempGrad(predictions_eggs_constant, LPtemps_dtr9, Temp.xs)
predictions_eggs_rs12 <- RSCalcTempGrad(predictions_eggs_constant, LPtemps_dtr12, Temp.xs)

predictions_gamma_rs9 <- RSCalcTempGrad(predictions_gamma_constant, LPtemps_dtr9, Temp.xs)
predictions_gamma_rs12 <- RSCalcTempGrad(predictions_gamma_constant, LPtemps_dtr12, Temp.xs)

predictions_bc_rs9 <- RSCalcTempGrad(predictions_bc, LPtemps_dtr9, Temp.xs)
predictions_bc_rs12 <- RSCalcTempGrad(predictions_bc, LPtemps_dtr12, Temp.xs)

predictions_eip50_rs9 <- RSCalcTempGrad(predictions_eip50, LPtemps_dtr9, Temp.xs)
predictions_eip50_rs12 <- RSCalcTempGrad(predictions_eip50, LPtemps_dtr12, Temp.xs)

predictions_pea_rs9 <- RSCalcTempGrad(predictions_pea, LPtemps_dtr9, Temp.xs)
predictions_pea_rs12 <- RSCalcTempGrad(predictions_pea, LPtemps_dtr12, Temp.xs)

predictions_mdr_rs9 <- RSCalcTempGrad(predictions_mdr, LPtemps_dtr9, Temp.xs)
predictions_mdr_rs12 <- RSCalcTempGrad(predictions_mdr, LPtemps_dtr12, Temp.xs)


# Save trait RS output
save(predictions_bite_rate_rs9, file = "data-processed/predictions_bite_rate_rs9.Rdata")
save(predictions_bite_rate_rs12, file = "data-processed/predictions_bite_rate_rs12.Rdata")

save(predictions_lifespan_rs9, file = "data-processed/predictions_lifespan_rs9.Rdata")
save(predictions_lifespan_rs12, file = "data-processed/predictions_lifespan_rs12.Rdata")

save(predictions_eggs_rs9, file = "data-processed/predictions_eggs_rs9.Rdata")
save(predictions_eggs_rs12, file = "data-processed/predictions_eggs_rs12.Rdata")

save(predictions_gamma_rs9, file = "data-processed/predictions_gamma_rs9.Rdata")
save(predictions_gamma_rs12, file = "data-processed/predictions_gamma_rs12.Rdata")

save(predictions_bc_rs9, file = "data-processed/predictions_bc_rs9.Rdata")
save(predictions_bc_rs12, file = "data-processed/predictions_bc_rs12.Rdata")

save(predictions_eip50_rs9, file = "data-processed/predictions_eip50_rs9.Rdata")
save(predictions_eip50_rs12, file = "data-processed/predictions_eip50_rs12.Rdata")

save(predictions_pea_rs9, file = "data-processed/predictions_pea_rs9.Rdata")
save(predictions_pea_rs12, file = "data-processed/predictions_pea_rs12.Rdata")

save(predictions_mdr_rs9, file = "data-processed/predictions_mdr_rs9.Rdata")
save(predictions_mdr_rs12, file = "data-processed/predictions_mdr_rs12.Rdata")


##### Load saved trait RS output
load("data-processed/predictions_bite_rate_rs9.Rdata")
load("data-processed/predictions_bite_rate_rs12.Rdata")

load("data-processed/predictions_lifespan_rs9.Rdata")
load("data-processed/predictions_lifespan_rs12.Rdata")

load("data-processed/predictions_eggs_rs9.Rdata")
load("data-processed/predictions_eggs_rs12.Rdata")

load("data-processed/predictions_gamma_rs9.Rdata")
load("data-processed/predictions_gamma_rs12.Rdata")

load("data-processed/predictions_bc_rs9.Rdata")
load("data-processed/predictions_bc_rs12.Rdata")

load("data-processed/predictions_eip50_rs9.Rdata")
load("data-processed/predictions_eip50_rs12.Rdata")

load("data-processed/predictions_pea_rs9.Rdata")
load("data-processed/predictions_pea_rs12.Rdata")

load("data-processed/predictions_mdr_rs9.Rdata")
load("data-processed/predictions_mdr_rs12.Rdata")


##########
###### 4. Define S(T)/R0 function
##########

# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48

# R0 formulation where: 1) mosquito density (M) depends on lifetime fecundity (B); 2) gamma (y) is substituted for exp^(-1/(lf*PDR))expression
R0eq = function(a, lf, B, y, bc, pEA, MDR) {
	M = B * pEA * MDR * (lf+ec.lf)
	R0 = (a^2 * bc * y * M * (lf+ec.lf) )^0.5
	return(R0)
}


##########
###### 5. Calculate versions of Suitability S(T)
##########

##### Summary of different versions of suitability calculations
#
# Five different versions of S(T) for specific pairwise comparisons:
# 	1) suit_const:			constant temperature TPCs
# 2) suit_empfluc:		empirically fit TPCs for a, lf, and B in fluctuating environments (other traits varying based on constant temperature TPCs)
# 3) suit_rs_traitsall:	rate summation on all trait TPCs
# 4) suit_rs_3traits:		rate summation on a, lf, and B (other traits varying based on constant temperature TPCs)
# 5) suit_rs_suit: 		rate summation on S(T) TPC
# 
# Pairwise comparisons and what they mean biologically or statistically:
# 1 vs 2 - Does fluctuating temperature (empirically) affect S(T)?
# 1 vs 3 - Does incorporating rate summation (on traits) affect predicted S(T)?
# 2 vs 4 - Does rate summation accurately predict the effect of fluctuating temperature on S(T)?
# 3 vs 5 - Does performing rate summation on traits vs. on S(T) TPC affect predicted S(T)?

# 1. All constant temperature TPCs
suit_const <- R0eq(predictions_bite_rate_constant, predictions_lifespan_constant, predictions_eggs_constant, predictions_gamma_constant, 
				   predictions_bc, predictions_pea, predictions_mdr)

# 2A. Empirically fit TPCs for a, lf, and B in dtr9, other traits constant temperature TPCs
suit_empfluc_dtr09 <- R0eq(predictions_bite_rate_dtr9, predictions_lifespan_dtr9, predictions_eggs_dtr9, predictions_gamma_dtr9,
						   predictions_bc, predictions_pea, predictions_mdr)

# 2B. Empirically fit TPCs for a, lf, and B in dtr12, other traits constant temperature TPCs
suit_empfluc_dtr12 <- R0eq(predictions_bite_rate_dtr12, predictions_lifespan_dtr12, predictions_eggs_dtr12, predictions_gamma_dtr12,
						   predictions_bc, predictions_pea, predictions_mdr)

# 3A. Rate summation on all traits assuming dtr9
suit_rs_traitsall_09 <- R0eq(predictions_bite_rate_rs9, predictions_lifespan_rs9, predictions_eggs_rs9, predictions_gamma_rs9,
							 predictions_bc_rs9, predictions_pea_rs9, predictions_mdr_rs9)

# 3B. Rate summation on all traits assuming dtr12
suit_rs_traitsall_12 <- R0eq(predictions_bite_rate_rs12, predictions_lifespan_rs12, predictions_eggs_rs12, predictions_gamma_rs12,
							 predictions_bc_rs12, predictions_pea_rs12, predictions_mdr_rs12)

# 4A. Rate summation on a, lf, and B assuming dtr9, other traits constant temperature TPCs
suit_rs_3traits_09 <- R0eq(predictions_bite_rate_rs9, predictions_lifespan_rs9, predictions_eggs_rs9, predictions_gamma_rs9,
						   predictions_bc, predictions_pea, predictions_mdr)

# 4B. Rate summation on a, lf, and B assuming dtr12, other traits constant temperature TPCs
suit_rs_3traits_12 <- R0eq(predictions_bite_rate_rs12, predictions_lifespan_rs12, predictions_eggs_rs12, predictions_gamma_rs12,
						   predictions_bc, predictions_pea, predictions_mdr)

# 5A. Rate summation on suit_const assuming dtr9
suit_rs_suit_09 <- RSCalcTempGrad(suit_const, LPtemps_dtr9, Temp.xs)

# 5A. Rate summation on suit_const assuming dtr12
suit_rs_suit_12 <- RSCalcTempGrad(suit_const, LPtemps_dtr12, Temp.xs)

# Save RS output for suitability (model 5)
save(suit_rs_suit_09, file = "data-processed/predictions_rs_suit_rs9.Rdata")
save(suit_rs_suit_12, file = "data-processed/predictions_rs_suit_rs12.Rdata")


##########
###### 6. Process rate summation calculation output for plotting
##########

#### Calculate TPC posterior summary data (means & quantiles)

predictions_bite_rate_rs9_summary <- calcPostQuants(predictions_bite_rate_rs9, "bite_rate_rs9", Temp.xs)
predictions_bite_rate_rs12_summary <- calcPostQuants(predictions_bite_rate_rs12, "bite_rate_rs12", Temp.xs)

predictions_lifespan_rs9_summary <- calcPostQuants(predictions_lifespan_rs9, "lifespan_rs9", Temp.xs)
predictions_lifespan_rs12_summary <- calcPostQuants(predictions_lifespan_rs12, "lifespan_rs12", Temp.xs)

predictions_eggs_rs9_summary <- calcPostQuants(predictions_eggs_rs9, "eggs_rs9", Temp.xs)
predictions_eggs_rs12_summary <- calcPostQuants(predictions_eggs_rs12, "eggs_rs12", Temp.xs)

predictions_gamma_rs9_summary <- calcPostQuants(predictions_gamma_rs9, "gamma_rs9", Temp.xs)
predictions_gamma_rs12_summary <- calcPostQuants(predictions_gamma_rs12, "gamma_rs12", Temp.xs)

predictions_bc_rs9_summary <- calcPostQuants(predictions_bc_rs9, "bc_rs9", Temp.xs)
predictions_bc_rs12_summary <- calcPostQuants(predictions_bc_rs12, "bc_rs12", Temp.xs)

predictions_eip50_rs9_summary <- calcPostQuants(predictions_eip50_rs9, "eip50_rs9", Temp.xs)
predictions_eip50_rs12_summary <- calcPostQuants(predictions_eip50_rs12, "eip50_rs12", Temp.xs)

predictions_pea_rs9_summary <- calcPostQuants(predictions_pea_rs9, "pea_rs9", Temp.xs)
predictions_pea_rs12_summary <- calcPostQuants(predictions_pea_rs12, "pea_rs12", Temp.xs)

predictions_mdr_rs9_summary <- calcPostQuants(predictions_mdr_rs9, "mdr_rs9", Temp.xs)
predictions_mdr_rs12_summary <- calcPostQuants(predictions_mdr_rs12, "mdr_rs12", Temp.xs)


### Get summary statistics of Tmin, Tmax, Topt, and Tbreadth

params_bite_rate_rs9_summary <- extractDerivedTPC(predictions_bite_rate_rs9, "bite_rate_rs9", Temp.xs)
params_bite_rate_rs12_summary <- extractDerivedTPC(predictions_bite_rate_rs12, "bite_rate_rs12", Temp.xs)

params_lifespan_rs9_summary <- extractDerivedTPC(predictions_lifespan_rs9, "lifespan_rs9", Temp.xs)
params_lifespan_rs12_summary <- extractDerivedTPC(predictions_lifespan_rs12, "lifespan_rs12", Temp.xs)

params_eggs_rs9_summary <- extractDerivedTPC(predictions_eggs_rs9, "eggs_rs9", Temp.xs)
params_eggs_rs12_summary <- extractDerivedTPC(predictions_eggs_rs12, "eggs_rs12", Temp.xs)


##########
###### 7. Join RS output with corresponding empirical TPCs
##########

### Biting Rate (a) 

### load empirical TPC summaries
predictions_bite_rate_constant_summary <- read_csv("data-processed/predictions_bite_rate_constant_summary.csv")
predictions_bite_rate_dtr9_summary <- read_csv("data-processed/predictions_bite_rate_dtr9_summary.csv")
predictions_bite_rate_dtr12_summary <- read_csv("data-processed/predictions_bite_rate_dtr12_summary.csv")

params_bite_rate_constant_summary <- read.csv("data-processed/params_bite_rate_constant_summary.csv")
params_bite_rate_dtr9_summary <- read.csv("data-processed/params_bite_rate_dtr9_summary.csv")
params_bite_rate_dtr12_summary <- read.csv("data-processed/params_bite_rate_dtr12_summary.csv")

### combine empirically measured and rate summation treatments (separately for dtr 9 and for dtr12) - plus constant for comparison
dtr9_rs9_bite_rate_predictions <- bind_rows(predictions_bite_rate_constant_summary, predictions_bite_rate_dtr9_summary, predictions_bite_rate_rs9_summary)
dtr12_rs12_bite_rate_predictions <- bind_rows(predictions_bite_rate_constant_summary, predictions_bite_rate_dtr12_summary, predictions_bite_rate_rs12_summary)

dtr9_rs9_bite_rate_params <- bind_rows(params_bite_rate_constant_summary, params_bite_rate_dtr9_summary, params_bite_rate_rs9_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
dtr12_rs12_bite_rate_params <- bind_rows(params_bite_rate_constant_summary, params_bite_rate_dtr12_summary, params_bite_rate_rs12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### Lifespan (lf) 

### load empirical TPC summaries
predictions_lifespan_constant_summary <- read_csv("data-processed/predictions_lifespan_constant_summary.csv")
predictions_lifespan_dtr9_summary <- read_csv("data-processed/predictions_lifespan_dtr9_summary.csv")
predictions_lifespan_dtr12_summary <- read_csv("data-processed/predictions_lifespan_dtr12_summary.csv")

params_lifespan_constant_summary <- read.csv("data-processed/params_lifespan_constant_summary.csv")
params_lifespan_dtr9_summary <- read.csv("data-processed/params_lifespan_dtr9_summary.csv")
params_lifespan_dtr12_summary <- read.csv("data-processed/params_lifespan_dtr12_summary.csv")

### combine empirically measured and rate summation treatments (separately for dtr 9 and for dtr12) - plus constant for comparison
dtr9_rs9_lifespan_predictions <- bind_rows(predictions_lifespan_constant_summary, predictions_lifespan_dtr9_summary, predictions_lifespan_rs9_summary)
dtr12_rs12_lifespan_predictions <- bind_rows(predictions_lifespan_constant_summary, predictions_lifespan_dtr12_summary, predictions_lifespan_rs12_summary)

dtr9_rs9_lifespan_params <- bind_rows(params_lifespan_constant_summary, params_lifespan_dtr9_summary, params_lifespan_rs9_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
dtr12_rs12_lifespan_params <- bind_rows(params_lifespan_constant_summary, params_lifespan_dtr12_summary, params_lifespan_rs12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### Lifetime eggs (B) 

### load empirical TPC summaries
predictions_eggs_constant_summary <- read_csv("data-processed/predictions_eggs_constant_summary.csv")
predictions_eggs_dtr9_summary <- read_csv("data-processed/predictions_eggs_dtr9_summary.csv")
predictions_eggs_dtr12_summary <- read_csv("data-processed/predictions_eggs_dtr12_summary.csv")

params_eggs_constant_summary <- read.csv("data-processed/params_eggs_constant_summary.csv")
params_eggs_dtr9_summary <- read.csv("data-processed/params_eggs_dtr9_summary.csv")
params_eggs_dtr12_summary <- read.csv("data-processed/params_eggs_dtr12_summary.csv")

### combine empirically measured and rate summation treatments (separately for dtr 9 and for dtr12) - plus constant for comparison
dtr9_rs9_eggs_predictions <- bind_rows(predictions_eggs_constant_summary, predictions_eggs_dtr9_summary, predictions_eggs_rs9_summary)
dtr12_rs12_eggs_predictions <- bind_rows(predictions_eggs_constant_summary, predictions_eggs_dtr12_summary, predictions_eggs_rs12_summary)

dtr9_rs9_eggs_params <- bind_rows(params_eggs_constant_summary, params_eggs_dtr9_summary, params_eggs_rs9_summary) %>%
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
dtr12_rs12_eggs_params <- bind_rows(params_eggs_constant_summary, params_eggs_dtr12_summary, params_eggs_rs12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

##########
###### 8. Calculate percent change compared to constant temperature model
##########

# Calculate % reduction in maximum value at Topt for each fluctuation model vs constant temperature model
1 - max(predictions_bite_rate_dtr9_summary$median)/max(predictions_bite_rate_constant_summary$median) # 25.1%
1 - max(predictions_bite_rate_dtr12_summary$median)/max(predictions_bite_rate_constant_summary$median) # 23.5%
1 - max(predictions_bite_rate_rs9_summary$median)/max(predictions_bite_rate_constant_summary$median) # 5.1%
1 - max(predictions_bite_rate_rs12_summary$median)/max(predictions_bite_rate_constant_summary$median) # 8.9%

1 - max(predictions_lifespan_dtr9_summary$median)/max(predictions_lifespan_constant_summary$median) # 2.7%
1 - max(predictions_lifespan_dtr12_summary$median)/max(predictions_lifespan_constant_summary$median) # -10.7% - actually higher!
1 - max(predictions_lifespan_rs9_summary$median)/max(predictions_lifespan_constant_summary$median) # 3.0%
1 - max(predictions_lifespan_rs12_summary$median)/max(predictions_lifespan_constant_summary$median) # 5.2%

1 - max(predictions_eggs_dtr9_summary$median)/max(predictions_eggs_constant_summary$median) # 7.9%
1 - max(predictions_eggs_dtr12_summary$median)/max(predictions_eggs_constant_summary$median) # 14.8%
1 - max(predictions_eggs_rs9_summary$median)/max(predictions_eggs_constant_summary$median) # 6.5%
1 - max(predictions_eggs_rs12_summary$median)/max(predictions_eggs_constant_summary$median) # 11.5%

# Calculate % reduction in empirical fluctuating model vs rate summation fluctuating model
1 - max(predictions_bite_rate_dtr9_summary$median)/max(predictions_bite_rate_rs9_summary$median) # 21.1%
1 - max(predictions_bite_rate_dtr12_summary$median)/max(predictions_bite_rate_rs12_summary$median) # 16.1%

1 - max(predictions_lifespan_dtr9_summary$median)/max(predictions_lifespan_rs9_summary$median) # -0.2%
1 - max(predictions_lifespan_dtr12_summary$median)/max(predictions_lifespan_rs12_summary$median) # -16.8% - actually higher!

1 - max(predictions_eggs_dtr9_summary$median)/max(predictions_eggs_rs9_summary$median) # 1.5%
1 - max(predictions_eggs_dtr12_summary$median)/max(predictions_eggs_rs12_summary$median) # 3.7%


##########
###### 9. Manuscript Figure 3
##########

### Biting Rate (a) panels

bite_rate_rs9_plot <- dtr9_rs9_bite_rate_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = -13.5, y = 0.25, size = 7, angle = 90, label = "DTR 9") +
	geom_text(x = 25, y = 0.6, size = 7, label = "Bite rate (a)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab(parse(text = "Bite~rate~(day^-1)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2.5, 1, 1, 2.5), "lines"), axis.title.y = element_text(vjust = -2)) +
	coord_cartesian(clip = "off") +
	annotate("text", x = 0, y = 0.53, label = "A", size = 5)

params_plot_bite_rate_rs9 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr9_rs9_bite_rate_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "DTR 9", "RS 9")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "D", size = 5)

bite_rate_rs12_plot <- dtr12_rs12_bite_rate_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = -13.5, y = 0.25, size = 7, angle = 90, label = "DTR 12") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp12, ct_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab(parse(text = "Bite~rate~(day^-1)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(1, 1, 1, 2.5), "lines"), axis.title.y = element_text(vjust = -2)) +
	coord_cartesian(clip = "off") +
	annotate("text", x = 0, y = 0.53, label = "G", size = 5)

params_plot_bite_rate_rs12 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr12_rs12_bite_rate_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "DTR 12", "RS 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "J", size = 5)

bite_rate9_plot <- bite_rate_rs9_plot / params_plot_bite_rate_rs9 + plot_layout(heights = c(3, 0.5))
bite_rate12_plot <- bite_rate_rs12_plot / params_plot_bite_rate_rs12 + plot_layout(heights = c(3, 0.5))

### Lifespan (lf) panels

lifespan_rs9_plot <- dtr9_rs9_lifespan_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = 25, y = 44, size = 7, label = "Lifespan (lf)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fitted", "Rate summation")) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -1)) +
	coord_cartesian(clip = "off") +
	theme(plot.margin = unit(c(2.5, 1, 1, 1), "lines")) +
	annotate("text", x = 0, y = 39, label = "B", size = 5)

params_plot_lifespan_rs9 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr9_rs9_lifespan_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "DTR 9", "RS 9")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "E", size = 5)

lifespan_rs12_plot <- dtr12_rs12_lifespan_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp12, ct_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fitted", "Rate summation")) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -1)) + 
	annotate("text", x = 0, y = 44, label = "H", size = 5)

params_plot_lifespan_rs12 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr12_rs12_lifespan_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "DTR 12", "RS 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "K", size = 5)

lifespan9_plot <- lifespan_rs9_plot / params_plot_lifespan_rs9 + plot_layout(heights = c(3, 0.5))
lifespan12_plot <- lifespan_rs12_plot / params_plot_lifespan_rs12 + plot_layout(heights = c(3, 0.5))

### Lifetime eggs (B) panels

eggs_rs9_plot <- dtr9_rs9_eggs_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = 25, y = 495, size = 7, label = "Lifetime eggs (B)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab("Lifetime eggs") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -2)) +
	coord_cartesian(clip = "off") +
	theme(plot.margin = unit(c(2.5, 1, 1, 1), "lines")) +
	annotate("text", x = 0, y = 440, label = "C", size = 5)

params_plot_eggs_rs9 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr9_rs9_eggs_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "DTR 9", "RS 9")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "F", size = 5)

eggs_rs12_plot <- dtr12_rs12_eggs_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp12, ct_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab("Lifetime eggs") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -2)) + 
	annotate("text", x = 0, y = 440, label = "I", size = 5)

params_plot_eggs_rs12 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr12_rs12_eggs_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "DTR 12", "RS 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "L", size = 5)

eggs9_plot <- eggs_rs9_plot / params_plot_eggs_rs9 + plot_layout(heights = c(3, 0.5))
eggs12_plot <- eggs_rs12_plot / params_plot_eggs_rs12 + plot_layout(heights = c(3, 0.5))

### Combine panels together

Fig3_plots <- wrap_plots(bite_rate9_plot, lifespan9_plot, eggs9_plot,
						 bite_rate12_plot, lifespan12_plot, eggs12_plot, ncol = 3)
ggsave('figures/Fig3_RS_trait_TPCs.pdf', Fig3_plots, width = 15, height = 12)


##########
###### 10. Process Suitability output
##########

# Calculate quantiles for plotting
suit_const_summary <- calcPostQuants(suit_const, "suit_const", Temp.xs)

suit_empfluc_dtr09_summary <- calcPostQuants(suit_empfluc_dtr09, "suit_empfluc_dtr09", Temp.xs)
suit_empfluc_dtr12_summary <- calcPostQuants(suit_empfluc_dtr12, "suit_empfluc_dtr12", Temp.xs)

suit_rs_traitsall_09_summary <- calcPostQuants(suit_rs_traitsall_09, "suit_rs_traitsall_09", Temp.xs)
suit_rs_traitsall_12_summary <- calcPostQuants(suit_rs_traitsall_12, "suit_rs_traitsall_12", Temp.xs)

suit_rs_3traits_09_summary <- calcPostQuants(suit_rs_3traits_09, "suit_rs_3traits_09", Temp.xs)
suit_rs_3traits_12_summary <- calcPostQuants(suit_rs_3traits_12, "suit_rs_3traits_12", Temp.xs)

suit_rs_suit_09_summary <- calcPostQuants(suit_rs_suit_09, "suit_rs_suit_09", Temp.xs)
suit_rs_suit_12_summary <- calcPostQuants(suit_rs_suit_12, "suit_rs_suit_12", Temp.xs)

# Calculate maximum predicted R0 value for each model
max(suit_const_summary$median)
max(suit_empfluc_dtr09_summary$median)
max(suit_empfluc_dtr12_summary$median)
max(suit_rs_traitsall_09_summary$median)
max(suit_rs_traitsall_12_summary$median)
max(suit_rs_3traits_09_summary$median)
max(suit_rs_3traits_12_summary$median)
max(suit_rs_suit_09_summary$median)
max(suit_rs_suit_12_summary$median)

# Calculate % reduction in maximum R0 at Topt for each fluctuation model vs constant temperature model
1 - max(suit_empfluc_dtr09_summary$median)/max(suit_const_summary$median) # 32.0%
1 - max(suit_empfluc_dtr12_summary$median)/max(suit_const_summary$median) # 33.8%
1 - max(suit_rs_traitsall_09_summary$median)/max(suit_const_summary$median) # 18.1%
1 - max(suit_rs_traitsall_12_summary$median)/max(suit_const_summary$median) # 30.6%
1 - max(suit_rs_3traits_09_summary$median)/max(suit_const_summary$median) # 10.0%
1 - max(suit_rs_3traits_12_summary$median)/max(suit_const_summary$median) # 17.1%
1 - max(suit_rs_suit_09_summary$median)/max(suit_const_summary$median) # 19.9 %
1 - max(suit_rs_suit_12_summary$median)/max(suit_const_summary$median) # 32.0 %

# Calculate the maximum value of the median S(T) from the constant temperature model
suitability_scale <- max(suit_const_summary$median)

# Get summary statistics of Tmin, Tmax, Topt, and Tbreadth 
params_suit_const <- extractDerivedTPC(suit_const, "suit_const", Temp.xs)

params_suit_empfluc_dtr09 <- extractDerivedTPC(suit_empfluc_dtr09, "suit_empfluc_dtr09", Temp.xs)
params_suit_empfluc_dtr12 <- extractDerivedTPC(suit_empfluc_dtr12, "suit_empfluc_dtr12", Temp.xs)

params_suit_rs_traitsall_09 <- extractDerivedTPC(suit_rs_traitsall_09, "suit_rs_traitsall_09", Temp.xs)
params_suit_rs_traitsall_12 <- extractDerivedTPC(suit_rs_traitsall_12, "suit_rs_traitsall_12", Temp.xs)

params_suit_rs_3traits_09 <- extractDerivedTPC(suit_rs_3traits_09, "suit_rs_3traits_09", Temp.xs)
params_suit_rs_3traits_12 <- extractDerivedTPC(suit_rs_3traits_12, "suit_rs_3traits_12", Temp.xs)

params_suit_rs_suit_09 <- extractDerivedTPC(suit_rs_suit_09, "suit_rs_suit_09", Temp.xs)
params_suit_rs_suit_12 <- extractDerivedTPC(suit_rs_suit_12, "suit_rs_suit_12", Temp.xs)


##########
###### 11. Manuscript Figure 4
##########

##### Summary of different versions of suitability calculations
#
# Five different versions of S(T) for specific pairwise comparisons:
# 	1) suit_const:			constant temperature TPCs
#	2) suit_empfluc:		empirically fit TPCs for a, lf, and B in fluctuating environments (other traits varying based on constant temperature TPCs)
#	3) suit_rs_traitsall:	rate summation on all trait TPCs
#	4) suit_rs_3traits:		rate summation on a, lf, and B (other traits varying based on constant temperature TPCs)
#	5) suit_rs_suit: 		rate summation on S(T) TPC
#
# Pairwise comparisons and what they mean biologically or statistically:
#	1 vs 2 - Does fluctuating temperature (empirically) affect S(T)?
#	1 vs 3 - Does incorporating rate summation (on traits) affect predicted S(T)?
#	2 vs 4 - Does rate summation accurately predict the effect of fluctuating temperature on S(T)?
#	3 vs 5 - Does performing rate summation on traits vs. on S(T) TPC affect predicted S(T)?


################################ Combining treatments for different comparisons

#### Comparison #1: combine all empirical treatments
emp_suit_predictions <- bind_rows(suit_const_summary, suit_empfluc_dtr09_summary, suit_empfluc_dtr12_summary)
emp_suit_params <- bind_rows(params_suit_const, params_suit_empfluc_dtr09, params_suit_empfluc_dtr12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

#### Comparison #2: combine rate summation on 3 traits and empirical dtr treatments
rsaccuracy_suit_predictions <- bind_rows(suit_rs_3traits_09_summary, suit_rs_3traits_12_summary, suit_empfluc_dtr09_summary, suit_empfluc_dtr12_summary)
rsaccuracy_suit_params <- bind_rows(params_suit_rs_3traits_09, params_suit_rs_3traits_12, params_suit_empfluc_dtr09, params_suit_empfluc_dtr12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

#### Comparison #3: combine rate summation on all traits and rate summation on suitability
rslevel_suit_predictions <- bind_rows(suit_rs_traitsall_09_summary, suit_rs_traitsall_12_summary, suit_rs_suit_09_summary, suit_rs_suit_12_summary)
rslevel_suit_params <- bind_rows(params_suit_rs_traitsall_09, params_suit_rs_traitsall_12, params_suit_rs_suit_09, params_suit_rs_suit_12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

# #### Comparison #4: combine rate summation on all traits and constant
rseffect_suit_predictions <- bind_rows(suit_const_summary, suit_rs_traitsall_09_summary, suit_rs_traitsall_12_summary)
rseffect_suit_params <- bind_rows(params_suit_const, params_suit_rs_traitsall_09, params_suit_rs_traitsall_12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

#### Combine all suitability calculations - for supplemental figure
all_suit_predictions <- bind_rows(suit_const_summary, suit_empfluc_dtr09_summary, suit_empfluc_dtr12_summary, suit_rs_3traits_09_summary, suit_rs_3traits_12_summary,
								  suit_rs_traitsall_09_summary, suit_rs_traitsall_12_summary, suit_rs_suit_09_summary, suit_rs_suit_12_summary)


#### Comparison #1: combine all empirical treatments

emp_suit_plot <- emp_suit_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("1: Constant", "2: Empirical DTR 9", "2: Empirical DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("1: Constant", "2: Empirical DTR 9", "2: Empirical DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1), labels = c("1: Constant", "2: Empirical DTR 9", "2: Empirical DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 67, label = "A", size = 5)

params_emp_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = emp_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "E", size = 5)

emp_suit_plot_combined <- emp_suit_plot / params_emp_suit_plot + plot_layout(heights = c(3, 0.75))

#### Comparison #2: combine rate summation on 3 traits and empirical dtr treatments

rsaccuracy_suit_plot <- rsaccuracy_suit_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_emp09, c_emp12, c_rstraits09, c_rstraits12), labels = c("2: Empirical DTR 9", "2: Empirical DTR 12", "3: RS on 3 traits - DTR 9", "3: RS on 3 traits - DTR 12")) +
	scale_fill_manual(values = c(ct_emp09, ct_emp12, ct_rstraits09, ct_rstraits12), labels = c("2: Empirical DTR 9", "2: Empirical DTR 12", "3: RS on 3 traits - DTR 9", "3: RS on 3 traits - DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 2, 2), labels = c("2: Empirical DTR 9", "2: Empirical DTR 12", "3: RS on 3 traits - DTR 9", "3: RS on 3 traits - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 60, label = "B", size = 5)

params_rsaccuracy_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rsaccuracy_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_emp09, c_emp12, c_rstraits09, c_rstraits12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "F", size = 5)

rsaccuracy_suit_plot_combined <- rsaccuracy_suit_plot / params_rsaccuracy_suit_plot + plot_layout(heights = c(3, 0.75))

#### Comparison #3: combine rate summation on all traits and rate summation on suitability

rslevel_suit_plot <- rslevel_suit_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_rssuit09, c_rssuit12, c_rstraits09, c_rstraits12), labels = c("5: RS on S(T) - DTR 9", "5: RS on S(T) - DTR 12", "4: RS on all traits - DTR 9", "4: RS on all traits - DTR 12")) +
	scale_fill_manual(values = c(ct_rssuit09, ct_rssuit12, ct_rstraits09, ct_rstraits12), labels = c("5: RS on S(T) - DTR 9", "5: RS on S(T) - DTR 12", "4: RS on all traits - DTR 9", "4: RS on all traits - DTR 12")) +
	scale_linetype_manual(values = c(3, 3, 2, 2), labels = c("5: RS on S(T) - DTR 9", "5: RS on S(T) - DTR 12", "4: RS on all traits - DTR 9", "4: RS on all traits - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 54, label = "C", size = 5)

params_rslevel_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rslevel_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_rssuit09, c_rssuit12, c_rstraits09, c_rstraits12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "G", size = 5)

rslevel_suit_plot_combined <- rslevel_suit_plot / params_rslevel_suit_plot + plot_layout(heights = c(3, 0.75))

#### Panel D: Lines for all models, without CIs

rseffect_suit_plot <- rseffect_suit_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("1: Constant", "4: RS all traits - DTR 9", "4: RS all traits - DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("1: Constant", "4: RS all traits - DTR 9", "4: RS all traits - DTR 12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("1: Constant", "4: RS all traits - DTR 9", "4: RS all traits - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 67, label = "D", size = 5)

params_rseffect_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rseffect_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "H", size = 5)

rseffect_suit_plot_combined <- rseffect_suit_plot / params_rseffect_suit_plot + plot_layout(heights = c(3, 0.75))

#### Combining plots

Fig4_plots <- wrap_plots(emp_suit_plot_combined, rsaccuracy_suit_plot_combined, rslevel_suit_plot_combined, rseffect_suit_plot_combined, nrow = 2)
ggsave('figures/Fig4_Suitability.pdf', Fig4_plots, width = 12, height = 12)


# All models plotted together
all_suit_plot <- all_suit_predictions %>% 
	ggplot() +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12, c_emp09, c_emp12, c_emp09, c_emp12, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12", "RS 3 - DTR 9", "RS 3 - DTR 12", "RS all - DTR 9", "RS all - DTR 12", "RS S(T) - DTR 9", "RS S(T) - DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1, 2, 2, 2, 2, 3, 3), labels = c("Constant", "DTR 9", "DTR 12", "RS 3 - DTR 9", "RS 3 - DTR 12", "RS all - DTR 9", "RS all - DTR 12", "RS S(T) - DTR 9", "RS S(T) - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.15, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 67, label = "A", size = 5)


##########
###### 12. Processing suitability output for mapping
##########

head(suit_const_summary)

# Scale each model's lower CI contour to itself get relative R0/S(T) for mapping
suitabilityOutputForMapping_lowerCI <- data.frame(temp = suit_const_summary$temperature,
												  constant = suit_const_summary$lowerCI/max(suit_const_summary$lowerCI),
												  empDTR12 = suit_empfluc_dtr12_summary$lowerCI/max(suit_empfluc_dtr12_summary$lowerCI),
												  rsAllTraitsDTR12 = suit_rs_traitsall_12_summary$lowerCI/max(suit_rs_traitsall_12_summary$lowerCI),
												  rsSuitDTR12 = suit_rs_suit_12_summary$lowerCI/max(suit_rs_suit_12_summary$lowerCI))
#Check output
plot(constant ~ temp, data = suitabilityOutputForMapping_lowerCI, type = "l")
lines(empDTR12 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "darkorange")
lines(rsAllTraitsDTR12 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "blue")
lines(rsSuitDTR12 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "purple")

write.csv(suitabilityOutputForMapping_lowerCI, "data-processed/RateSummationProjectR0forMapping_lowerCI.csv")
