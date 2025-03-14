## Rate Summation Project
## Script for sensitivity and uncertainty analyses
## Written by Marta Shocket in 2022-24
## Based on code originally written by Marta Shocket in 2018 from code written by Leah Johnson in 2015
## Previously updated by Erin Mordecai in 2018 and Kerri Miazgowicz in 2019-2020

##	Table of contents:
##
##	1. Set-up workspace
##  2. Define S(T)/R0 function and JAGS Settings
##  3. Generate hourly temperatures with the Parton-Logan model
##  4. Calculate trait means across temperature gradient
##  5. Sensitivity analysis #1 - partial derivatives
##  6. Sensitivity analysis #2 - holding single parameters constant
##  7. Uncertainty analysis


##########
###### 1. Set-up workspace
##########

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
load("saved-posteriors/bc.Rdata")

Temp.xs <- seq(0, 45, 0.1) # temperature gradient

##### Extract matrices of predicted trait values from JAGS models
# NOTE: These need to be matrices *NOT* be data frames to use colMeans() & the SA/UA functions below,
# unlike in the rate summation calculations that use dplyr and must take data frames
predictions_bite_rate_constant <- model_bite_rate_constant$BUGSoutput$sims.list$z.trait.mu.pred
predictions_bite_rate_dtr9 <- model_bite_rate_dtr9$BUGSoutput$sims.list$z.trait.mu.pred
predictions_bite_rate_dtr12 <- model_bite_rate_dtr12$BUGSoutput$sims.list$z.trait.mu.pred

predictions_lifespan_constant <- model_lifespan_constant$BUGSoutput$sims.list$z.trait.mu.pred
predictions_lifespan_dtr9 <- model_lifespan_dtr9$BUGSoutput$sims.list$z.trait.mu.pred
predictions_lifespan_dtr12 <- model_lifespan_dtr12$BUGSoutput$sims.list$z.trait.mu.pred

predictions_eggs_constant <- model_eggs_constant$BUGSoutput$sims.list$z.trait.mu.pred
predictions_eggs_dtr9 <- model_eggs_dtr9$BUGSoutput$sims.list$z.trait.mu.pred
predictions_eggs_dtr12 <- model_eggs_dtr12$BUGSoutput$sims.list$z.trait.mu.pred

predictions_gamma_constant <- model_gamma_constant$BUGSoutput$sims.list$z.trait.mu.pred
predictions_gamma_dtr9 <- model_gamma_dtr9$BUGSoutput$sims.list$z.trait.mu.pred
predictions_gamma_dtr12 <- model_gamma_dtr12$BUGSoutput$sims.list$z.trait.mu.pred

predictions_bc <- model_bc$BUGSoutput$sims.list$z.trait.mu.pred
predictions_pea <- model_pea$BUGSoutput$sims.list$z.trait.mu.pred
predictions_mdr <- model_mdr$BUGSoutput$sims.list$z.trait.mu.pred

##### Load saved data frames of predicted trait values from RS calculations
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

##### Convert data frames of predicted trait values from RS calculations to matrices
predictions_bite_rate_rs9 <- as.matrix(predictions_bite_rate_rs9)
predictions_bite_rate_rs12 <- as.matrix(predictions_bite_rate_rs12)

predictions_lifespan_rs9 <- as.matrix(predictions_lifespan_rs9)
predictions_lifespan_rs12 <- as.matrix(predictions_lifespan_rs12)

predictions_eggs_rs9 <- as.matrix(predictions_eggs_rs9)
predictions_eggs_rs12 <- as.matrix(predictions_eggs_rs12)

predictions_gamma_rs9 <- as.matrix(predictions_gamma_rs9)
predictions_gamma_rs12 <- as.matrix(predictions_gamma_rs12)

predictions_bc_rs9 <- as.matrix(predictions_bc_rs9)
predictions_bc_rs12 <- as.matrix(predictions_bc_rs12)

predictions_eip50_rs9 <- as.matrix(predictions_eip50_rs9)
predictions_eip50_rs12 <- as.matrix(predictions_eip50_rs12)

predictions_pea_rs9 <- as.matrix(predictions_pea_rs9)
predictions_pea_rs12 <- as.matrix(predictions_pea_rs12)

predictions_mdr_rs9 <- as.matrix(predictions_mdr_rs9)
predictions_mdr_rs12 <- as.matrix(predictions_mdr_rs12)


##########
###### 2. Define S(T)/R0 function and JAGS settings
##########

##### Small constants to keep values from being numerically zero
# For S(T)/RO equation, assume minimum survival time is half an hour
ec.lf = 1/48
# A generic constant for all parameters in sensitivity analysis #1
ec <- 0.000001

# R0 formulation where: 1) mosquito density (M) depends on lifetime fecundity; 2) gamma (y) is substituted for exp^(-1/(lf*PDR))expression
R0eq = function(a, bc, lf, y, B, pEA, MDR) {
	M = B * pEA * MDR * (lf+ec.lf)
	R0 = (a^2 * bc * y * M * (lf+ec.lf) )^0.5
	return(R0)
}

# Temperature levels and # MCMC steps
N.Temp.xs <- length(Temp.xs)
nMCMC <- 7500


##########
###### 3. Generate hourly temperatures with the Parton-Logan model
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
###### 4. Calculate trait means across temperature gradient
##########

m_bite_rate_constant <- colMeans(predictions_bite_rate_constant)
m_bite_rate_dtr9 <- colMeans(predictions_bite_rate_dtr9)
m_bite_rate_dtr12 <- colMeans(predictions_bite_rate_dtr12)
m_bite_rate_rs9 <- colMeans(predictions_bite_rate_rs9)
m_bite_rate_rs12 <- colMeans(predictions_bite_rate_rs12)

m_lifespan_constant <- colMeans(predictions_lifespan_constant)
m_lifespan_dtr9 <- colMeans(predictions_lifespan_dtr9)
m_lifespan_dtr12 <- colMeans(predictions_lifespan_dtr12)
m_lifespan_rs9 <- colMeans(predictions_lifespan_rs9)
m_lifespan_rs12 <- colMeans(predictions_lifespan_rs12)

m_eggs_constant <- colMeans(predictions_eggs_constant)
m_eggs_dtr9 <- colMeans(predictions_eggs_dtr9)
m_eggs_dtr12 <- colMeans(predictions_eggs_dtr12)
m_eggs_rs9 <- colMeans(predictions_eggs_rs9)
m_eggs_rs12 <- colMeans(predictions_eggs_rs12)

m_gamma_constant <- colMeans(predictions_gamma_constant)
m_gamma_dtr9 <- colMeans(predictions_gamma_dtr9)
m_gamma_dtr12 <- colMeans(predictions_gamma_dtr12)
m_gamma_rs9 <- colMeans(predictions_gamma_rs9)
m_gamma_rs12 <- colMeans(predictions_gamma_rs12)

m_bc <- colMeans(predictions_bc)
m_bc_rs9 <- colMeans(predictions_bc_rs9)
m_bc_rs12 <- colMeans(predictions_bc_rs12)

m_pea <- colMeans(predictions_pea)
m_pea_rs9 <- colMeans(predictions_pea_rs9)
m_pea_rs12 <- colMeans(predictions_pea_rs12)

m_mdr <- colMeans(predictions_mdr)
m_mdr_rs9 <- colMeans(predictions_mdr_rs9)
m_mdr_rs12 <- colMeans(predictions_mdr_rs12)


##########
###### 5. Sensitivity Analysis #1 - partial derivatives
##########

# NOTE: This method can only be applied to the empirical models fit in JAGS (not the rate summation models)
# because it uses derivatives based on the quadratic/Briere functions the fitted Tmin, Tmax and q parameters

###################################### Constant

# Calculate sensitivity using partial derivatives
SA1_constant <- SensitivityAnalysis1_pd(model_bite_rate_constant, model_bc, model_lifespan_constant, model_gamma_constant, model_eggs_constant, model_pea, model_mdr,
										   m_bite_rate_constant, m_bc, m_lifespan_constant, m_gamma_constant, m_eggs_constant, m_pea, m_mdr)

# Get sensitivity posteriors for each parameter and summarize them
dR0da_constant		<- calcPostQuants(as.data.frame(SA1_constant[[1]]), "SA1_constant_a", Temp.xs)
dR0dbc_constant		<- calcPostQuants(as.data.frame(SA1_constant[[2]]), "SA1_constant_bc", Temp.xs)
dR0dlf_constant		<- calcPostQuants(as.data.frame(SA1_constant[[3]]), "SA1_constant_lf", Temp.xs)
dR0dgamma_constant	<- calcPostQuants(as.data.frame(SA1_constant[[4]]), "SA1_constant_gamma", Temp.xs)
dR0dB_constant		<- calcPostQuants(as.data.frame(SA1_constant[[5]]), "SA1_constant_B", Temp.xs)
dR0dpEA_constant 	<- calcPostQuants(as.data.frame(SA1_constant[[6]]), "SA1_constant_pEA", Temp.xs)
dR0dMDR_constant 	<- calcPostQuants(as.data.frame(SA1_constant[[7]]), "SA1_constant_MDR", Temp.xs)
dR0dT_constant		<- calcPostQuants(as.data.frame(SA1_constant[[8]]), "SA1_constant_R0", Temp.xs)

###################################### Empirical DTR 9 

# Calculate sensitivity using partial derivatives
SA1_dtr9 <- SensitivityAnalysis1_pd(model_bite_rate_dtr9, model_bc, model_lifespan_dtr9, model_gamma_dtr9, model_eggs_dtr9, model_pea, model_mdr,
										m_bite_rate_dtr9, m_bc, m_lifespan_dtr9, m_gamma_dtr9, m_eggs_dtr9, m_pea, m_mdr)

# Get sensitivity posteriors for each parameter and summarize them
dR0da_dtr9		<- calcPostQuants(as.data.frame(SA1_dtr9[[1]]), "SA1_dtr9_a", Temp.xs)
dR0dbc_dtr9		<- calcPostQuants(as.data.frame(SA1_dtr9[[2]]), "SA1_dtr9_bc", Temp.xs)
dR0dlf_dtr9		<- calcPostQuants(as.data.frame(SA1_dtr9[[3]]), "SA1_dtr9_lf", Temp.xs)
dR0dgamma_dtr9	<- calcPostQuants(as.data.frame(SA1_dtr9[[4]]), "SA1_dtr9_gamma", Temp.xs)
dR0dB_dtr9		<- calcPostQuants(as.data.frame(SA1_dtr9[[5]]), "SA1_dtr9_B", Temp.xs)
dR0dpEA_dtr9 	<- calcPostQuants(as.data.frame(SA1_dtr9[[6]]), "SA1_dtr9_pEA", Temp.xs)
dR0dMDR_dtr9 	<- calcPostQuants(as.data.frame(SA1_dtr9[[7]]), "SA1_dtr9_MDR", Temp.xs)
dR0dT_dtr9		<- calcPostQuants(as.data.frame(SA1_dtr9[[8]]), "SA1_dtr9_R0", Temp.xs)

###################################### Empirical DTR 12

# Calculate sensitivity using partial derivatives
SA1_dtr12 <- SensitivityAnalysis1_pd(model_bite_rate_dtr12, model_bc, model_lifespan_dtr12, model_gamma_dtr12, model_eggs_dtr12, model_pea, model_mdr,
										m_bite_rate_dtr12, m_bc, m_lifespan_dtr12, m_gamma_dtr12, m_eggs_dtr12, m_pea, m_mdr)

# Get sensitivity posteriors for each parameter and summarize them
dR0da_dtr12		<- calcPostQuants(as.data.frame(SA1_dtr12[[1]]), "SA1_dtr12_a", Temp.xs)
dR0dbc_dtr12	<- calcPostQuants(as.data.frame(SA1_dtr12[[2]]), "SA1_dtr12_bc", Temp.xs)
dR0dlf_dtr12	<- calcPostQuants(as.data.frame(SA1_dtr12[[3]]), "SA1_dtr12_lf", Temp.xs)
dR0dgamma_dtr12	<- calcPostQuants(as.data.frame(SA1_dtr12[[4]]), "SA1_dtr12_gamma", Temp.xs)
dR0dB_dtr12		<- calcPostQuants(as.data.frame(SA1_dtr12[[5]]), "SA1_dtr12_B", Temp.xs)
dR0dpEA_dtr12 	<- calcPostQuants(as.data.frame(SA1_dtr12[[6]]), "SA1_dtr12_pEA", Temp.xs)
dR0dMDR_dtr12 	<- calcPostQuants(as.data.frame(SA1_dtr12[[7]]), "SA1_dtr12_MDR", Temp.xs)
dR0dT_dtr12		<- calcPostQuants(as.data.frame(SA1_dtr12[[8]]), "SA1_dtr12_R0", Temp.xs)


##########
###### 6. Sensitivity Analysis #2 - holding single parameters constant
##########

# NOTE: This method can only be applied to all models, including the rate summation models

###################################### 1. Constant

SA2_constant <- SensitivityAnalysis2_hspc(predictions_bite_rate_constant, predictions_bc, predictions_lifespan_constant, predictions_gamma_constant, predictions_eggs_constant, predictions_pea, predictions_mdr)

# Get posterior quantiles for plotting
SA2_constant_a_out <- calcPostQuants(as.data.frame(SA2_constant[[1]]), "SA2_constant_a", Temp.xs)
SA2_constant_bc_out <- calcPostQuants(as.data.frame(SA2_constant[[2]]), "SA2_constant_bc", Temp.xs)
SA2_constant_lf_out <- calcPostQuants(as.data.frame(SA2_constant[[3]]), "SA2_constant_lf", Temp.xs)
SA2_constant_gamma_out <- calcPostQuants(as.data.frame(SA2_constant[[4]]), "SA2_constant_gamma", Temp.xs)
SA2_constant_B_out <- calcPostQuants(as.data.frame(SA2_constant[[5]]), "SA2_constant_B", Temp.xs)
SA2_constant_pEA_out <- calcPostQuants(as.data.frame(SA2_constant[[6]]), "SA2_constant_pEA", Temp.xs)
SA2_constant_MDR_out <- calcPostQuants(as.data.frame(SA2_constant[[7]]), "SA2_constant_MDR", Temp.xs)
SA2_constant_all_out <- calcPostQuants(as.data.frame(SA2_constant[[8]]), "SA2_constant_all", Temp.xs)

###################################### 2A. Empirical DTR 9

SA2_dtr9 <- SensitivityAnalysis2_hspc(predictions_bite_rate_dtr9, predictions_bc, predictions_lifespan_dtr9, predictions_gamma_dtr9, predictions_eggs_dtr9, predictions_pea, predictions_mdr)

# Get posterior quantiles for plotting
SA2_dtr9_a_out <- calcPostQuants(as.data.frame(SA2_dtr9[[1]]), "SA2_dtr9_a", Temp.xs)
SA2_dtr9_bc_out <- calcPostQuants(as.data.frame(SA2_dtr9[[2]]), "SA2_dtr9_bc", Temp.xs)
SA2_dtr9_lf_out <- calcPostQuants(as.data.frame(SA2_dtr9[[3]]), "SA2_dtr9_lf", Temp.xs)
SA2_dtr9_gamma_out <- calcPostQuants(as.data.frame(SA2_dtr9[[4]]), "SA2_dtr9_gamma", Temp.xs)
SA2_dtr9_B_out <- calcPostQuants(as.data.frame(SA2_dtr9[[5]]), "SA2_dtr9_B", Temp.xs)
SA2_dtr9_pEA_out <- calcPostQuants(as.data.frame(SA2_dtr9[[6]]), "SA2_dtr9_pEA", Temp.xs)
SA2_dtr9_MDR_out <- calcPostQuants(as.data.frame(SA2_dtr9[[7]]), "SA2_dtr9_MDR", Temp.xs)
SA2_dtr9_all_out <- calcPostQuants(as.data.frame(SA2_dtr9[[8]]), "SA2_dtr9_all", Temp.xs)

###################################### 2B. Empirical DTR 12

# Sensitivities by holding each parameter constant
SA2_dtr12 <- SensitivityAnalysis2_hspc(predictions_bite_rate_dtr12, predictions_bc, predictions_lifespan_dtr12, predictions_gamma_dtr12, predictions_eggs_dtr12, predictions_pea, predictions_mdr)

# Get posterior quantiles for plotting
SA2_dtr12_a_out <- calcPostQuants(as.data.frame(SA2_dtr12[[1]]), "SA2_dtr12_a", Temp.xs)
SA2_dtr12_bc_out <- calcPostQuants(as.data.frame(SA2_dtr12[[2]]), "SA2_dtr12_bc", Temp.xs)
SA2_dtr12_lf_out <- calcPostQuants(as.data.frame(SA2_dtr12[[3]]), "SA2_dtr12_lf", Temp.xs)
SA2_dtr12_gamma_out <- calcPostQuants(as.data.frame(SA2_dtr12[[4]]), "SA2_dtr12_gamma", Temp.xs)
SA2_dtr12_B_out <- calcPostQuants(as.data.frame(SA2_dtr12[[5]]), "SA2_dtr12_B", Temp.xs)
SA2_dtr12_pEA_out <- calcPostQuants(as.data.frame(SA2_dtr12[[6]]), "SA2_dtr12_pEA", Temp.xs)
SA2_dtr12_MDR_out <- calcPostQuants(as.data.frame(SA2_dtr12[[7]]), "SA2_dtr12_MDR", Temp.xs)
SA2_dtr12_all_out <- calcPostQuants(as.data.frame(SA2_dtr12[[8]]), "SA2_dtr12_all", Temp.xs)

###################################### 3A. RS 9 - 3 traits

SA2_rs9_3traits <- SensitivityAnalysis2_hspc(predictions_bite_rate_rs9, predictions_bc, predictions_lifespan_rs9, predictions_gamma_rs9, predictions_eggs_rs9, predictions_pea, predictions_mdr)

# Get posterior quantiles for plotting
SA2_rs9_3traits_a_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[1]]), "SA2_rs9_3traits_a", Temp.xs)
SA2_rs9_3traits_bc_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[2]]), "SA2_rs9_3traits_bc", Temp.xs)
SA2_rs9_3traits_lf_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[3]]), "SA2_rs9_3traits_lf", Temp.xs)
SA2_rs9_3traits_gamma_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[4]]), "SA2_rs9_3traits_gamma", Temp.xs)
SA2_rs9_3traits_B_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[5]]), "SA2_rs9_3traits_B", Temp.xs)
SA2_rs9_3traits_pEA_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[6]]), "SA2_rs9_3traits_pEA", Temp.xs)
SA2_rs9_3traits_MDR_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[7]]), "SA2_rs9_3traits_MDR", Temp.xs)
SA2_rs9_3traits_all_out <- calcPostQuants(as.data.frame(SA2_rs9_3traits[[8]]), "SA2_rs9_3traits_all", Temp.xs)

###################################### 3B. RS 12 - 3 traits

# Sensitivities by holding each parameter constant
SA2_rs12_3traits <- SensitivityAnalysis2_hspc(predictions_bite_rate_rs12, predictions_bc, predictions_lifespan_rs12, predictions_gamma_rs12, predictions_eggs_rs12, predictions_pea, predictions_mdr)

# Get posterior quantiles for plotting
SA2_rs12_3traits_a_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[1]]), "SA2_rs12_3traits_a", Temp.xs)
SA2_rs12_3traits_bc_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[2]]), "SA2_rs12_3traits_bc", Temp.xs)
SA2_rs12_3traits_lf_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[3]]), "SA2_rs12_3traits_lf", Temp.xs)
SA2_rs12_3traits_gamma_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[4]]), "SA2_rs12_3traits_gamma", Temp.xs)
SA2_rs12_3traits_B_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[5]]), "SA2_rs12_3traits_B", Temp.xs)
SA2_rs12_3traits_pEA_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[6]]), "SA2_rs12_3traits_pEA", Temp.xs)
SA2_rs12_3traits_MDR_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[7]]), "SA2_rs12_3traits_MDR", Temp.xs)
SA2_rs12_3traits_all_out <- calcPostQuants(as.data.frame(SA2_rs12_3traits[[8]]), "SA2_rs12_3traits_all", Temp.xs)

###################################### 4A. RS 9 - all traits

SA2_rs9_alltraits <- SensitivityAnalysis2_hspc(predictions_bite_rate_rs9, predictions_bc_rs9, predictions_lifespan_rs9, predictions_gamma_rs9, predictions_eggs_rs9, predictions_pea_rs9, predictions_mdr_rs9)

# Get posterior quantiles for plotting
SA2_rs9_alltraits_a_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[1]]), "SA2_rs9_alltraits_a", Temp.xs)
SA2_rs9_alltraits_bc_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[2]]), "SA2_rs9_alltraits_bc", Temp.xs)
SA2_rs9_alltraits_lf_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[3]]), "SA2_rs9_alltraits_lf", Temp.xs)
SA2_rs9_alltraits_gamma_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[4]]), "SA2_rs9_alltraits_gamma", Temp.xs)
SA2_rs9_alltraits_B_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[5]]), "SA2_rs9_alltraits_B", Temp.xs)
SA2_rs9_alltraits_pEA_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[6]]), "SA2_rs9_alltraits_pEA", Temp.xs)
SA2_rs9_alltraits_MDR_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[7]]), "SA2_rs9_alltraits_MDR", Temp.xs)
SA2_rs9_alltraits_all_out <- calcPostQuants(as.data.frame(SA2_rs9_alltraits[[8]]), "SA2_rs9_alltraits_all", Temp.xs)

###################################### 4B. RS 9 - all traits

# Sensitivities by holding each parameter constant
SA2_rs12_alltraits <- SensitivityAnalysis2_hspc(predictions_bite_rate_rs12, predictions_bc_rs9, predictions_lifespan_rs12, predictions_gamma_rs12, predictions_eggs_rs12, predictions_pea_rs9, predictions_mdr_rs9)

# Get posterior quantiles for plotting
SA2_rs12_alltraits_a_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[1]]), "SA2_rs12_alltraits_a", Temp.xs)
SA2_rs12_alltraits_bc_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[2]]), "SA2_rs12_alltraits_bc", Temp.xs)
SA2_rs12_alltraits_lf_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[3]]), "SA2_rs12_alltraits_lf", Temp.xs)
SA2_rs12_alltraits_gamma_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[4]]), "SA2_rs12_alltraits_gamma", Temp.xs)
SA2_rs12_alltraits_B_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[5]]), "SA2_rs12_alltraits_B", Temp.xs)
SA2_rs12_alltraits_pEA_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[6]]), "SA2_rs12_alltraits_pEA", Temp.xs)
SA2_rs12_alltraits_MDR_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[7]]), "SA2_rs12_alltraits_MDR", Temp.xs)
SA2_rs12_alltraits_all_out <- calcPostQuants(as.data.frame(SA2_rs12_alltraits[[8]]), "SA2_rs12_alltraits_all", Temp.xs)

###################################### 5A. RS on suitability - 9

# Perform rate summation on the outputs from Constant Model 1
SA2_a_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[1]]), LPtemps_dtr9, Temp.xs)
SA2_bc_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[2]]), LPtemps_dtr9, Temp.xs)
SA2_lf_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[3]]), LPtemps_dtr9, Temp.xs)
SA2_gamma_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[4]]), LPtemps_dtr9, Temp.xs)
SA2_B_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[5]]), LPtemps_dtr9, Temp.xs)
SA2_pEA_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[6]]), LPtemps_dtr9, Temp.xs)
SA2_MDR_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[7]]), LPtemps_dtr9, Temp.xs)
SA2_all_model5_rs9 <- RSCalcTempGrad(data.frame(SA2_constant[[8]]), LPtemps_dtr9, Temp.xs)

# Get posterior quantiles for plotting
SA2_rs9_model5_a_out <- calcPostQuants(SA2_a_model5_rs9, "SA2_rs9_model5_a", Temp.xs)
SA2_rs9_model5_bc_out <- calcPostQuants(SA2_bc_model5_rs9, "SA2_rs9_model5_bc", Temp.xs)
SA2_rs9_model5_lf_out <- calcPostQuants(SA2_lf_model5_rs9, "SA2_rs9_model5_lf", Temp.xs)
SA2_rs9_model5_gamma_out <- calcPostQuants(SA2_gamma_model5_rs9, "SA2_rs9_model5_gamma", Temp.xs)
SA2_rs9_model5_B_out <- calcPostQuants(SA2_B_model5_rs9, "SA2_rs9_model5_B", Temp.xs)
SA2_rs9_model5_pEA_out <- calcPostQuants(SA2_pEA_model5_rs9, "SA2_rs9_model5_pEA", Temp.xs)
SA2_rs9_model5_MDR_out <- calcPostQuants(SA2_MDR_model5_rs9, "SA2_rs9_model5_MDR", Temp.xs)
SA2_rs9_model5_all_out <- calcPostQuants(SA2_all_model5_rs9, "SA2_rs9_model5_all", Temp.xs)

###################################### 5A. RS on suitability - 12

# Perform rate summation on the outputs from Constant Model 1
SA2_a_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[1]]), LPtemps_dtr12, Temp.xs)
SA2_bc_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[2]]), LPtemps_dtr12, Temp.xs)
SA2_lf_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[3]]), LPtemps_dtr12, Temp.xs)
SA2_gamma_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[4]]), LPtemps_dtr12, Temp.xs)
SA2_B_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[5]]), LPtemps_dtr12, Temp.xs)
SA2_pEA_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[6]]), LPtemps_dtr12, Temp.xs)
SA2_MDR_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[7]]), LPtemps_dtr12, Temp.xs)
SA2_all_model5_rs12 <- RSCalcTempGrad(data.frame(SA2_constant[[8]]), LPtemps_dtr12, Temp.xs)

# Get posterior quantiles for plotting
SA2_rs12_model5_a_out <- calcPostQuants(SA2_a_model5_rs12, "SA2_rs12_model5_a", Temp.xs)
SA2_rs12_model5_bc_out <- calcPostQuants(SA2_bc_model5_rs12, "SA2_rs12_model5_bc", Temp.xs)
SA2_rs12_model5_lf_out <- calcPostQuants(SA2_lf_model5_rs12, "SA2_rs12_model5_lf", Temp.xs)
SA2_rs12_model5_gamma_out <- calcPostQuants(SA2_gamma_model5_rs12, "SA2_rs12_model5_gamma", Temp.xs)
SA2_rs12_model5_B_out <- calcPostQuants(SA2_B_model5_rs12, "SA2_rs12_model5_B", Temp.xs)
SA2_rs12_model5_pEA_out <- calcPostQuants(SA2_pEA_model5_rs12, "SA2_rs12_model5_pEA", Temp.xs)
SA2_rs12_model5_MDR_out <- calcPostQuants(SA2_MDR_model5_rs12, "SA2_rs12_model5_MDR", Temp.xs)
SA2_rs12_model5_all_out <- calcPostQuants(SA2_all_model5_rs12, "SA2_rs12_model5_all", Temp.xs)


##########
###### 7. Uncertainty Analysis
##########

# Constant
UA_constant <- UncertaintyAnalysis(predictions_bite_rate_constant, predictions_bc, predictions_lifespan_constant, predictions_gamma_constant, predictions_eggs_constant, predictions_pea, predictions_mdr,
								   m_bite_rate_constant, m_bc, m_lifespan_constant, m_gamma_constant, m_eggs_constant, m_pea, m_mdr)

# DTR 9 
UA_dtr9 <- UncertaintyAnalysis(predictions_bite_rate_dtr9, predictions_bc, predictions_lifespan_dtr9, predictions_gamma_dtr9, predictions_eggs_dtr9, predictions_pea, predictions_mdr,
							   m_bite_rate_dtr9, m_bc, m_lifespan_dtr9, m_gamma_dtr9, m_eggs_dtr9, m_pea, m_mdr)

# DTR 12
UA_dtr12 <- UncertaintyAnalysis(predictions_bite_rate_dtr12, predictions_bc, predictions_lifespan_dtr12, predictions_gamma_dtr12, predictions_eggs_dtr12, predictions_pea, predictions_mdr,
								m_bite_rate_dtr12, m_bc, m_lifespan_dtr12, m_gamma_dtr12, m_eggs_dtr12, m_pea, m_mdr)

# RS 9 - 3 traits 
UA_rs9_3traits <- UncertaintyAnalysis(predictions_bite_rate_rs9, predictions_bc, predictions_lifespan_rs9, predictions_gamma_rs9, predictions_eggs_rs9, predictions_pea, predictions_mdr,
									  m_bite_rate_rs9, m_bc, m_lifespan_rs9, m_gamma_rs9, m_eggs_rs9, m_pea, m_mdr)

# RS 12 - 3 traits 
UA_rs12_3traits <- UncertaintyAnalysis(predictions_bite_rate_rs12, predictions_bc, predictions_lifespan_rs12, predictions_gamma_rs12, predictions_eggs_rs12, predictions_pea, predictions_mdr,
									   m_bite_rate_rs12, m_bc, m_lifespan_rs12, m_gamma_rs12, m_eggs_rs12, m_pea, m_mdr)

# RS 9 - all traits 
UA_rs9_alltraits <- UncertaintyAnalysis(predictions_bite_rate_rs9, predictions_bc_rs9, predictions_lifespan_rs9, predictions_gamma_rs9, predictions_eggs_rs9, predictions_pea_rs9, predictions_mdr_rs9,
										m_bite_rate_rs9, m_bc_rs9, m_lifespan_rs9, m_gamma_rs9, m_eggs_rs9, m_pea_rs9, m_mdr_rs9)

# RS 12 - all traits 
UA_rs12_alltraits <- UncertaintyAnalysis(predictions_bite_rate_rs12, predictions_bc_rs12, predictions_lifespan_rs12, predictions_gamma_rs12, predictions_eggs_rs12, predictions_pea_rs12, predictions_mdr_rs12,
										 m_bite_rate_rs12, m_bc_rs12, m_lifespan_rs12, m_gamma_rs12, m_eggs_rs12, m_pea_rs12, m_mdr_rs12)

#### NOTE - the UA for model 5 versions takes ~10 min each to run because it does 8 RS calculations

# S(T)-level RS - DTR 9
UA_rs9_ST <- UncertaintyAnalysis_RSModel5(predictions_bite_rate_constant, predictions_bc, predictions_lifespan_constant, predictions_gamma_constant, predictions_eggs_constant, predictions_pea, predictions_mdr,
											m_bite_rate_constant, m_bc, m_lifespan_constant, m_gamma_constant, m_eggs_constant, m_pea, m_mdr, LPtemps_dtr9, Temp.xs)

# S(T)-level RS - DTR 12
UA_rs12_ST <- UncertaintyAnalysis_RSModel5(predictions_bite_rate_constant, predictions_bc, predictions_lifespan_constant, predictions_gamma_constant, predictions_eggs_constant, predictions_pea, predictions_mdr,
											m_bite_rate_constant, m_bc, m_lifespan_constant, m_gamma_constant, m_eggs_constant, m_pea, m_mdr, LPtemps_dtr12, Temp.xs)


##########
###### 6. Plot Results
##########

###################################### Sensitivity Analysis #1 - partial derivatives

# PDF Export size: 13 x 5 inches

par(mfrow = c(1,3))

# Plot posterior median sensitivities - constant
plot(Temp.xs, dR0dT_constant$median, type = "l", lwd = 2, col = 1, xlab = expression(paste("Temperature (",degree,"C)")), xlim = c(14, 37), ylab = "Relative sensitivity", main = "Model 1: Constant T", cex.lab = 1.4)
abline(h = 0, col = "gray80")
lines(Temp.xs, dR0da_constant$median, lwd = 2, col = 2)
lines(Temp.xs, dR0dbc_constant$median, lwd = 2, col = 3)
lines(Temp.xs, dR0dlf_constant$median, lwd = 2, col = 4)
lines(Temp.xs, dR0dgamma_constant$median, lwd = 2, col = 5)
lines(Temp.xs, dR0dB_constant$median, lwd = 2, col = 6)
lines(Temp.xs, dR0dpEA_constant$median, lwd = 2, col = 7)
lines(Temp.xs, dR0dMDR_constant$median, lwd = 2, col = 8)
legend('topright', legend = c(expression(R[0]), "a", "bc", "lf", expression(gamma), "B", expression(p[EA]), "MDR"), col = c(1:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "A", cex = 1.5)

# Plot posterior median sensitivities - DTR 9
plot(Temp.xs, dR0dT_dtr9$median, type = "l", lwd = 2, col = 1, xlab = expression(paste("Temperature (",degree,"C)")), xlim = c(14, 37), ylab = "", main = "Model 2: Empirical DTR 9", cex.lab = 1.4)
abline(h = 0, col = "gray80")
lines(Temp.xs, dR0da_dtr9$median, lwd = 2, col = 2)
lines(Temp.xs, dR0dbc_dtr9$median, lwd = 2, col = 3)
lines(Temp.xs, dR0dlf_dtr9$median, lwd = 2, col = 4)
lines(Temp.xs, dR0dgamma_dtr9$median, lwd = 2, col = 5)
lines(Temp.xs, dR0dB_dtr9$median, lwd = 2, col = 6)
lines(Temp.xs, dR0dpEA_dtr9$median, lwd = 2, col = 7)
lines(Temp.xs, dR0dMDR_dtr9$median, lwd = 2, col = 8)
legend('topright', legend = c(expression(R[0]), "a", "bc", "lf", expression(gamma), "B", expression(p[EA]), "MDR"), col = c(1:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "B", cex = 1.5)

# Plot posterior median sensitivities - DTR 12
plot(Temp.xs, dR0dT_dtr12$median, type = "l", lwd = 2, col = 1, xlab = expression(paste("Temperature (",degree,"C)")), xlim = c(14, 37), ylab = "", main = "Model 2: Empirical DTR 12", cex.lab = 1.4)
abline(h = 0, col = "gray80")
lines(Temp.xs, dR0da_dtr12$median, lwd = 2, col = 2)
lines(Temp.xs, dR0dbc_dtr12$median, lwd = 2, col = 3)
lines(Temp.xs, dR0dlf_dtr12$median, lwd = 2, col = 4)
lines(Temp.xs, dR0dgamma_dtr12$median, lwd = 2, col = 5)
lines(Temp.xs, dR0dB_dtr12$median, lwd = 2, col = 6)
lines(Temp.xs, dR0dpEA_dtr12$median, lwd = 2, col = 7)
lines(Temp.xs, dR0dMDR_dtr12$median, lwd = 2, col = 8)
legend('topright', legend = c(expression(R[0]), "a", "bc", "lf", expression(gamma), "B", expression(p[EA]), "MDR"), col = c(1:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "C", cex = 1.5)


###################################### Sensitivity Analysis #2 - holding single parameters constant

###################################### Plot results

# PDF Export size: 13 x 13 inches

par(mfrow = c(3,3))

# Plot R0 curves with single traits held constant - constant
plot(SA2_constant_all_out$mean / max(SA2_constant_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(14,37),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 1: Constant T", cex.lab = 1.4)
lines(SA2_constant_a_out$mean / max(SA2_constant_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_constant_bc_out$mean / max(SA2_constant_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_constant_lf_out$mean / max(SA2_constant_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_constant_gamma_out$mean / max(SA2_constant_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_constant_B_out$mean / max(SA2_constant_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_constant_pEA_out$mean / max(SA2_constant_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_constant_MDR_out$mean / max(SA2_constant_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
# legend("topleft", legend = "A", bty = "n", adj = 1.5)
legend(x = 13.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 17.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 13.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "A", cex = 1.5)

# Plot R0 curves with single traits held constant - DTR 9
plot(SA2_dtr9_all_out$mean / max(SA2_dtr9_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(14,37),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 2: Empirical DTR 9", cex.lab = 1.4)
lines(SA2_dtr9_a_out$mean / max(SA2_dtr9_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_dtr9_bc_out$mean / max(SA2_dtr9_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_dtr9_lf_out$mean / max(SA2_dtr9_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_dtr9_gamma_out$mean / max(SA2_dtr9_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_dtr9_B_out$mean / max(SA2_dtr9_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_dtr9_pEA_out$mean / max(SA2_dtr9_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_dtr9_MDR_out$mean / max(SA2_dtr9_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 13.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 17.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 13.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "B", cex = 1.5)

# Plot R0 curves with single traits held constant - DTR 12
plot(SA2_dtr12_all_out$mean / max(SA2_dtr12_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(14,37),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 2: Empirical DTR 12", cex.lab = 1.4)
lines(SA2_dtr12_a_out$mean / max(SA2_dtr12_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_dtr12_bc_out$mean / max(SA2_dtr12_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_dtr12_lf_out$mean / max(SA2_dtr12_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_dtr12_gamma_out$mean / max(SA2_dtr12_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_dtr12_B_out$mean / max(SA2_dtr12_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_dtr12_pEA_out$mean / max(SA2_dtr12_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_dtr12_MDR_out$mean / max(SA2_dtr12_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 13.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 17.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 13.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "C", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 9 3 traits
plot(SA2_rs9_3traits_all_out$mean / max(SA2_rs9_3traits_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(14,37),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 3: Trait-level RS DTR 9 - 3 traits", cex.lab = 1.4)
lines(SA2_rs9_3traits_a_out$mean / max(SA2_rs9_3traits_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs9_3traits_bc_out$mean / max(SA2_rs9_3traits_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs9_3traits_lf_out$mean / max(SA2_rs9_3traits_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs9_3traits_gamma_out$mean / max(SA2_rs9_3traits_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs9_3traits_B_out$mean / max(SA2_rs9_3traits_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs9_3traits_pEA_out$mean / max(SA2_rs9_3traits_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs9_3traits_MDR_out$mean / max(SA2_rs9_3traits_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
legend(x = 13.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 17.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 13.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "D", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 12 3 traits
plot(SA2_rs12_3traits_all_out$mean / max(SA2_rs12_3traits_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(14,37),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 3: Trait-level RS DTR 12 - 3 traits", cex.lab = 1.4)
lines(SA2_rs12_3traits_a_out$mean / max(SA2_rs12_3traits_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs12_3traits_bc_out$mean / max(SA2_rs12_3traits_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs12_3traits_lf_out$mean / max(SA2_rs12_3traits_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs12_3traits_gamma_out$mean / max(SA2_rs12_3traits_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs12_3traits_B_out$mean / max(SA2_rs12_3traits_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs12_3traits_pEA_out$mean / max(SA2_rs12_3traits_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs12_3traits_MDR_out$mean / max(SA2_rs12_3traits_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 13.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 17.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 13.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "E", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 9 all traits
plot(SA2_rs9_alltraits_all_out$mean / max(SA2_rs9_alltraits_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(12,40),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 4: Trait-level RS DTR 9 - all traits", cex.lab = 1.4)
lines(SA2_rs9_alltraits_a_out$mean / max(SA2_rs9_alltraits_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs9_alltraits_bc_out$mean / max(SA2_rs9_alltraits_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs9_alltraits_lf_out$mean / max(SA2_rs9_alltraits_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs9_alltraits_gamma_out$mean / max(SA2_rs9_alltraits_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs9_alltraits_B_out$mean / max(SA2_rs9_alltraits_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs9_alltraits_pEA_out$mean / max(SA2_rs9_alltraits_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs9_alltraits_MDR_out$mean / max(SA2_rs9_alltraits_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 11.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 16.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 11.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "F", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 12 all traits
plot(SA2_rs12_alltraits_all_out$mean / max(SA2_rs12_alltraits_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(12,40),
	 xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 4: Trait-level RS DTR 12 - all traits", cex.lab = 1.4)
lines(SA2_rs12_alltraits_a_out$mean / max(SA2_rs12_alltraits_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs12_alltraits_bc_out$mean / max(SA2_rs12_alltraits_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs12_alltraits_lf_out$mean / max(SA2_rs12_alltraits_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs12_alltraits_gamma_out$mean / max(SA2_rs12_alltraits_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs12_alltraits_B_out$mean / max(SA2_rs12_alltraits_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs12_alltraits_pEA_out$mean / max(SA2_rs12_alltraits_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs12_alltraits_MDR_out$mean / max(SA2_rs12_alltraits_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
legend(x = 11.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 16.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 11.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "G", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 9 model 5
plot(SA2_rs9_model5_all_out$mean / max(SA2_rs9_model5_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(8,42),
	 xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 5:  S(T)-level RS DTR 9", cex.lab = 1.4)
lines(SA2_rs9_model5_a_out$mean / max(SA2_rs9_model5_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs9_model5_bc_out$mean / max(SA2_rs9_model5_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs9_model5_lf_out$mean / max(SA2_rs9_model5_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs9_model5_gamma_out$mean / max(SA2_rs9_model5_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs9_model5_B_out$mean / max(SA2_rs9_model5_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs9_model5_pEA_out$mean / max(SA2_rs9_model5_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs9_model5_MDR_out$mean / max(SA2_rs9_model5_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 7, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 13, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 7, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "H", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 12 model 5
plot(SA2_rs12_model5_all_out$mean / max(SA2_rs12_model5_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 1.1), xlim = c(8,42),
	 xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 5: S(T)-level RS DTR 12", cex.lab = 1.4)
lines(SA2_rs12_model5_a_out$mean / max(SA2_rs12_model5_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs12_model5_bc_out$mean / max(SA2_rs12_model5_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs12_model5_lf_out$mean / max(SA2_rs12_model5_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs12_model5_gamma_out$mean / max(SA2_rs12_model5_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs12_model5_B_out$mean / max(SA2_rs12_model5_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs12_model5_pEA_out$mean / max(SA2_rs12_model5_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs12_model5_MDR_out$mean / max(SA2_rs12_model5_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 7, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 13, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 7, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "I", cex = 1.5)


###################################### Uncertainty Analysis

# PDF Export size: 13 x 13 inches

par(mfrow = c(3,3))

# Plot uncertainty widths - constant
plot(Temp.xs, UA_constant$a/(UA_constant$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1), xlim = c(15, 36), cex.lab = 1.4, main = "Model 1: Constant T",
	 xlab="", ylab="Width of quantiles")
lines(Temp.xs, UA_constant$bc/(UA_constant$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_constant$lf/(UA_constant$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_constant$gamma/(UA_constant$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_constant$B/(UA_constant$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_constant$pEA/(UA_constant$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_constant$MDR/(UA_constant$R0+ec), col=8, lwd=2)
legend(18, 1.03, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "A", cex = 1.5)

# Plot uncertainty widths - DTR 9
plot(Temp.xs, UA_dtr9$a/(UA_dtr9$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1), xlim = c(15, 36), cex.lab = 1.4, main = "Model 2: Empirical DTR 9",
	 xlab="", ylab="")
lines(Temp.xs, UA_dtr9$bc/(UA_dtr9$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_dtr9$lf/(UA_dtr9$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_dtr9$gamma/(UA_dtr9$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_dtr9$B/(UA_dtr9$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_dtr9$pEA/(UA_dtr9$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_dtr9$MDR/(UA_dtr9$R0+ec), col=8, lwd=2)
legend(18, 1.03, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "B", cex = 1.5)

# Plot uncertainty widths - DTR 12
plot(Temp.xs, UA_dtr12$a/(UA_dtr12$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1), xlim = c(15, 36), cex.lab = 1.4, main = "Model 2: Empirical DTR 12",
	 xlab="", ylab="")
lines(Temp.xs, UA_dtr12$bc/(UA_dtr12$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_dtr12$lf/(UA_dtr12$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_dtr12$gamma/(UA_dtr12$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_dtr12$B/(UA_dtr12$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_dtr12$pEA/(UA_dtr12$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_dtr12$MDR/(UA_dtr12$R0+ec), col=8, lwd=2)
legend(18, 1.03, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "C", cex = 1.5)

# Plot uncertainty widths - RS 9 3 traits
plot(Temp.xs, UA_rs9_3traits$a/(UA_rs9_3traits$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1.05), xlim = c(15, 36), cex.lab = 1.4, main = "Model 3: Trait-level RS DTR 9 - 3 traits",
	 xlab="", ylab="Width of quantiles")
lines(Temp.xs, UA_rs9_3traits$bc/(UA_rs9_3traits$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_rs9_3traits$lf/(UA_rs9_3traits$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_rs9_3traits$gamma/(UA_rs9_3traits$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_rs9_3traits$B/(UA_rs9_3traits$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_rs9_3traits$pEA/(UA_rs9_3traits$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_rs9_3traits$MDR/(UA_rs9_3traits$R0+ec), col=8, lwd=2)
legend(18, 1.08, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "D", cex = 1.5)

# Plot uncertainty widths - RS 12 3 traits
plot(Temp.xs, UA_rs12_3traits$a/(UA_rs12_3traits$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1.05), xlim = c(15, 36), cex.lab = 1.4, main = "Model 3: Trait-level RS DTR 12 - 3 traits",
	 xlab="", ylab="")
lines(Temp.xs, UA_rs12_3traits$bc/(UA_rs12_3traits$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_rs12_3traits$lf/(UA_rs12_3traits$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_rs12_3traits$gamma/(UA_rs12_3traits$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_rs12_3traits$B/(UA_rs12_3traits$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_rs12_3traits$pEA/(UA_rs12_3traits$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_rs12_3traits$MDR/(UA_rs12_3traits$R0+ec), col=8, lwd=2)
legend(18, 1.08, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "E", cex = 1.5)

# Plot uncertainty widths - RS 9 all traits
plot(Temp.xs, UA_rs9_alltraits$a/(UA_rs9_alltraits$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1.05), xlim = c(10, 40), cex.lab = 1.4, main = "Model 4: Trait-level RS DTR 9 - all traits",
	 xlab="", ylab="")
lines(Temp.xs, UA_rs9_alltraits$bc/(UA_rs9_alltraits$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_rs9_alltraits$lf/(UA_rs9_alltraits$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_rs9_alltraits$gamma/(UA_rs9_alltraits$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_rs9_alltraits$B/(UA_rs9_alltraits$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_rs9_alltraits$pEA/(UA_rs9_alltraits$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_rs9_alltraits$MDR/(UA_rs9_alltraits$R0+ec), col=8, lwd=2)
legend(32, 1.08, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "F", cex = 1.5)

# Plot uncertainty widths - RS 12 all traits
plot(Temp.xs, UA_rs12_alltraits$a/(UA_rs12_alltraits$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1.05), xlim = c(8, 41), cex.lab = 1.4, main = "Model 4: Trait-level RS DTR 12 - all traits",
	 xlab=expression(paste("Temperature (",degree,"C)")), ylab="Width of quantiles")
lines(Temp.xs, UA_rs12_alltraits$bc/(UA_rs12_alltraits$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_rs12_alltraits$lf/(UA_rs12_alltraits$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_rs12_alltraits$gamma/(UA_rs12_alltraits$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_rs12_alltraits$B/(UA_rs12_alltraits$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_rs12_alltraits$pEA/(UA_rs12_alltraits$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_rs12_alltraits$MDR/(UA_rs12_alltraits$R0+ec), col=8, lwd=2)
legend(32, 1.08, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "G", cex = 1.5)

# Plot uncertainty widths - RS 9 S(T)
plot(Temp.xs, UA_rs9_ST$a/(UA_rs9_ST$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1.05), xlim = c(8, 41), cex.lab = 1.4, main = "Model 5:  S(T)-level RS DTR 9",
	 xlab=expression(paste("Temperature (",degree,"C)")), ylab="")
lines(Temp.xs, UA_rs9_ST$bc/(UA_rs9_ST$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_rs9_ST$lf/(UA_rs9_ST$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_rs9_ST$gamma/(UA_rs9_ST$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_rs9_ST$B/(UA_rs9_ST$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_rs9_ST$pEA/(UA_rs9_ST$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_rs9_ST$MDR/(UA_rs9_ST$R0+ec), col=8, lwd=2)
legend(32, 1.08, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "H", cex = 1.5)

# Plot uncertainty widths - RS 12 S(T)
plot(Temp.xs, UA_rs12_ST$a/(UA_rs12_ST$R0+ec), col=2, type="l", lwd=2, ylim = c(0,1.05), xlim = c(8, 41), cex.lab = 1.4, main = "Model 5:  S(T)-level RS DTR 12",
	 xlab=expression(paste("Temperature (",degree,"C)")), ylab="")
lines(Temp.xs, UA_rs12_ST$bc/(UA_rs12_ST$R0+ec), col=3, lwd=2)
lines(Temp.xs, UA_rs12_ST$lf/(UA_rs12_ST$R0+ec), col=4, lwd=2)
lines(Temp.xs, UA_rs12_ST$gamma/(UA_rs12_ST$R0+ec), col=5, lwd=2)
lines(Temp.xs, UA_rs12_ST$B/(UA_rs12_ST$R0+ec), col=6, lwd=2)
lines(Temp.xs, UA_rs12_ST$pEA/(UA_rs12_ST$R0+ec), col=7, lwd=2)
lines(Temp.xs, UA_rs12_ST$MDR/(UA_rs12_ST$R0+ec), col=8, lwd=2)
legend(32, 1.08, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
legend("topleft", bty = "n", legend = "I", cex = 1.5)


###################################### Zooming in Tmax and Tmin - Sensitivity Analysis 2 for models 4 and 5

par(mfrow = c(2,2))

# Plot R0 curves with single traits held constant - RS 9 all traits
plot(SA2_rs9_alltraits_all_out$mean / max(SA2_rs9_alltraits_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 0.015), xlim = c(39,41),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 4: Trait-level RS DTR 9 - all traits", cex.lab = 1.4)
lines(SA2_rs9_alltraits_a_out$mean / max(SA2_rs9_alltraits_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs9_alltraits_bc_out$mean / max(SA2_rs9_alltraits_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs9_alltraits_lf_out$mean / max(SA2_rs9_alltraits_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs9_alltraits_gamma_out$mean / max(SA2_rs9_alltraits_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs9_alltraits_B_out$mean / max(SA2_rs9_alltraits_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs9_alltraits_pEA_out$mean / max(SA2_rs9_alltraits_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs9_alltraits_MDR_out$mean / max(SA2_rs9_alltraits_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 11.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
text(x = 16.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 11.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "A", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 12 all traits
plot(SA2_rs12_alltraits_all_out$mean / max(SA2_rs12_alltraits_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 0.01), xlim = c(39.5,41.5),
	 xlab = "", cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 4: Trait-level RS DTR 12 - all traits", cex.lab = 1.4)
lines(SA2_rs12_alltraits_a_out$mean / max(SA2_rs12_alltraits_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs12_alltraits_bc_out$mean / max(SA2_rs12_alltraits_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs12_alltraits_lf_out$mean / max(SA2_rs12_alltraits_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs12_alltraits_gamma_out$mean / max(SA2_rs12_alltraits_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs12_alltraits_B_out$mean / max(SA2_rs12_alltraits_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs12_alltraits_pEA_out$mean / max(SA2_rs12_alltraits_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs12_alltraits_MDR_out$mean / max(SA2_rs12_alltraits_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 11.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 16.5, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 11.5, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "B", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 9 model 5
plot(SA2_rs9_model5_all_out$mean / max(SA2_rs9_model5_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 0.025), xlim = c(39.5,41.5),
	 xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 5:  S(T)-level RS DTR 9", cex.lab = 1.4)
lines(SA2_rs9_model5_a_out$mean / max(SA2_rs9_model5_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs9_model5_bc_out$mean / max(SA2_rs9_model5_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs9_model5_lf_out$mean / max(SA2_rs9_model5_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs9_model5_gamma_out$mean / max(SA2_rs9_model5_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs9_model5_B_out$mean / max(SA2_rs9_model5_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs9_model5_pEA_out$mean / max(SA2_rs9_model5_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs9_model5_MDR_out$mean / max(SA2_rs9_model5_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 7, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
text(x = 13, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 7, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "C", cex = 1.5)

# Plot R0 curves with single traits held constant - RS 12 model 5
plot(SA2_rs12_model5_all_out$mean / max(SA2_rs12_model5_all_out$mean) ~ Temp.xs, type = "l", ylim = c(0, 0.075), xlim = c(40.5,42.5),
	 xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2, main = "Model 5: S(T)-level RS DTR 12", cex.lab = 1.4)
lines(SA2_rs12_model5_a_out$mean / max(SA2_rs12_model5_a_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(SA2_rs12_model5_bc_out$mean / max(SA2_rs12_model5_bc_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(SA2_rs12_model5_lf_out$mean / max(SA2_rs12_model5_lf_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(SA2_rs12_model5_gamma_out$mean / max(SA2_rs12_model5_gamma_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(SA2_rs12_model5_B_out$mean / max(SA2_rs12_model5_B_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(SA2_rs12_model5_pEA_out$mean / max(SA2_rs12_model5_pEA_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(SA2_rs12_model5_MDR_out$mean / max(SA2_rs12_model5_MDR_out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
legend(x = 7, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 1)
text(x = 13, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 7, y = 0.85, col = c(2:8), lwd = 1.5,
	   lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
legend("topleft", bty = "n", legend = "D", cex = 1.5)