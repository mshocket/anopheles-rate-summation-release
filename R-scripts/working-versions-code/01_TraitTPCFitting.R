## Rate Summation Project
## Script for fitting trait thermal performance curves (TPCs) using JAGS
## Written by Marta Shocket & Joey Bernhardt in 2022-2024
## Modified from code provided by Kerri Miazgowicz in 2020, which was modified from code provided by Marta Shocket in 2018

##	Table of contents:
##
##	1. Set-up workspace
##  2. Shared settings for all models
##	3. Fitting trait TPCs - individual data from fluctuating temperature experiment: biting rate (a), lifespan (lf), & lifetime eggs (B)
##	4. Fitting trait TPCs - group mean data from previous experiments: vector competence (bc), extrinsic incubation period (EIP), 
##		juvenile survival (pEA), and mosquito development rate (MDR)
##	5. Fitting trait TPCs - gamma - block mean data - calculated by combining fluctuating temperature experiment survival data with
##		EIP data from previous experiments
##	6. Process and save TPC model output for plotting

# NOTE: mcmcplot() commands are currently commented out - inspecting the MCMC chains is important for diagnostics when initially 
# fitting models. However, it also takes time and opens up a browser window, and it not necessary to just reproduce this analysis.


##########
###### 1. Set-up workspace 
##########

##### Load JAGS libraries
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits

##### Load functions
source("R-scripts/working-versions-code/00_RSProjectFunctions.R")

##### Create JAGS models
# NOTES:
# - Running the code below writes a .txt file for each model to your current working directory.
#   	It must be saved there in order for JAGS to find it.
# - The source() command only needs to be run **once** to create the .txt files, which can be reused
#		in later fitting sessions
# - The models include a section for 'derived quantities and predictions' that calculates the 
#   	trait across a temperature gradient for each saved sample in the MCMC chain. This output 
#   	is what we use later on to calculate R0. 

# Create .txt files of JAGS models - un-comment the line below the first time you use this script
# source("R-scripts/working-versions-code/00_JAGSModels.R")

##### Load raw trait data from fluctuating temperature experiment and change column from "Treatment" to "temp"
data.constant <- read.csv("data/constant.individual.trait.csv")
names(data.constant)[names(data.constant) == 'Treatment'] <- 'temp'
data.fluc9 <- read.csv("data/fluc9.individual.trait.csv")
names(data.fluc9)[names(data.fluc9) == 'Treatment'] <- 'temp'
data.fluc12 <- read.csv("data/fluc12.individual.trait.csv")
names(data.fluc12)[names(data.fluc12) == 'Treatment'] <- 'temp'

##### Load data from previous studies for other traits
data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")
data.bc.EIP <- read.csv("data/ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50

##### Load gamma calculations
data.gamma.constant <- read.csv("data/gamma_values_constant.csv")
data.gamma.dtr9 <- read.csv("data/gamma_values_dtr9.csv")
data.gamma.dtr12 <- read.csv("data/gamma_values_dtr12.csv")

##### Counting sample sizes
# Constant
data.constant.ss <- data.constant |>
	group_by(temp) |>
	summarise(ss = length(bite.rate)) |>
	ungroup()

# DTR 9 
data.fluc9.ss <- data.fluc9 |>
	group_by(temp) |>
	summarise(ss = length(bite.rate)) |>
	ungroup()

# DTR 12
data.fluc12.ss <- data.fluc12 |>
	group_by(temp) |>
	summarise(ss = length(bite.rate)) |>
	ungroup()

##### Make posterior output folder
dir.create("saved-posteriors", showWarnings = FALSE)
dir.create("data-processed", showWarnings = FALSE)

##########
###### 2. Shared settings for all models
##########

##### MCMC Settings
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Derived Quantity Settings
Temp.xs <- seq(0, 45, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

# Save the temperature sequence for future analyses
save(Temp.xs, file = "saved-posteriors/temps.Rdata")



##########
###### 3. Fitting traits - individual data from fluctuating temperature experiment: biting rate (a), lifespan (lf), & lifetime eggs (B)
##########


###################################################################################### Model fitting: bite rate (a)


###################################### Settings for fitting normal / truncated normal distribution

#####  inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Bite rate at constant temperature ---------------------------------------

##### Get data
data_bite_rate_constant <- with(data.constant, data.frame('T' = temp, 'trait' = bite.rate))

##### Organize Data for JAGS
trait <- data_bite_rate_constant$trait
N.obs <- length(trait)
temp <- data_bite_rate_constant$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_bite_rate_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_bite_rate_constantq <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
								 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_constant$BUGSoutput$summary[1:5,]
#mcmcplot(model_bite_rate_constant)

# DIC for Briere is lower
model_bite_rate_constant$BUGSoutput$DIC
model_bite_rate_constantq$BUGSoutput$DIC

# ##### Save model output 
save(model_bite_rate_constant, file = "saved-posteriors/constant_biterate.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_constant, 
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Bite rate at DTR 9 ---------------------------------------

##### Get data
data_bite_rate_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = bite.rate))

##### Organize Data for JAGS
trait <- data_bite_rate_dtr9$trait
N.obs <- length(trait)
temp <- data_bite_rate_dtr9$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_bite_rate_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_bite_rate_dtr9q <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
								 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_dtr9$BUGSoutput$summary[1:5,]
#mcmcplot(model_bite_rate_dtr9)

# DIC for Briere is lower
model_bite_rate_dtr9$BUGSoutput$DIC
model_bite_rate_dtr9q$BUGSoutput$DIC

# ##### Save model output 
save(model_bite_rate_dtr9, file = "saved-posteriors/dtr9_biterate.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_dtr9,
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Bite rate at DTR 12 ---------------------------------------

##### Get data
data_bite_rate_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = bite.rate))

##### Organize Data for JAGS
trait <- data_bite_rate_dtr12$trait
N.obs <- length(trait)
temp <- data_bite_rate_dtr12$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_bite_rate_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_bite_rate_dtr12q <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
							  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_dtr12$BUGSoutput$summary[1:5,]
#mcmcplot(model_bite_rate_dtr12)

# DIC for Briere is lower
model_bite_rate_dtr12$BUGSoutput$DIC
model_bite_rate_dtr12q$BUGSoutput$DIC

# ##### Save model output
save(model_bite_rate_dtr12, file = "saved-posteriors/dtr12_biterate.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_dtr12,
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Bite rate all data (for model comparison via DIC) ---------------------------------------

##### Get data
data_bite_rate_all <- rbind(data_bite_rate_constant, data_bite_rate_dtr9, data_bite_rate_dtr12)
#data_bite_rate_all <- rbind(data_bite_rate_dtr9, data_bite_rate_dtr12)

##### Organize Data for JAGS
trait <- data_bite_rate_all$trait
N.obs <- length(trait)
temp <- data_bite_rate_all$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_bite_rate_all <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
								 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_all$BUGSoutput$summary[1:5,]
#mcmcplot(model_bite_rate_all)

# ##### Save model output 
save(model_bite_rate_all, file = "saved-posteriors/all_biterate.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_all, 
	 ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################################################################### Model fitting: lifespan (lf)


###################################### Settings for fitting gamma distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.ra = 0.003)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.ra", "z.trait.mu.pred")


##### Lifespan at constant temperature ---------------------------------------

##### Get data
data_lifespan_constant <- with(data.constant, data.frame('T' = temp, 'trait' = lifespan))

##### Organize Data for JAGS
trait <- data_lifespan_constant$trait
N.obs <- length(trait)
temp <- data_lifespan_constant$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_lifespan_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_gamma.txt",
				  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_lifespan_constantb <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_gamma.txt",
								n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_lifespan_constant$BUGSoutput$summary[1:5,]
#mcmcplot(model_lifespan_constant)

# DIC for quadratic is lower
model_lifespan_constant$BUGSoutput$DIC
model_lifespan_constantb$BUGSoutput$DIC

##### Save model output
save(model_lifespan_constant, file = "saved-posteriors/constant_lifespan.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_lifespan_constant,
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_lifespan_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifespan at DTR 9 ---------------------------------------

##### Get data
data_lifespan_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifespan))

##### Organize Data for JAGS
trait <- data_lifespan_dtr9$trait
N.obs <- length(trait)
temp <- data_lifespan_dtr9$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_lifespan_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_gamma.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_lifespan_dtr9b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_gamma.txt",
							n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_lifespan_dtr9$BUGSoutput$summary[1:5,]
#mcmcplot(model_lifespan_dtr9)

# DIC for quadratic is lower
model_lifespan_dtr9$BUGSoutput$DIC
model_lifespan_dtr9b$BUGSoutput$DIC

##### Save model output
save(model_lifespan_dtr9, file = "saved-posteriors/dtr9_lifespan.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_lifespan_dtr9,
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_lifespan_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifespan at DTR 12 ---------------------------------------

##### Get data
data_lifespan_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifespan))

##### Organize Data for JAGS
trait <- data_lifespan_dtr12$trait
N.obs <- length(trait)
temp <- data_lifespan_dtr12$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_lifespan_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_gamma.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_lifespan_dtr12b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_gamma.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_lifespan_dtr12$BUGSoutput$summary[1:5,]
#mcmcplot(model_lifespan_dtr12)

# DIC for quadratic is lower
model_lifespan_dtr12$BUGSoutput$DIC
model_lifespan_dtr12b$BUGSoutput$DIC

##### Save model output
save(model_lifespan_dtr12, file = "saved-posteriors/dtr12_lifespan.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_lifespan_dtr12,
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifespan all data (for model comparison via DIC) ---------------------------------------

##### Get data
data_lifespan_all <- rbind(data_lifespan_constant, data_lifespan_dtr9, data_lifespan_dtr12)
#data_lifespan_all <- rbind(data_lifespan_dtr9, data_lifespan_dtr12)

##### Organize Data for JAGS
trait <- data_lifespan_all$trait
N.obs <- length(trait)
temp <- data_lifespan_all$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_lifespan_all <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_gamma.txt",
								n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_lifespan_all$BUGSoutput$summary[1:5,]
#mcmcplot(model_lifespan_all)

##### Save model output
save(model_lifespan_all, file = "saved-posteriors/all_lifespan.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_lifespan_all,
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_lifespan_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################################################################### Model fitting: lifetime eggs (B)


###################################### Settings for fitting negative binomial distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.r = 1)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.r", "z.trait.mu.pred")

##### Lifetime eggs at constant -----------------------------------

##### Get data
data_eggs_constant <- with(data.constant, data.frame('T' = temp, 'trait' = lifetime.eggs))

##### Organize Data for JAGS
trait <- data_eggs_constant$trait
N.obs <- length(trait)
temp <- data_eggs_constant$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_eggs_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_negbin.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_eggs_constant_b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_negbin.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())


##### Examine model output & run diagnostics
model_eggs_constant$BUGSoutput$summary[1:5,]
#mcmcplot(model_eggs_constant)

# DIC for quadratic is lower
model_eggs_constant$BUGSoutput$DIC
model_eggs_constant_b$BUGSoutput$DIC

##### Save model output
save(model_eggs_constant, file = "saved-posteriors/constant_eggs.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_constant,
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_constant,
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_constant_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_constant_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_constant_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifetime eggs at DTR 9 ------------------------------------

##### Get data
data_eggs_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_eggs_dtr9$trait
N.obs <- length(trait)
temp <- data_eggs_dtr9$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_eggs_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_negbin.txt",
						 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_eggs_dtr9_b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_negbin.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_eggs_dtr9$BUGSoutput$summary[1:5,]
#mcmcplot(model_eggs_dtr9)

# DIC for Briere is slightly lower
model_eggs_dtr9$BUGSoutput$DIC
model_eggs_dtr9_b$BUGSoutput$DIC

##### Save model output
save(model_eggs_dtr9, file = "saved-posteriors/dtr9_eggs.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr9,
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr9,
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_dtr9_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr9_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr9_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifetime eggs at DTR 12 ------------------------------------

##### Get data
data_eggs_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifetime.eggs)) 

##### Organize Data for JAGS
trait <- data_eggs_dtr12$trait
N.obs <- length(trait)
temp <- data_eggs_dtr12$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_eggs_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_negbin.txt",
						  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_eggs_dtr12_b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_negbin.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())


##### Examine model output & run diagnostics
model_eggs_dtr12$BUGSoutput$summary[1:5,]
#mcmcplot(model_eggs_dtr12)

# DIC for quadratic is lower
model_eggs_dtr12$BUGSoutput$DIC
model_eggs_dtr12_b$BUGSoutput$DIC

##### Save model output
save(model_eggs_dtr12, file = "saved-posteriors/dtr12_eggs.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr12,
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr12,
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_dtr12_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr12_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr12_b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifetime eggs all data (for model comparison via DIC)  -----------------------------------

##### Get data
data_eggs_all <- rbind(data_eggs_constant, data_eggs_dtr9, data_eggs_dtr12)
#data_eggs_all <- rbind(data_eggs_dtr9, data_eggs_dtr12)

##### Organize Data for JAGS
trait <- data_eggs_all$trait
N.obs <- length(trait)
temp <- data_eggs_all$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_eggs_all <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_negbin.txt",
							n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Run JAGS - **select correct model file**
model_eggs_all_b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_negbin.txt",
							  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_eggs_all$BUGSoutput$summary[1:5,]
#mcmcplot(model_eggs_constant)

# DIC for quadratic is lower
model_eggs_all$BUGSoutput$DIC
model_eggs_all_b$BUGSoutput$DIC

##### Save model output
save(model_eggs_all, file = "saved-posteriors/all_eggs.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_all,
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##########
###### 4. Fitting traits - group mean data from previous experiments: vector competence (bc), extrinic incubation period (EIP), juvenile survival (pEA), and mosquito development rate (MDR)
##########

#### The fluctating temperature experiment measured bite rate (a), adult mosquito lifespan (lf), lifetime egg production (B). 
# To calculate S(T) we need addtional parameters: vector competence (bc; the proportion of infectious mosquitoes),
#   parasite development rate or extrinsic incubation period (PDR = 1/EIP), probability of egg to adult survival (pEA), 
#   and mosquito development rate (MDR)


###################################### Settings for fitting normal / truncated normal distribution

#####  inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")



####################################### Model fitting: EIP50 ----------------------------------------------------

############ EIP50 from Shapiro et al. 2017 Plos Biology
# Get data
data_eip50 <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = 1/EIP50)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_eip50$trait
N.obs <- length(trait)
temp <- data_eip50$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_eip50 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_eip50$BUGSoutput$summary[1:5,]
#mcmcplot(model_eip50)

# Save model output
save(model_eip50, file = "saved-posteriors/EIP50.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,0.2), data = data_eip50,
     ylab = "EIP-50", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


####################################### Model fitting: Vector competence ----------------------------------------

############ Vector competence from Shapiro et al. 2017 Plos Biology
# Get data
data_bc <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_bc$trait
N.obs <- length(trait)
temp <- data_bc$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_bc <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_bc$BUGSoutput$summary[1:5,]
#mcmcplot(model_bc)

# Save model output
save(model_bc, file = "saved-posteriors/bc.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_bc,
     ylab = "Vector competence", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


####################################### Model fitting: pEA ------------------------------------------------------

############ pEA - Quadratic from Paaijmans 2013 Global Climate Change
# Get data
data_pea <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = Pea)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_pea$trait
N.obs <- length(trait)
temp <- data_pea$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_pea <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_pea$BUGSoutput$summary[1:5,]
#mcmcplot(model_pea)

# Save model output
save(model_pea, file = "saved-posteriors/pea.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_pea,
     ylab = "Pea", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


####################################### Model fitting: MDR ------------------------------------------------------

############ Mosqutio development rate from Paaijmans 2013 Global Climate Change
# Get data
data_mdr <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = MDR)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_mdr$trait
N.obs <- length(trait)
temp <- data_mdr$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_mdr <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_mdr$BUGSoutput$summary[1:5,]
#mcmcplot(model_mdr)

# Save model output
save(model_mdr, file = "saved-posteriors/MDR.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_mdr,
     ylab = "MDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 5. Fitting traits - gamma - block mean data - calculated by combining fluctuating temperature experiment survival data with EIP data from previous experiments
##########


#### Gamma - the proportion of mosquitoes surviving the EIP50 - is calculated based on Kaplan-Meier survival estimates
#  from our experiment (with constant and fluctuating temperatures) and the EIP50 values estimated for our
#  temperature treatments based on the TPC fit to data from previous studies (above, constant temperature only).


###################################################################################### Model fitting: Gamma


###################################### Settings for fitting normal / truncated normal distribution

#####  inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Gamma at constant -----------------------------------

# Get data
data_gamma_constant <- with(data.gamma.constant, data.frame('T' = temp, 'trait' = gamma)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_gamma_constant$trait
N.obs <- length(trait)
temp <- data_gamma_constant$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_gamma_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
					n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_gamma_constant$BUGSoutput$summary[1:5,]
#mcmcplot(model_gamma_constant)

# Save model output
save(model_gamma_constant, file = "saved-posteriors/constant_gamma.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_gamma_constant,
	 ylab = "Gamma", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_gamma_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Gamma at dtr9 -----------------------------------

# Get data
data_gamma_dtr9 <- with(data.gamma.dtr9, data.frame('T' = temp, 'trait' = gamma)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_gamma_dtr9$trait
N.obs <- length(trait)
temp <- data_gamma_dtr9$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_gamma_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_gamma_dtr9$BUGSoutput$summary[1:5,]
#mcmcplot(model_gamma_dtr9)

# Save model output
save(model_gamma_dtr9, file = "saved-posteriors/dtr9_gamma.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_gamma_dtr9,
	 ylab = "Gamma", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_gamma_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Gamma at dtr12 -----------------------------------

# Get data
data_gamma_dtr12 <- with(data.gamma.dtr12, data.frame('T' = temp, 'trait' = gamma)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_gamma_dtr12$trait
N.obs <- length(trait)
temp <- data_gamma_dtr12$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_gamma_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_gamma_dtr12$BUGSoutput$summary[1:5,]
#mcmcplot(model_gamma_dtr12)

# Save model output
save(model_gamma_dtr12, file = "saved-posteriors/dtr12_gamma.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_gamma_dtr12,
	 ylab = "Gamma", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_gamma_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Gamma all data (for model comparison via DIC) -----------------------------------

# Get data
data_gamma_all <- rbind(data_gamma_constant, data_gamma_dtr9, data_gamma_dtr12)
data_gamma_all <- rbind(data_gamma_dtr9, data_gamma_dtr12)

# Organize Data for JAGS
trait <- data_gamma_all$trait
N.obs <- length(trait)
temp <- data_gamma_all$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_gamma_all <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
						  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_gamma_all$BUGSoutput$summary[1:5,]
#mcmcplot(model_gamma_all)

# Save model output
save(model_gamma_all, file = "saved-posteriors/all_gamma.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_gamma_all,
	 ylab = "Gamma", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_gamma_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_all$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##########
###### 6. Process and save model output for plotting
##########

############## 1A. Bite rate: constant ---------------------------------------------------------------

# Analyze TPC model
bite_rate_constant_TPC_analysis <- extractTPC(model_bite_rate_constant, "bite_rate_constant", Temp.xs)
predictions_bite_rate_constant_summary <- bite_rate_constant_TPC_analysis[[1]]
params_bite_rate_constant_summary <- bite_rate_constant_TPC_analysis[[2]]
params_bite_rate_constant_fullposts <- bite_rate_constant_TPC_analysis[[3]]

write_csv(predictions_bite_rate_constant_summary, "data-processed/predictions_bite_rate_constant_summary.csv")
write_csv(params_bite_rate_constant_summary, "data-processed/params_bite_rate_constant_summary.csv")
write_csv(params_bite_rate_constant_fullposts, "data-processed/params_bite_rate_constant_fullposts.csv")

# Process trait data for plotting
data_bite_rate_constant_summary <- processTraitData(data_bite_rate_constant, "bite rate", "bite_rate_constant")

write_csv(data_bite_rate_constant_summary, "data-processed/data_bite_rate_constant_summary.csv")


############## 1B. Bite rate: DTR 9 ---------------------------------------------------------------

# Analyze TPC model
bite_rate_dtr9_TPC_analysis <- extractTPC(model_bite_rate_dtr9, "bite_rate_dtr9", Temp.xs)
predictions_bite_rate_dtr9_summary <- bite_rate_dtr9_TPC_analysis[[1]]
params_bite_rate_dtr9_summary <- bite_rate_dtr9_TPC_analysis[[2]]
params_bite_rate_dtr9_fullposts <- bite_rate_dtr9_TPC_analysis[[3]]

write_csv(predictions_bite_rate_dtr9_summary, "data-processed/predictions_bite_rate_dtr9_summary.csv")
write_csv(params_bite_rate_dtr9_summary, "data-processed/params_bite_rate_dtr9_summary.csv")
write_csv(params_bite_rate_dtr9_fullposts, "data-processed/params_bite_rate_dtr9_fullposts.csv")

# Process trait data for plotting
data_bite_rate_dtr9_summary <- processTraitData(data_bite_rate_dtr9, "bite rate", "bite_rate_dtr9")

write_csv(data_bite_rate_dtr9_summary, "data-processed/data_bite_rate_dtr9_summary.csv")


############## 1C. Bite rate: DTR 12 ---------------------------------------------------------------

# Analyze TPC model
bite_rate_dtr12_TPC_analysis <- extractTPC(model_bite_rate_dtr12, "bite_rate_dtr12", Temp.xs)
predictions_bite_rate_dtr12_summary <- bite_rate_dtr12_TPC_analysis[[1]]
params_bite_rate_dtr12_summary <- bite_rate_dtr12_TPC_analysis[[2]]
params_bite_rate_dtr12_fullposts <- bite_rate_dtr12_TPC_analysis[[3]]

write_csv(predictions_bite_rate_dtr12_summary, "data-processed/predictions_bite_rate_dtr12_summary.csv")
write_csv(params_bite_rate_dtr12_summary, "data-processed/params_bite_rate_dtr12_summary.csv")
write_csv(params_bite_rate_dtr12_fullposts, "data-processed/params_bite_rate_dtr12_fullposts.csv")

# Process trait data for plotting
data_bite_rate_dtr12_summary <- processTraitData(data_bite_rate_dtr12, "bite rate", "bite_rate_dtr12")

write_csv(data_bite_rate_dtr12_summary, "data-processed/data_bite_rate_dtr12_summary.csv")


############## 2A. Lifespan: constant ---------------------------------------------------------------

# Analyze TPC model
lifespan_constant_TPC_analysis <- extractTPC(model_lifespan_constant, "lifespan_constant", Temp.xs)
predictions_lifespan_constant_summary <- lifespan_constant_TPC_analysis[[1]]
params_lifespan_constant_summary <- lifespan_constant_TPC_analysis[[2]]
params_lifespan_constant_fullposts <- lifespan_constant_TPC_analysis[[3]]

write_csv(predictions_lifespan_constant_summary, "data-processed/predictions_lifespan_constant_summary.csv")
write_csv(params_lifespan_constant_summary, "data-processed/params_lifespan_constant_summary.csv")
write_csv(params_lifespan_constant_fullposts, "data-processed/params_lifespan_constant_fullposts.csv")

# Process trait data for plotting
data_lifespan_constant_summary <- processTraitData(data_lifespan_constant, "lifespan", "lifespan_constant")

write_csv(data_lifespan_constant_summary, "data-processed/data_lifespan_constant_summary.csv")


############## 2B. Lifespan: DTR 9 ---------------------------------------------------------------

# Analyze TPC model
lifespan_dtr9_TPC_analysis <- extractTPC(model_lifespan_dtr9, "lifespan_dtr9", Temp.xs)
predictions_lifespan_dtr9_summary <- lifespan_dtr9_TPC_analysis[[1]]
params_lifespan_dtr9_summary <- lifespan_dtr9_TPC_analysis[[2]]
params_lifespan_dtr9_fullposts <- lifespan_dtr9_TPC_analysis[[3]]

write_csv(predictions_lifespan_dtr9_summary, "data-processed/predictions_lifespan_dtr9_summary.csv")
write_csv(params_lifespan_dtr9_summary, "data-processed/params_lifespan_dtr9_summary.csv")
write_csv(params_lifespan_dtr9_fullposts, "data-processed/params_lifespan_dtr9_fullposts.csv")

# Process trait data for plotting
data_lifespan_dtr9_summary <- processTraitData(data_lifespan_dtr9, "lifespan", "lifespan_dtr9")

write_csv(data_lifespan_dtr9_summary, "data-processed/data_lifespan_dtr9_summary.csv")


############## 2C. Lifespan: DTR 12 ---------------------------------------------------------------

# Analyze TPC model
lifespan_dtr12_TPC_analysis <- extractTPC(model_lifespan_dtr12, "lifespan_dtr12", Temp.xs)
predictions_lifespan_dtr12_summary <- lifespan_dtr12_TPC_analysis[[1]]
params_lifespan_dtr12_summary <- lifespan_dtr12_TPC_analysis[[2]]
params_lifespan_dtr12_fullposts <- lifespan_dtr12_TPC_analysis[[3]]

write_csv(predictions_lifespan_dtr12_summary, "data-processed/predictions_lifespan_dtr12_summary.csv")
write_csv(params_lifespan_dtr12_summary, "data-processed/params_lifespan_dtr12_summary.csv")
write_csv(params_lifespan_dtr12_fullposts, "data-processed/params_lifespan_dtr12_fullposts.csv")

# Process trait data for plotting
data_lifespan_dtr12_summary <- processTraitData(data_lifespan_dtr12, "lifespan", "lifespan_dtr12")

write_csv(data_lifespan_dtr12_summary, "data-processed/data_lifespan_dtr12_summary.csv")


############## 3A. Lifetime Eggs: constant ---------------------------------------------------------------

# Analyze TPC model
eggs_constant_TPC_analysis <- extractTPC(model_eggs_constant, "eggs_constant", Temp.xs)
predictions_eggs_constant_summary <- eggs_constant_TPC_analysis[[1]]
params_eggs_constant_summary <- eggs_constant_TPC_analysis[[2]]
params_eggs_constant_fullposts <- eggs_constant_TPC_analysis[[3]]

write_csv(predictions_eggs_constant_summary, "data-processed/predictions_eggs_constant_summary.csv")
write_csv(params_eggs_constant_summary, "data-processed/params_eggs_constant_summary.csv")
write_csv(params_eggs_constant_fullposts, "data-processed/params_eggs_constant_fullposts.csv")

# Process trait data for plotting
data_eggs_constant_summary <- processTraitData(data_eggs_constant, "eggs", "eggs_constant")

write_csv(data_eggs_constant_summary, "data-processed/data_eggs_constant_summary.csv")


############## 3B. Lifetime Eggs: DTR 9 ---------------------------------------------------------------

# Analyze TPC model
eggs_dtr9_TPC_analysis <- extractTPC(model_eggs_dtr9, "eggs_dtr9", Temp.xs)
predictions_eggs_dtr9_summary <- eggs_dtr9_TPC_analysis[[1]]
params_eggs_dtr9_summary <- eggs_dtr9_TPC_analysis[[2]]
params_eggs_dtr9_fullposts <- eggs_dtr9_TPC_analysis[[3]]

write_csv(predictions_eggs_dtr9_summary, "data-processed/predictions_eggs_dtr9_summary.csv")
write_csv(params_eggs_dtr9_summary, "data-processed/params_eggs_dtr9_summary.csv")
write_csv(params_eggs_dtr9_fullposts, "data-processed/params_eggs_dtr9_fullposts.csv")

# Process trait data for plotting
data_eggs_dtr9_summary <- processTraitData(data_eggs_dtr9, "eggs", "eggs_dtr9")

write_csv(data_eggs_dtr9_summary, "data-processed/data_eggs_dtr9_summary.csv")


############## 3C. Lifetime Eggs: DTR 12 ---------------------------------------------------------------

# Analyze TPC model
eggs_dtr12_TPC_analysis <- extractTPC(model_eggs_dtr12, "eggs_dtr12", Temp.xs)
predictions_eggs_dtr12_summary <- eggs_dtr12_TPC_analysis[[1]]
params_eggs_dtr12_summary <- eggs_dtr12_TPC_analysis[[2]]
params_eggs_dtr12_fullposts <- eggs_dtr12_TPC_analysis[[3]]

write_csv(predictions_eggs_dtr12_summary, "data-processed/predictions_eggs_dtr12_summary.csv")
write_csv(params_eggs_dtr12_summary, "data-processed/params_eggs_dtr12_summary.csv")
write_csv(params_eggs_dtr12_fullposts, "data-processed/params_eggs_dtr12_fullposts.csv")

# Process trait data for plotting
data_eggs_dtr12_summary <- processTraitData(data_eggs_dtr12, "eggs", "eggs_dtr12")

write_csv(data_eggs_dtr12_summary, "data-processed/data_eggs_dtr12_summary.csv")


############## 4. vector competence (bc) -------------------------------------------------------

# Analyze TPC model
bc_TPC_analysis <- extractTPC(model_bc, "bc", Temp.xs)
predictions_bc_summary <- bc_TPC_analysis[[1]]
params_bc_summary <- bc_TPC_analysis[[2]]
params_bc_fullposts <- bc_TPC_analysis[[3]]

write_csv(predictions_bc_summary, "data-processed/predictions_bc_summary.csv")
write_csv(params_bc_summary, "data-processed/params_bc_summary.csv")
write_csv(params_bc_fullposts, "data-processed/params_bc_fullposts.csv")

# Process trait data for plotting
data_bc_summary <- processTraitData(data_bc, "bc", "bc")

write_csv(data_bc_summary, "data-processed/data_bc_summary.csv")


############## 5. EIP50 ---------------------------------------------------------------

# Analyze TPC model
eip50_TPC_analysis <- extractTPC(model_eip50, "eip50", Temp.xs)
predictions_eip50_summary <- eip50_TPC_analysis[[1]]
params_eip50_summary <- eip50_TPC_analysis[[2]]
params_eip50_fullposts <- eip50_TPC_analysis[[3]]

write_csv(predictions_eip50_summary, "data-processed/predictions_eip50_summary.csv")
write_csv(params_eip50_summary, "data-processed/params_eip50_summary.csv")
write_csv(params_eip50_fullposts, "data-processed/params_eip50_fullposts.csv")

# Process trait data for plotting
data_eip50_summary <- processTraitData(data_eip50, "eip50", "eip50")

write_csv(data_eip50_summary, "data-processed/data_eip50_summary.csv")


############## 6. pEA -----------------------------------------------------------------

# Analyze TPC model
pea_TPC_analysis <- extractTPC(model_pea, "pea", Temp.xs)
predictions_pea_summary <- pea_TPC_analysis[[1]]
params_pea_summary <- pea_TPC_analysis[[2]]
params_pea_fullposts <- pea_TPC_analysis[[3]]

write_csv(predictions_pea_summary, "data-processed/predictions_pea_summary.csv")
write_csv(params_pea_summary, "data-processed/params_pea_summary.csv")
write_csv(params_pea_fullposts, "data-processed/params_pea_fullposts.csv")

# Process trait data for plotting
data_pea_summary <- processTraitData(data_pea, "pea", "pea")

write_csv(data_pea_summary, "data-processed/data_pea_summary.csv")


############## 7. MDR -----------------------------------------------------------------

# Analyze TPC model
mdr_TPC_analysis <- extractTPC(model_mdr, "mdr", Temp.xs)
predictions_mdr_summary <- mdr_TPC_analysis[[1]]
params_mdr_summary <- mdr_TPC_analysis[[2]]
params_mdr_fullposts <- mdr_TPC_analysis[[3]]

write_csv(predictions_mdr_summary, "data-processed/predictions_mdr_summary.csv")
write_csv(params_mdr_summary, "data-processed/params_mdr_summary.csv")
write_csv(params_mdr_fullposts, "data-processed/params_mdr_fullposts.csv")

# Process trait data for plotting
data_mdr_summary <- processTraitData(data_mdr, "mdr", "mdr")

write_csv(data_mdr_summary, "data-processed/data_mdr_summary.csv")


############## 8A. Gamma: constant ---------------------------------------------------------------

# Analyze TPC model
gamma_constant_TPC_analysis <- extractTPC(model_gamma_constant, "gamma_constant", Temp.xs)
predictions_gamma_constant_summary <- gamma_constant_TPC_analysis[[1]]
params_gamma_constant_summary <- gamma_constant_TPC_analysis[[2]]
params_gamma_constant_fullposts <- gamma_constant_TPC_analysis[[3]]

write_csv(predictions_gamma_constant_summary, "data-processed/predictions_gamma_constant_summary.csv")
write_csv(params_gamma_constant_summary, "data-processed/params_gamma_constant_summary.csv")
write_csv(params_gamma_constant_fullposts, "data-processed/params_gamma_constant_fullposts.csv")

# Process trait data for plotting
data_gamma_constant_summary <- processTraitData(data_gamma_constant, "gamma", "gamma_constant")

write_csv(data_gamma_constant_summary, "data-processed/data_gamma_constant_summary.csv")


############## 8B. Gamma: DTR 9 ---------------------------------------------------------------

# Analyze TPC model
gamma_dtr9_TPC_analysis <- extractTPC(model_gamma_dtr9, "gamma_dtr9", Temp.xs)
predictions_gamma_dtr9_summary <- gamma_dtr9_TPC_analysis[[1]]
params_gamma_dtr9_summary <- gamma_dtr9_TPC_analysis[[2]]
params_gamma_dtr9_fullposts <- gamma_dtr9_TPC_analysis[[3]]

write_csv(predictions_gamma_dtr9_summary, "data-processed/predictions_gamma_dtr9_summary.csv")
write_csv(params_gamma_dtr9_summary, "data-processed/params_gamma_dtr9_summary.csv")
write_csv(params_gamma_dtr9_fullposts, "data-processed/params_gamma_dtr9_fullposts.csv")

# Process trait data for plotting
data_gamma_dtr9_summary <- processTraitData(data_gamma_dtr9, "gamma", "gamma_dtr9")

write_csv(data_gamma_dtr9_summary, "data-processed/data_gamma_dtr9_summary.csv")


############## 8C. Gamma: DTR 12 ---------------------------------------------------------------

# Analyze TPC model
gamma_dtr12_TPC_analysis <- extractTPC(model_gamma_dtr12, "gamma_dtr12", Temp.xs)
predictions_gamma_dtr12_summary <- gamma_dtr12_TPC_analysis[[1]]
params_gamma_dtr12_summary <- gamma_dtr12_TPC_analysis[[2]]
params_gamma_dtr12_fullposts <- gamma_dtr12_TPC_analysis[[3]]

write_csv(predictions_gamma_dtr12_summary, "data-processed/predictions_gamma_dtr12_summary.csv")
write_csv(params_gamma_dtr12_summary, "data-processed/params_gamma_dtr12_summary.csv")
write_csv(params_gamma_dtr12_fullposts, "data-processed/params_gamma_dtr12_fullposts.csv")

# Process trait data for plotting
data_gamma_dtr12_summary <- processTraitData(data_gamma_dtr12, "gamma", "gamma_dtr12")

write_csv(data_gamma_dtr12_summary, "data-processed/data_gamma_dtr12_summary.csv")
