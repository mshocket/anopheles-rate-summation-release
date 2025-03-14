## Anopheles Rate Summation Project - Demo Analysis
## Written by Marta Shocket in 2024

##	Table of contents:
##
##	1. Set-up workspace
##  2. Fit a trait TPC
##	3. Perform rate summation calculation on the TPC
##	4. Process the results
##	5. Plot the results


##########
###### 1. Set-up workspace 
##########

##### Load JAGS libraries
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('progress')

##### Set working directory
# You will need to set your working directory to the top-level anopheles-rate-summation directory.

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
# - The text files should already be present in the top-level working directory, so the source() line below is commented out

# Create .txt files of JAGS models - un-comment the line below if you need to regenerate the model files
# source("R-scripts/working-versions-code/00_JAGSModels.R")

##### Load trait data from constant temperatures to use as an example, change one column name
data.constant <- read.csv("data/constant.individual.trait.csv")
names(data.constant)[names(data.constant) == 'Treatment'] <- 'temp'


##########
###### 2. Fit a trait TPC
##########

###################################### Settings for JAGS 

##### General settings

##### MCMC Settings
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Derived Quantity Settings
Temp.xs <- seq(0, 45, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

##### Settings specific for fitting a normal / truncated normal distribution

#####  inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


###################################### Fit TPC for bite rate at constant temperature

##### Get data
data_bite_rate_constant <- with(data.constant, data.frame('T' = temp, 'trait' = bite.rate))

##### Organize Data for JAGS
trait <- data_bite_rate_constant$trait
N.obs <- length(trait)
temp <- data_bite_rate_constant$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file** - Briere function, truncated normal distribution
model_bite_rate_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
								 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_constant$BUGSoutput$summary[1:5,]
mcmcplot(model_bite_rate_constant)

# Extract the DIC for future model comparisons
model_bite_rate_constant$BUGSoutput$DIC

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_constant, 
	 ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 3. Perform rate summation calculation on the TPC
##########

######################################  Generate hourly temperatures for a fluctuating environment with the Parton-Logan model

# Temperature gradient that matches TPC predictions
Temp.xs <- seq(0, 45, 0.1)

# Generate hourly temperature sequences across the temperature gradient
LPtemps_dtr9 <- LoganPartonCalc(dtr = 9, Temp.xs)
LPtemps_dtr12 <- LoganPartonCalc(dtr = 12, Temp.xs)

# Set negative temperature values to 0, since TPCs stop at 0 on the low end
LPtemps_dtr9[LPtemps_dtr9 < 0 ] <- 0
LPtemps_dtr12[LPtemps_dtr12 < 0 ] <- 0

# Set temperature values > 45 to 45, since TPCs stop at 45 on the high end
LPtemps_dtr9[LPtemps_dtr9 > 45 ] <- 45
LPtemps_dtr12[LPtemps_dtr12 > 45 ] <- 45

######################################  Apply rate summation

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_bite_rate_constant <- as.data.frame(model_bite_rate_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_bite_rate_rs9 <- RSCalcTempGrad(predictions_bite_rate_constant, LPtemps_dtr9, Temp.xs)
predictions_bite_rate_rs12 <- RSCalcTempGrad(predictions_bite_rate_constant, LPtemps_dtr12, Temp.xs)


##########
###### 4. Process the results
##########

# Process the observed trait data for plotting as points
data_bite_rate_constant_summary <- processTraitData(data_bite_rate_constant, "bite rate", "bite_rate_constant")

### Process the TPC model
bite_rate_constant_TPC_analysis <- extractTPC(model_bite_rate_constant, "bite_rate_constant", Temp.xs)
# Extract a data frame of mean, median, lower CI, and upper CI across the temperature gradient
predictions_bite_rate_constant_summary <- bite_rate_constant_TPC_analysis[[1]]
# Get summary statistics of Tmin, Tmax, Topt, and Tbreadth and other parameters
params_bite_rate_constant_summary <- bite_rate_constant_TPC_analysis[[2]]

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_bite_rate_rs9_summary <- calcPostQuants(predictions_bite_rate_rs9, "bite_rate_rs09", Temp.xs)
predictions_bite_rate_rs12_summary <- calcPostQuants(predictions_bite_rate_rs12, "bite_rate_rs12", Temp.xs)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_bite_rate_rs9_summary <- extractDerivedTPC(predictions_bite_rate_rs9, "bite_rate_rs09", Temp.xs)
params_bite_rate_rs12_summary <- extractDerivedTPC(predictions_bite_rate_rs12, "bite_rate_rs12", Temp.xs)


##########
###### 5. Plot the results
##########

### combine empirically measured constant and rate summation fluctuation treatments
bite_rate_demo_predictions <- bind_rows(predictions_bite_rate_constant_summary, predictions_bite_rate_rs9_summary, predictions_bite_rate_rs12_summary)

bite_rate_rs9_plot <- bite_rate_demo_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Bite~rate~(day^-1)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2.5, 1, 1, 2.5), "lines"), axis.title.y = element_text(vjust = -2)) +
	coord_cartesian(clip = "off") +
	annotate("text", x = 0, y = 0.53, label = "A", size = 5)
