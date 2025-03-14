## Rate Summation Project
## Functions to load for use in other scripts
## Written by Marta Shocket in 2022-24

##	Table of contents:
##
##	1. Functions to process TPC model output
##			A. Function to calculate TPC posterior summary statistics across a temperature gradient
##			B. Function to extract full posterior distributions for 3 mean-defining TPC parameters & calculate Tbreadth
##			C. Function to calculate Topt
##			D. Function to calculate Tmin, Max, and Tbreadth for derived TPCs 
##			E. Wrapper function to calculate summary data for and extract parameter posteriors from JAGS fitted TPCs
##			F. Wrapper function to calculate summary data for derived TPCs 
##	2. Function to process trait data for plotting
##	3. Function to generate hourly temperatures with the Parton-Logan model
##  4. Function to perform rate summation calculations
##  5. Function to calculate relative suitability
##  6. Functions for Sensitivity Analysis
##			A. Function for derivative of Briere thermal response
##			B. Function for derivative of quadratic thermal response
##			C. Function for sensitivity analysis #1 - partial derivatives
##  7. Custom colors for plotting


########################################### 0. Load tidyverse library

library(tidyverse)

########################################### 1. Functions to process TPC model output

###### A. Function to calculate TPC posterior summary statistics across a temperature gradient
calcPostQuants <- function(TPC_predictions, trait_treatment_name, temp_gradient) {
	
	# Reassign column names to the temperature gradient
	colnames(TPC_predictions) <- temp_gradient
	
	output <- TPC_predictions %>% 
		mutate(iteration = rownames(.)) %>% # add column with iteration number 
		pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>% # convert to long format (3 columns: iteration, temp, & trait value)
		group_by(temperature) %>%
		summarise(lowerCI = quantile(trait_value, probs = 0.025),
				  upperCI = quantile(trait_value, probs = 0.975),
				  mean = mean(trait_value),
				  median = median(trait_value)) %>% 
		mutate(temperature = as.numeric(temperature)) %>% # make temperature numeric
		arrange(temperature) %>% # re-order rows by ascending temperature (since grouping made it categorical, it is ordered alphabetical)
		mutate(treatment = trait_treatment_name) # add column with variable + treatment name
	
	return(output) # return output
	
}

###### B. Function to extract full posterior distributions for 3 mean-defining TPC parameters & calculate Tbreadth
getTPCParamFullPosts <- function (TPC_model, trait_treatment_name) {
	
	output <- data.frame(iteration = seq(1,length(TPC_model$BUGSoutput$sims.list$cf.T0),1),
						 cf.T0 = TPC_model$BUGSoutput$sims.list$cf.T0[,1],
						 cf.Tm = TPC_model$BUGSoutput$sims.list$cf.Tm[,1],
						 cf.q = TPC_model$BUGSoutput$sims.list$cf.q[,1],
						 Tbreadth = (TPC_model$BUGSoutput$sims.list$cf.Tm[,1] - TPC_model$BUGSoutput$sims.list$cf.T0[,1]),
						 traitTreatmentName = trait_treatment_name)
	
	return(output)
}

###### C. Function to calculate Topt
calcToptQuants <- function(TPC_predictions, trait_treatment_name, temp_gradient) {
	
	# Reassign column names to the temperature gradient
	colnames(TPC_predictions) <- temp_gradient
	
	output <- TPC_predictions %>% 
		mutate(iteration = rownames(.)) %>% # add column with iteration number 
		pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>% # convert to long format (3 columns: iteration, temp, & trait value)
		group_by(iteration) %>% 
		slice_max(order_by = trait_value, n = 1) %>% # for each iteration, select row with highest value for the trait
		ungroup() %>% 
		mutate(temperature = as.numeric(temperature)) %>% # make temperature numeric
		summarise(mean = mean(temperature),
				  sd = sd(temperature),
				  lowerCI = quantile(temperature, probs = 0.025),
				  upperCI = quantile(temperature, probs = 0.975),
				  median = median(temperature)) %>% 
		mutate(term = "Topt") %>% # add column with calculation type
		mutate(treatment = trait_treatment_name) # add column with variable + treatment name
	
	return(output) # return output
	
}

###### D. Function to calculate posteriors for Tmin, Tmax, and Tbreadth for derived TPCs 
calcDerivedTPCParamPosteriors <- function(TPC_predictions, temp_gradient) {
	
	# Create output data frame
	output_length <- nrow(TPC_predictions)
	output <- data.frame(cf.T0 = numeric(output_length),
						 cf.Tm = numeric(output_length),
						 Tbreadth = numeric(output_length))
	
	# Get length of temp gradient to compare to final index of Tmax
	length_temp_gradient <- length(temp_gradient)
	
	# Loop through each row (MCMC step of input traits)
	for (i in 1:output_length) {
		
		# Create vector of indices where R0 > 0 and get its length
		#index_list <- which(TPC_predictions[i, ] > 0)
		index_list <- which(TPC_predictions[i, ] > 0.001)
		length_index_list <- length(index_list)
		
		# Get Tmin using the index *before first value* in index_list, unless the first index is 1 (i.e. if Tmin is the lowest value in the temperature gradient)
		ifelse(index_list[1] == 1,
			   output$cf.T0[i] <- temp_gradient[1],
			   output$cf.T0[i] <- temp_gradient[index_list[1] - 1])
		
		# Get Tmax using the index *after last value* in index_list, unless the last index is the length of the temperature gradient (i.e. if Tmax is the highest value in the temperature gradient)
		ifelse(index_list[length_index_list] == length_temp_gradient,
			   output$cf.Tm[i] <- temp_gradient[index_list[length_index_list]],
			   output$cf.Tm[i] <- temp_gradient[index_list[length_index_list] + 1])
		
		# Calculate Tbreadth from Tmin and Tmax
		output$Tbreadth[i] <- output$cf.Tm[i] - output$cf.T0[i]
	}
	
	output # return
	
}


###### E. Wrapper function to calculate summary data for and extract parameter posteriors from JAGS fitted TPCs
extractTPC <- function(TPC_model, trait_treatment_name, temp_gradient) {
	
	# Extract predicted trait values over the temperature gradient
	TPC_predictions <- as.data.frame(TPC_model$BUGSoutput$sims.list$z.trait.mu.pred)
	
	# Calculate TPC posterior summary statistics (means & quantiles)
	TPC_pred_summary <- calcPostQuants(TPC_predictions, trait_treatment_name, temp_gradient)
	
	# Extract full posterior distribution for 3 mean-defining TPC parameters + calculate Tbreadth for each iteration
	TPC_param_full_posts <- getTPCParamFullPosts(TPC_model, trait_treatment_name)
	
	# Calculate Tbreadth summary statistics (mean, sd, & quantiles)
	Tbreadth_summary <- data.frame(term = "Tbreadth",
								   mean = mean(TPC_param_full_posts$Tbreadth),
								   sd = sd(TPC_param_full_posts$Tbreadth),
								   lowerCI = quantile(TPC_param_full_posts$Tbreadth, 0.025)[[1]],
								   median =  quantile(TPC_param_full_posts$Tbreadth, 0.5)[[1]],
								   upperCI = quantile(TPC_param_full_posts$Tbreadth, 0.975)[[1]],
								   treatment = trait_treatment_name)
	
	# Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
	Topt_summary <- calcToptQuants(TPC_predictions, trait_treatment_name, temp_gradient)

	# Pull out parameter summary from the fitted model 
	TPC_param_summary <- as.data.frame(TPC_model$BUGSoutput$summary[1:5,]) %>%
		rownames_to_column(var = "term") %>%
		rename(lowerCI = `2.5%`, median = `50%`, upperCI = `97.5%`) %>% # Rename columns so they are easier to reference & can merge with Topt quantiles
		mutate(treatment = trait_treatment_name)

	# Remove unwanted columns (25% quantile, 75% quantile, Rhat, and n.eff)
	TPC_param_summary <- TPC_param_summary %>%
		select(term, mean, sd, lowerCI, median, upperCI, treatment)
	
	# Add Topt and Tbreadth to parameters summary data frame
	TPC_param_summary_all <- bind_rows(TPC_param_summary, Topt_summary, Tbreadth_summary)
	
	# Bundle output in a list
	output_list <- list(TPC_pred_summary, TPC_param_summary_all, TPC_param_full_posts)
	
	return(output_list) # return output
	
}

###### F. Wrapper function to calculate summary data for derived TPCs 
extractDerivedTPC <- function(TPC_predictions, trait_treatment_name, temp_gradient) {
	
	# Calculate Tmin, Tmax, and Tbreadth posteriors
	TPC_param_full_posts <- calcDerivedTPCParamPosteriors(TPC_predictions, temp_gradient)
	
	# Calculate Tmin, Tmax, and Tbreadth summary statistics (mean, sd, & quantiles)
	Tmin_summary <- data.frame(term = "cf.T0",
								mean = mean(TPC_param_full_posts$cf.T0),
								sd = sd(TPC_param_full_posts$cf.T0),
								lowerCI = quantile(TPC_param_full_posts$cf.T0, 0.025)[[1]],
								median =  quantile(TPC_param_full_posts$cf.T0, 0.5)[[1]],
								upperCI = quantile(TPC_param_full_posts$cf.T0, 0.975)[[1]],
								treatment = trait_treatment_name)
	
	Tmax_summary <- data.frame(term = "cf.Tm",
								mean = mean(TPC_param_full_posts$cf.Tm),
								sd = sd(TPC_param_full_posts$cf.Tm),
								lowerCI = quantile(TPC_param_full_posts$cf.Tm, 0.025)[[1]],
								median =  quantile(TPC_param_full_posts$cf.Tm, 0.5)[[1]],
								upperCI = quantile(TPC_param_full_posts$cf.Tm, 0.975)[[1]],
								treatment = trait_treatment_name)
	
	Tbreadth_summary <- data.frame(term = "Tbreadth",
								   mean = mean(TPC_param_full_posts$Tbreadth),
								   sd = sd(TPC_param_full_posts$Tbreadth),
								   lowerCI = quantile(TPC_param_full_posts$Tbreadth, 0.025)[[1]],
								   median =  quantile(TPC_param_full_posts$Tbreadth, 0.5)[[1]],
								   upperCI = quantile(TPC_param_full_posts$Tbreadth, 0.975)[[1]],
								   treatment = trait_treatment_name)
	
	# Calculate Topt for each iteration and calculate summary statistics (mean, sd, & quantiles)
	Topt_summary <- calcToptQuants(TPC_predictions, trait_treatment_name, temp_gradient)
	
	# Add Topt and Tbreadth to parameters summary data frame
	TPC_param_summary_all <- bind_rows(Tmin_summary, Tmax_summary, Topt_summary, Tbreadth_summary)
	
	return(TPC_param_summary_all) # return output
	
}


########################################### 2. Function to process trait data for plotting
processTraitData <- function (data_input, trait_name, trait_treatment_name) {
	
	output <- data_input %>%
		rename("temperature" = "T") %>% 
		group_by(temperature) %>% 
		summarise(mean = mean(trait),
				  std_error = sd(trait)/sqrt(n())) %>% 
		mutate(trait = trait_name) %>% 
		mutate(treatment = trait_treatment_name)
	
	return(output)
}


########################################### 3. Function to generate hourly temperatures with the Parton-Logan model

### Formula taken from the Murdock Lab Fluctuation Program Calculator Excel Sheet (Creator: Krijn Paaijmans)
### A sinusoial relationship during the day (starting at sunrise) followed by an exponential decay at night (start at sunset)

# NOTE: For the rate summation calculations below to work, the temperature gradient in the hourly temperatures generated here
#		by the L-P function must match the temperature gradient in the TPC predictions (i.e., in 01_TraitTPCFitting.R)

LoganPartonCalc <- function(dtr, mean_temp) {
	
	##################### Set parameters
	#      NOTE: Uses military time 0/24 = midnight; 1 = 1 AM; 14 = 2 PM, etc.
	#      NOTE: For the exp decay function to work, hours must continue to count up (past 23) and do not reset at midnight
	
	# Set parameters for day/night length, start times, and sin wave amplitude
	day_length <- 12
	night_length <- 24 - day_length
	time_sunrise = 6
	time_sunset = time_sunrise + day_length
	amplitude <- dtr/2
	
	# Set other constants
	nocturnal_constant_tau = 4
	time_Tmax_afterNoon = 1.5
	correction_factor <- -0.0575824  # Coefficient for difference btw mean and median temp is specific for day_length = 12
	
	##################### Create hour & temperature data frames for doing math
	
	# Create vectors of hours based on sunrise/sunset times
	sin_hours <- seq(time_sunrise, (time_sunset - 1)  , 1) # hours that follow sin function (starts at sunrise)
	exp_hours <- seq(time_sunset,  (time_sunrise + 23), 1) # hours that follow exp decay function (starts at sunset)
	
	# Create empty matrices for hours input
	sin_hour_input <- matrix(nrow = length(sin_hours), ncol = length(mean_temp))
	exp_hour_input <- matrix(nrow = length(exp_hours), ncol = length(mean_temp))
	
	# Fill them with the appropriate hours
	sin_hour_input[, 1:length(mean_temp)] <- sin_hours
	exp_hour_input[, 1:length(mean_temp)] <- exp_hours
	
	# Turn mean_temps into a matrix with one row, so we can rep() it below
	mean_temp_input <- matrix(nrow = 1, ncol = length(mean_temp))
	mean_temp_input[1, ] <- mean_temp
	
	# Create (filled) matrices for temp input
	sin_mean_temp_input <- mean_temp_input[rep(1, length(sin_hours)), ]
	exp_mean_temp_input <- mean_temp_input[rep(1, length(exp_hours)), ]
	
	# Create (filled) matrices for Tmedian, Tmin, Tmax, and Tsunset - need two versions/sizes for all except Tsunset 
	sin_median_temp <- sin_mean_temp_input - correction_factor*dtr
	sin_Tmin <- sin_median_temp - amplitude
	sin_Tmax <- sin_median_temp + amplitude
	
	exp_median_temp <- exp_mean_temp_input - correction_factor*dtr
	exp_Tmin <- exp_median_temp - amplitude
	exp_Tmax <- exp_median_temp + amplitude
	exp_Tsunset = exp_Tmin + (exp_Tmax-exp_Tmin) * sin(pi*day_length / (day_length + 2*time_Tmax_afterNoon) ) 
	
	##################### Do the math & process the output
	
	# Calculate sin function with matrices
	sin_output <- sin_Tmin + (sin_Tmax-sin_Tmin) * sin(pi*(sin_hour_input - 12 + day_length/2) / (day_length + 2*time_Tmax_afterNoon))
	
	# Calculate exponential decay function with matrices
	exp_output <- (exp_Tmin - exp_Tsunset*exp(-night_length/nocturnal_constant_tau) + (exp_Tsunset-exp_Tmin)*exp(-(exp_hour_input-time_sunset)/nocturnal_constant_tau)) / 
		(1 - exp(-night_length/nocturnal_constant_tau))  
	
	# Combine sin and exponential function outputs into one matrix
	Logan_Parton_output <- rbind(sin_output, exp_output)
	
	# Round to nearest 0.1C and convert to df
	Logan_Parton_output_rounded = round(as.data.frame(Logan_Parton_output), 1)
	
	# Add appropriate row and column names
	rownames(Logan_Parton_output_rounded) <- c(sin_hours, exp_hours)
	colnames(Logan_Parton_output_rounded) <- mean_temp
	
	# Return output
	return(Logan_Parton_output_rounded)
	
} # End function


########################################### 4. Function to perform rate summation calculations

# Arguments:
#	1) TPC_predictions = predicted trait values (from constant temperature TPC) across the temperature gradient for every MCMC iteration 
#			df: 451 cols (temperature gradient) x 7500 rows (MCMC iterations)
#   2) timetemps_df = hourly temperatures across the mean temperature gradient for a given dtr treatment
#			df: 451 cols (mean temperature gradient) x 24 rows (24 hourly temperature values)
#   3) temp_grad = temperature gradient values
#			vector: 451 elements
#
# Returns: 
#	1) rs_calc_gradient = trait values predicted by rate summation across a mean temperature gradient for every MCMC iteration
#			df: 451 cols (temperature gradient) x 7500 rows (MCMC iterations)
#


RSCalcTempGrad <- function(TPC_predictions, timetemps_df, temp_grad) {
	
	# Name the columns so the pivot_longer() below will work
	colnames(TPC_predictions) <- temp_grad
	
	# Pivot constant temperature TPC trait predictions into long format:
	#		3 cols (iteration, temperature, & trait_value) x 3.4 million rows (7500 MCMC iterations x 451 temperatures)
	# Note: this is the longest step, do it as few times as possible
	predictions_long <- TPC_predictions %>%
		mutate(iteration = rownames(.)) %>% # add column with iteration number 
		pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value")  %>%
		mutate(temperature = as.numeric(temperature)) # make temperature numeric
	
	# Create empty data frame to store rate summation calculations
	rs_calc_gradient <- as.data.frame(matrix(nrow = nrow(TPC_predictions), ncol = ncol(timetemps_df)))
	
	# Add in a progress bar because this function takes a way
	pb <- progress_bar$new(
		format = "[:bar] :percent eta: :eta",
		total = ncol(timetemps_df)
	)
	
	# Loop through each mean temperature along the gradient
	for (i in 1:ncol(timetemps_df)) {
		
		# Increment progress bar
		pb$tick()
		
		# Pull out the hourly temperatures for that mean temperature as a df with one col
		# NOTE: the column name must match the pivoted trait predictions column for the left_join to work
		temp_seq <- data.frame(temperature = timetemps_df[, i])
		
		# Left join temperature sequence and trait prediction:
		#		for every row in temp_seq, it adds every trait value + iteration combo from predictions_long with that temperature value
		#		-> df with 3 cols (temperature, iteration, trait value) and 180,000 rows (7500 MCMC x 24 hourly temperatures)
		sj <- left_join(temp_seq, predictions_long, by = "temperature", relationship = "many-to-many")
		
		# Rate summation calculation averages trait values at all hourly temperatures for a given MCMC iteration
		#		-> df with 2 cols (iteration, rs_calc) and 7500 rows (MCMC)
		rs_calc <- sj %>%
			group_by(iteration) %>%
			summarise(rs_calc = mean(trait_value))
		
		# Store rate summation calculation in the corresponding column
		rs_calc_gradient[, i] <- pull(rs_calc[, 2])
		
	} # End loop through mean temperature gradient
	
	return(rs_calc_gradient)
	
} # End function


########################################### 5. Function to calculate relative suitability

calcRelSuit <- function(suitability_calc, comparison_suitability_calc) {
	
	suitability_scale <- max(comparison_suitability_calc$median)
	
	output <- tibble(temperature = suitability_calc$temperature,
					 lowerCI = suitability_calc$lowerCI / suitability_scale,
					 upperCI = suitability_calc$upperCI / suitability_scale,
					 mean = suitability_calc$mean / suitability_scale,
					 median = suitability_calc$median / suitability_scale,
					 treatment = suitability_calc$treatment)
	
	return(output)
	
}


########################################### 6. Functions for Sensitivity Analysis

###### A. Function for derivative of Briere thermal response
d_briere = function(t, T0, Tm, q) {
	
	b <- c()
	
	for (i in 1:length(t)) {
		if (t[i]>T0 && t[i]<Tm) {b[i] <- (q*(-5*(t[i]^2) + 3*t[i]*T0 + 4*t[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-t[i])))}
		else {b[i] <- 0}
	}
	
	b # return output
	
}

###### B. Function for derivative of quadratic thermal response
d_quad = function(t, T0, Tm, q){
	
	b <- c()
	
	for (i in 1:length(t)){
		if (t[i]>T0 && t[i]<Tm) {b[i] <- -1*q*(2*t[i] - T0 - Tm)}
		else {b[i] <- 0}
	}
	
	b # return output
	
}

###### C. Function for sensitivity analysis #1 - partial derivatives

# Arguments: mod_x = the JAGS model for each trait (for the fitted TPC parameters - T0, Tm, and q);
#			 m_x = the mean value for each trait over the temperature gradient

SensitivityAnalysis1_pd = function(mod_a, mod_bc, mod_lf, mod_gamma, mod_B, mod_pea, mod_mdr,
								   m_a, m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr) {

	# Create matrices to hold results
	dR0.da <- dR0.dbc <- dR0.dlf <- dR0.dgamma <- dR0.dB <- dR0.dpEA <- dR0.dMDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)
	
	# Extract predicted trait values
	mod_a_preds <- mod_a$BUGSoutput$sims.list$z.trait.mu.pred
	mod_bc_preds <- mod_bc$BUGSoutput$sims.list$z.trait.mu.pred
	mod_lf_preds <- mod_lf$BUGSoutput$sims.list$z.trait.mu.pred
	mod_gamma_preds <- mod_gamma$BUGSoutput$sims.list$z.trait.mu.pred
	mod_B_preds <- mod_B$BUGSoutput$sims.list$z.trait.mu.pred
	mod_pea_preds <- mod_pea$BUGSoutput$sims.list$z.trait.mu.pred
	mod_mdr_preds <- mod_mdr$BUGSoutput$sims.list$z.trait.mu.pred
	
	# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
	for(i in 1:nMCMC){ # loop through MCMC steps
		
		# Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
		# The sims.list refers to the lists of fitted TPC parameters (T0, Tm, and q)
		da.dT <- d_briere(Temp.xs, 
						  mod_a$BUGSoutput$sims.list[[1]][i],
						  mod_a$BUGSoutput$sims.list[[2]][i],
						  mod_a$BUGSoutput$sims.list[[3]][i])
		dbc.dT <- d_quad(Temp.xs, 
						 mod_bc$BUGSoutput$sims.list[[1]][i],
						 mod_bc$BUGSoutput$sims.list[[2]][i],
						 mod_bc$BUGSoutput$sims.list[[3]][i])
		dlf.dT <- d_quad(Temp.xs, 
						 mod_lf$BUGSoutput$sims.list[[1]][i],
						 mod_lf$BUGSoutput$sims.list[[2]][i],
						 mod_lf$BUGSoutput$sims.list[[3]][i])
		dgamma.dT <- d_quad(Temp.xs,
							mod_gamma$BUGSoutput$sims.list[[1]][i],
							mod_gamma$BUGSoutput$sims.list[[2]][i],
							mod_gamma$BUGSoutput$sims.list[[3]][i])
		dB.dT <- d_quad(Temp.xs,
						  mod_B$BUGSoutput$sims.list[[1]][i],
						  mod_B$BUGSoutput$sims.list[[2]][i],
						  mod_B$BUGSoutput$sims.list[[3]][i])
		dpEA.dT <- d_quad(Temp.xs,
						  mod_pea$BUGSoutput$sims.list[[1]][i],
						  mod_pea$BUGSoutput$sims.list[[2]][i],
						  mod_pea$BUGSoutput$sims.list[[3]][i])
		dMDR.dT <- d_briere(Temp.xs,
							mod_mdr$BUGSoutput$sims.list[[1]][i],
							mod_mdr$BUGSoutput$sims.list[[2]][i],
							mod_mdr$BUGSoutput$sims.list[[3]][i])

		# Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
		# Most parameters are multiplied in R0 equation once, except lf and a (both squared - gamma is used instead of e^[-1/{lf*PDR})
	   	# See Mathematica notebook from Shocket et al. 2018 eLife for dR0/dy derivative calculations
		dR0.da[i, ] <- R0eq(mod_a_preds[i, ], m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr)/(mod_a_preds[i, ]+ec) * da.dT
		dR0.dbc[i, ] <- 1/2 * (R0eq(m_a, mod_bc_preds[i, ], m_lf, m_gamma, m_B, m_pea, m_mdr)/(mod_bc_preds[i, ]+ec) * dbc.dT)
		dR0.dlf[i, ] <- R0eq(m_a, m_bc, mod_lf_preds[i, ], m_gamma, m_B, m_pea, m_mdr)/(mod_lf_preds[i, ]+ec) * dlf.dT
		dR0.dgamma[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, mod_gamma_preds[i, ], m_B, m_pea, m_mdr)/(mod_gamma_preds[i, ]+ec) * dgamma.dT)
		dR0.dB[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, m_gamma, mod_B_preds[i, ], m_pea, m_mdr)/(mod_B_preds[i, ]+ec) * dB.dT)
		dR0.dpEA[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, m_gamma, m_B, mod_pea_preds[i, ], m_mdr)/(mod_pea_preds[i, ]+ec) * dpEA.dT)
		dR0.dMDR[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, m_gamma, m_B, m_pea, mod_mdr_preds[i, ])/(mod_mdr_preds[i, ]+ec) * dMDR.dT)
		dR0.dT[i, ] <-  dR0.da[i, ] + dR0.dbc[i, ] + dR0.dlf[i, ] + dR0.dgamma[i, ] + dR0.dB[i, ] + dR0.dpEA[i, ] + dR0.dMDR[i, ]
		
	} # end MCMC loop
	
	# Collect output in a list and return it
	SA1_list_out <- list(dR0.da, dR0.dbc, dR0.dlf, dR0.dgamma, dR0.dB, dR0.dpEA, dR0.dMDR, dR0.dT)
	SA1_list_out
	
} # end function

###### C. Function for sensitivity analysis #1 - partial derivatives

# Arguments: mod_x = the JAGS model for each trait (for the fitted TPC parameters - T0, Tm, and q);
#			 m_x = the mean value for each trait over the temperature gradient

SensitivityAnalysis1_pd = function(mod_a, mod_bc, mod_lf, mod_gamma, mod_B, mod_pea, mod_mdr,
								   m_a, m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr) {
	
	# Create matrices to hold results
	dR0.da <- dR0.dbc <- dR0.dlf <- dR0.dgamma <- dR0.dB <- dR0.dpEA <- dR0.dMDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)
	
	# Extract predicted trait values
	mod_a_preds <- mod_a$BUGSoutput$sims.list$z.trait.mu.pred
	mod_bc_preds <- mod_bc$BUGSoutput$sims.list$z.trait.mu.pred
	mod_lf_preds <- mod_lf$BUGSoutput$sims.list$z.trait.mu.pred
	mod_gamma_preds <- mod_gamma$BUGSoutput$sims.list$z.trait.mu.pred
	mod_B_preds <- mod_B$BUGSoutput$sims.list$z.trait.mu.pred
	mod_pea_preds <- mod_pea$BUGSoutput$sims.list$z.trait.mu.pred
	mod_mdr_preds <- mod_mdr$BUGSoutput$sims.list$z.trait.mu.pred
	
	# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
	for(i in 1:nMCMC){ # loop through MCMC steps
		
		# Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
		# The sims.list refers to the lists of fitted TPC parameters (T0, Tm, and q)
		da.dT <- d_briere(Temp.xs, 
						  mod_a$BUGSoutput$sims.list[[1]][i],
						  mod_a$BUGSoutput$sims.list[[2]][i],
						  mod_a$BUGSoutput$sims.list[[3]][i])
		dbc.dT <- d_quad(Temp.xs, 
						 mod_bc$BUGSoutput$sims.list[[1]][i],
						 mod_bc$BUGSoutput$sims.list[[2]][i],
						 mod_bc$BUGSoutput$sims.list[[3]][i])
		dlf.dT <- d_quad(Temp.xs, 
						 mod_lf$BUGSoutput$sims.list[[1]][i],
						 mod_lf$BUGSoutput$sims.list[[2]][i],
						 mod_lf$BUGSoutput$sims.list[[3]][i])
		dgamma.dT <- d_quad(Temp.xs,
							mod_gamma$BUGSoutput$sims.list[[1]][i],
							mod_gamma$BUGSoutput$sims.list[[2]][i],
							mod_gamma$BUGSoutput$sims.list[[3]][i])
		dB.dT <- d_quad(Temp.xs,
						mod_B$BUGSoutput$sims.list[[1]][i],
						mod_B$BUGSoutput$sims.list[[2]][i],
						mod_B$BUGSoutput$sims.list[[3]][i])
		dpEA.dT <- d_quad(Temp.xs,
						  mod_pea$BUGSoutput$sims.list[[1]][i],
						  mod_pea$BUGSoutput$sims.list[[2]][i],
						  mod_pea$BUGSoutput$sims.list[[3]][i])
		dMDR.dT <- d_briere(Temp.xs,
							mod_mdr$BUGSoutput$sims.list[[1]][i],
							mod_mdr$BUGSoutput$sims.list[[2]][i],
							mod_mdr$BUGSoutput$sims.list[[3]][i])
		
		# Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
		# Most parameters are multiplied in R0 equation once, except lf and a (both squared - gamma is used instead of e^[-1/{lf*PDR})
		# See Mathematica notebook from Shocket et al. 2018 eLife for dR0/dy derivative calculations
		dR0.da[i, ] <- R0eq(mod_a_preds[i, ], m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr)/(mod_a_preds[i, ]+ec) * da.dT
		dR0.dbc[i, ] <- 1/2 * (R0eq(m_a, mod_bc_preds[i, ], m_lf, m_gamma, m_B, m_pea, m_mdr)/(mod_bc_preds[i, ]+ec) * dbc.dT)
		dR0.dlf[i, ] <- R0eq(m_a, m_bc, mod_lf_preds[i, ], m_gamma, m_B, m_pea, m_mdr)/(mod_lf_preds[i, ]+ec) * dlf.dT
		dR0.dgamma[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, mod_gamma_preds[i, ], m_B, m_pea, m_mdr)/(mod_gamma_preds[i, ]+ec) * dgamma.dT)
		dR0.dB[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, m_gamma, mod_B_preds[i, ], m_pea, m_mdr)/(mod_B_preds[i, ]+ec) * dB.dT)
		dR0.dpEA[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, m_gamma, m_B, mod_pea_preds[i, ], m_mdr)/(mod_pea_preds[i, ]+ec) * dpEA.dT)
		dR0.dMDR[i, ] <- 1/2 * (R0eq(m_a, m_bc, m_lf, m_gamma, m_B, m_pea, mod_mdr_preds[i, ])/(mod_mdr_preds[i, ]+ec) * dMDR.dT)
		dR0.dT[i, ] <-  dR0.da[i, ] + dR0.dbc[i, ] + dR0.dlf[i, ] + dR0.dgamma[i, ] + dR0.dB[i, ] + dR0.dpEA[i, ] + dR0.dMDR[i, ]
		
	} # end MCMC loop
	
	# Collect output in a list and return it
	SA1_list_out <- list(dR0.da, dR0.dbc, dR0.dlf, dR0.dgamma, dR0.dB, dR0.dpEA, dR0.dMDR, dR0.dT)
	SA1_list_out
	
} # end function


###### D. Function for sensitivity analysis #2 - holding single parameters constant

# Arguments: mod_x_preds = matrix of predicted values for each trait - rows are MCMC iterations, columns are the temperature gradient;
#			 m_x = the mean value for each trait over the temperature gradient

SensitivityAnalysis2_hspc = function(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds) {
	
	# Create matrices to hold results
	dR0.da <- dR0.dbc <- dR0.dlf <- dR0.dgamma <- dR0.dB <- dR0.dpEA <- dR0.dMDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)
	
	# Calculate R0 holding each parameter constant
	SA2_a		<- R0eq(1, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
	SA2_bc		<- R0eq(mod_a_preds, 1, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
	SA2_lf		<- R0eq(mod_a_preds, mod_bc_preds, 1, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
	SA2_gamma	<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, 1, mod_B_preds, mod_pea_preds, mod_mdr_preds)
	SA2_B 		<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, 1, mod_pea_preds, mod_mdr_preds)
	SA2_pEA 	<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, 1, mod_mdr_preds)
	SA2_MDR 	<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, 1)
	SA2_all		<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
	
	# Collect output in a list and return it
	SA2_list_out <- list(SA2_a, SA2_bc, SA2_lf, SA2_gamma, SA2_B, SA2_pEA, SA2_MDR, SA2_all)
	SA2_list_out
	
} # end function


###### Original version using models instead of matrix predicted values as the input

# ###### D. Function for sensitivity analysis #2 - holding single parameters constant
# 
# # Arguments: mod_x = the JAGS model for each trait (for the fitted TPC parameters - T0, Tm, and q);
# #			 m_x = the mean value for each trait over the temperature gradient
# 
# SensitivityAnalysis2_hspc = function(mod_a, mod_bc, mod_lf, mod_gamma, mod_B, mod_pea, mod_mdr) {
# 	
# 	# Create matrices to hold results
# 	dR0.da <- dR0.dbc <- dR0.dlf <- dR0.dgamma <- dR0.dB <- dR0.dpEA <- dR0.dMDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)
# 	
# 	# Extract predicted trait values
# 	mod_a_preds <- mod_a$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_bc_preds <- mod_bc$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_lf_preds <- mod_lf$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_gamma_preds <- mod_gamma$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_B_preds <- mod_B$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_pea_preds <- mod_pea$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_mdr_preds <- mod_mdr$BUGSoutput$sims.list$z.trait.mu.pred
# 	
# 	# Calculate R0 holding each parameter constant
# 	SA2_a		<- R0eq(1, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
# 	SA2_bc		<- R0eq(mod_a_preds, 1, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
# 	SA2_lf		<- R0eq(mod_a_preds, mod_bc_preds, 1, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
# 	SA2_gamma	<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, 1, mod_B_preds, mod_pea_preds, mod_mdr_preds)
# 	SA2_B 		<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, 1, mod_pea_preds, mod_mdr_preds)
# 	SA2_pEA 	<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, 1, mod_mdr_preds)
# 	SA2_MDR 	<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, 1)
# 	SA2_all		<- R0eq(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds)
# 	
# 	# Collect output in a list and return it
# 	SA2_list_out <- list(SA2_a, SA2_bc, SA2_lf, SA2_gamma, SA2_B, SA2_pEA, SA2_MDR, SA2_all)
# 	SA2_list_out
# 	
# } # end function


###### E. Functions for uncertainty analysis

### Standard version of Uncertainty Analysis for Models 1-4

# Arguments: mod_x_preds = matrix of predicted values for each trait - rows are MCMC iterations, columns are the temperature gradient;
#			 m_x = the mean value for each trait over the temperature gradient

UncertaintyAnalysis = function(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds,
								   m_a, m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr) {
	
	# Create matrices to hold results
	R0.a <- R0.bc <- R0.lf <- R0.gamma <- R0.B <- R0.pEA <- R0.MDR <- R0.full <- matrix(NA, nMCMC, N.Temp.xs)
	
	# Calculate posterior samples for R0 across the temp gradient with all-but-one trait fixed to the posterior mean
	for(j in 1:nMCMC){ # loop through MCMC steps
		
		R0.a[j,] <- R0eq(mod_a_preds[j,], m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr)
		R0.bc[j,] <- R0eq(m_a, mod_bc_preds[j,], m_lf, m_gamma, m_B, m_pea, m_mdr)
		R0.lf[j,] <- R0eq(m_a, m_bc, mod_lf_preds[j,], m_gamma, m_B, m_pea, m_mdr)
		R0.gamma[j,] <- R0eq(m_a, m_bc, m_lf, mod_gamma_preds[j,], m_B, m_pea, m_mdr)
		R0.B[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, mod_B_preds[j,], m_pea, m_mdr)
		R0.pEA[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, m_B, mod_pea_preds[j,], m_mdr)
		R0.MDR[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, m_B, m_pea, mod_mdr_preds[j,])
		R0.full[j,] <- R0eq(mod_a_preds[j,], mod_bc_preds[j,], mod_lf_preds[j,], mod_gamma_preds[j,], mod_B_preds[j,], mod_pea_preds[j,], mod_mdr_preds[j,])
		
	}
	
	## Calculate the width of the inner 95% quantile for the posteriors from above
	a.q 	<- apply(R0.a, 2, FUN=quantile, probs=0.925) - apply(R0.a, 2, FUN=quantile, probs=0.025)
	bc.q	<- apply(R0.bc, 2, FUN=quantile, probs=0.925) - apply(R0.bc, 2, FUN=quantile, probs=0.025)
	lf.q	<- apply(R0.lf, 2, FUN=quantile, probs=0.925) - apply(R0.lf, 2, FUN=quantile, probs=0.025)
	gamma.q <- apply(R0.gamma, 2, FUN=quantile, probs=0.925) - apply(R0.gamma, 2, FUN=quantile, probs=0.025)
	B.q 	<- apply(R0.B, 2, FUN=quantile, probs=0.925) - apply(R0.B, 2, FUN=quantile, probs=0.025)
	pEA.q	<- apply(R0.pEA, 2, FUN=quantile, probs=0.925) - apply(R0.pEA, 2, FUN=quantile, probs=0.025)
	MDR.q	<- apply(R0.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.MDR, 2, FUN=quantile, probs=0.025)
	R0.q	<- apply(R0.full, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.full, 2, FUN=quantile, probs=0.025, na.rm=F)
	
	UA_df_out <- data.frame(a = a.q, bc = bc.q, lf = lf.q, gamma = gamma.q, B = B.q, pEA = pEA.q, MDR = MDR.q, R0 = R0.q)
	UA_df_out
	
} # end function

### Rate summation version Uncertainty Analysis for Models 5

# Arguments: mod_x_preds = matrix of predicted values for each trait - rows are MCMC iterations, columns are the temperature gradient;
#			 m_x = the mean value for each trait over the temperature gradient
#			 timetemps_df = hourly temperatures across the mean temperature gradient for a given dtr treatment
#					df: 451 cols (mean temperature gradient) x 24 rows (24 hourly temperature values)
#			 temp_grad = temperature gradient values
#					vector: 451 elements

UncertaintyAnalysis_RSModel5 = function(mod_a_preds, mod_bc_preds, mod_lf_preds, mod_gamma_preds, mod_B_preds, mod_pea_preds, mod_mdr_preds,
							   m_a, m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr, timetemps_df, temp_grad) {
	
	# Create matrices to hold results
	R0.a <- R0.bc <- R0.lf <- R0.gamma <- R0.B <- R0.pEA <- R0.MDR <- R0.full <- matrix(NA, nMCMC, N.Temp.xs)
	
	# Calculate posterior samples for R0 across the temp gradient with all-but-one trait fixed to the posterior mean
	for(j in 1:nMCMC){ # loop through MCMC steps
		
		R0.a[j,] <- R0eq(mod_a_preds[j,], m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr)
		R0.bc[j,] <- R0eq(m_a, mod_bc_preds[j,], m_lf, m_gamma, m_B, m_pea, m_mdr)
		R0.lf[j,] <- R0eq(m_a, m_bc, mod_lf_preds[j,], m_gamma, m_B, m_pea, m_mdr)
		R0.gamma[j,] <- R0eq(m_a, m_bc, m_lf, mod_gamma_preds[j,], m_B, m_pea, m_mdr)
		R0.B[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, mod_B_preds[j,], m_pea, m_mdr)
		R0.pEA[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, m_B, mod_pea_preds[j,], m_mdr)
		R0.MDR[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, m_B, m_pea, mod_mdr_preds[j,])
		R0.full[j,] <- R0eq(mod_a_preds[j,], mod_bc_preds[j,], mod_lf_preds[j,], mod_gamma_preds[j,], mod_B_preds[j,], mod_pea_preds[j,], mod_mdr_preds[j,])
		
	}
	
	# NOTE: Need to change the matrices to data frames for dplyr in the RSCalcTempGrad function
	R0.a_rs <- RSCalcTempGrad(as.data.frame(R0.a), timetemps_df, temp_grad)
	R0.bc_rs <- RSCalcTempGrad(as.data.frame(R0.bc), timetemps_df, temp_grad)
	R0.lf_rs <- RSCalcTempGrad(as.data.frame(R0.lf), timetemps_df, temp_grad)
	R0.gamma_rs <- RSCalcTempGrad(as.data.frame(R0.gamma), timetemps_df, temp_grad)
	R0.B_rs <- RSCalcTempGrad(as.data.frame(R0.B), timetemps_df, temp_grad)
	R0.pEA_rs <- RSCalcTempGrad(as.data.frame(R0.pEA), timetemps_df, temp_grad)
	R0.MDR_rs <- RSCalcTempGrad(as.data.frame(R0.MDR), timetemps_df, temp_grad)
	R0.full_rs <- RSCalcTempGrad(as.data.frame(R0.full), timetemps_df, temp_grad)
	
	## Calculate the width of the inner 95% quantile for the posteriors from above - works with data frames or matrices, so can stay the same
	a.q 	<- apply(R0.a_rs, 2, FUN=quantile, probs=0.925) - apply(R0.a_rs, 2, FUN=quantile, probs=0.025)
	bc.q	<- apply(R0.bc_rs, 2, FUN=quantile, probs=0.925) - apply(R0.bc_rs, 2, FUN=quantile, probs=0.025)
	lf.q	<- apply(R0.lf_rs, 2, FUN=quantile, probs=0.925) - apply(R0.lf_rs, 2, FUN=quantile, probs=0.025)
	gamma.q <- apply(R0.gamma_rs, 2, FUN=quantile, probs=0.925) - apply(R0.gamma_rs, 2, FUN=quantile, probs=0.025)
	B.q 	<- apply(R0.B_rs, 2, FUN=quantile, probs=0.925) - apply(R0.B_rs, 2, FUN=quantile, probs=0.025)
	pEA.q	<- apply(R0.pEA_rs, 2, FUN=quantile, probs=0.925) - apply(R0.pEA_rs, 2, FUN=quantile, probs=0.025)
	MDR.q	<- apply(R0.MDR_rs, 2, FUN=quantile, probs=0.925) - apply(R0.MDR_rs, 2, FUN=quantile, probs=0.025)
	R0.q	<- apply(R0.full_rs, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.full_rs, 2, FUN=quantile, probs=0.025, na.rm=F)
	
	UA_df_out <- data.frame(a = a.q, bc = bc.q, lf = lf.q, gamma = gamma.q, B = B.q, pEA = pEA.q, MDR = MDR.q, R0 = R0.q)
	UA_df_out
	
} # end function

###### Original version using models instead of matrix predicted values as the input

# ###### E. Function for uncertainty analysis
# 
# # Arguments: mod_x = the JAGS model for each trait (for the fitted TPC parameters - T0, Tm, and q);
# #			 m_x = the mean value for each trait over the temperature gradient
# 
# UncertaintyAnalysis = function(mod_a, mod_bc, mod_lf, mod_gamma, mod_B, mod_pea, mod_mdr,
# 							   m_a, m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr) {
# 	
# 	# Create matrices to hold results
# 	R0.a <- R0.bc <- R0.lf <- R0.gamma <- R0.B <- R0.pEA <- R0.MDR <- R0.full <- matrix(NA, nMCMC, N.Temp.xs)
# 	
# 	# Extract predicted trait values
# 	mod_a_preds <- mod_a$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_bc_preds <- mod_bc$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_lf_preds <- mod_lf$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_gamma_preds <- mod_gamma$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_B_preds <- mod_B$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_pea_preds <- mod_pea$BUGSoutput$sims.list$z.trait.mu.pred
# 	mod_mdr_preds <- mod_mdr$BUGSoutput$sims.list$z.trait.mu.pred
# 	
# 	# Calculate posterior samples for R0 across the temp gradient with all-but-one trait fixed to the posterior mean
# 	for(j in 1:nMCMC){ # loop through MCMC steps
# 		
# 		R0.a[j,] <- R0eq(mod_a_preds[j,], m_bc, m_lf, m_gamma, m_B, m_pea, m_mdr)
# 		R0.bc[j,] <- R0eq(m_a, mod_bc_preds[j,], m_lf, m_gamma, m_B, m_pea, m_mdr)
# 		R0.lf[j,] <- R0eq(m_a, m_bc, mod_lf_preds[j,], m_gamma, m_B, m_pea, m_mdr)
# 		R0.gamma[j,] <- R0eq(m_a, m_bc, m_lf, mod_gamma_preds[j,], m_B, m_pea, m_mdr)
# 		R0.B[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, mod_B_preds[j,], m_pea, m_mdr)
# 		R0.pEA[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, m_B, mod_pea_preds[j,], m_mdr)
# 		R0.MDR[j,] <- R0eq(m_a, m_bc, m_lf, m_gamma, m_B, m_pea, mod_mdr_preds[j,])
# 		R0.full[j,] <- R0eq(mod_a_preds[j,], mod_bc_preds[j,], mod_lf_preds[j,], mod_gamma_preds[j,], mod_B_preds[j,], mod_pea_preds[j,], mod_mdr_preds[j,])
# 		
# 	}
# 	
# 	## Calculate the width of the inner 95% quantile for the posteriors from above
# 	a.q 	<- apply(R0.a, 2, FUN=quantile, probs=0.925) - apply(R0.a, 2, FUN=quantile, probs=0.025)
# 	bc.q	<- apply(R0.bc, 2, FUN=quantile, probs=0.925) - apply(R0.bc, 2, FUN=quantile, probs=0.025)
# 	lf.q	<- apply(R0.lf, 2, FUN=quantile, probs=0.925) - apply(R0.lf, 2, FUN=quantile, probs=0.025)
# 	gamma.q <- apply(R0.gamma, 2, FUN=quantile, probs=0.925) - apply(R0.gamma, 2, FUN=quantile, probs=0.025)
# 	B.q 	<- apply(R0.B, 2, FUN=quantile, probs=0.925) - apply(R0.B, 2, FUN=quantile, probs=0.025)
# 	pEA.q	<- apply(R0.pEA, 2, FUN=quantile, probs=0.925) - apply(R0.pEA, 2, FUN=quantile, probs=0.025)
# 	MDR.q	<- apply(R0.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.MDR, 2, FUN=quantile, probs=0.025)
# 	R0.q	<- apply(R0.full, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.full, 2, FUN=quantile, probs=0.025, na.rm=F)
# 	
# 	UA_df_out <- data.frame(a = a.q, bc = bc.q, lf = lf.q, gamma = gamma.q, B = B.q, pEA = pEA.q, MDR = MDR.q, R0 = R0.q)
# 	UA_df_out
# 	
# } # end function


########################################### 7. Custom colors for plotting

# # We need to manually specify the color brewer colors to across different figures since there are different treatments in different panels
#   and modify them for different comparisons (rs vs. empirical fits)
#
# # Load the Color Brewer package
# library("RColorBrewer")
#
# # Get hex values for the automatically generated colors in Fig 1 using the eyedropper using GIMP:
#
# Try Purple-Green palette
# display.brewer.pal(9, "PRGn")
# brewer.pal(9, "PRGn")
# 
# # Convert the PRGn hex values to rgb
col2rgb("#762A83") # dark purple = c_emp12
col2rgb("#9970AB") # medium purple = c_rstraits12
col2rgb("#C2A5CF") # light purple = c_rssuit12
col2rgb("#A6DBA0") # light green = c_rssuit09
col2rgb("#5AAE61") # medium green = c_rstraits09
col2rgb("#1B7837") # dark green = c_emp09

# Old colours with different values for lines
# Store PRGn color rgb values for plotting points and lines
c_constant <-		rgb(190, 190, 190, max = 255, alpha = 255)
c_emp12 <- 			rgb(118, 42, 131, max = 255, alpha = 255)
c_rstraits12 <- 	rgb(153, 112, 171, max = 255, alpha = 255)
c_rssuit12 <-		rgb(194, 165, 207, max = 255, alpha = 255)
c_rssuit09 <-		rgb(166, 219, 160, max = 255, alpha = 255)
c_rstraits09 <-		rgb(90, 174, 97, max = 255, alpha = 255)
c_emp09 <-			rgb(27, 120, 55, max = 255, alpha = 255)

# Store PRGn color rgb values for plotting ribbons
ct_constant <-		rgb(190, 190, 190, max = 255, alpha = 95)
ct_emp12 <- 		rgb(118, 42, 131, max = 255, alpha = 95)
ct_rstraits12 <- 	rgb(153, 112, 171, max = 255, alpha = 95)
ct_rssuit12 <-		rgb(194, 165, 207, max = 255, alpha = 95)
ct_rssuit09 <-		rgb(166, 219, 160, max = 255, alpha = 95)
ct_rstraits09 <-	rgb(90, 174, 97, max = 255, alpha = 95)
ct_emp09 <-			rgb(27, 120, 55, max = 255, alpha = 95)
