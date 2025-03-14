## Rate Summation Project
## JAGS models for fitting thermal performance curves (TPCs)
## Written by Marta Shocket in 2022-24

# NOTES:
# - Running the code below writes a .txt file for each model to your current working directory.
#   	It must be saved there in order for JAGS to find it.
# - The models include a section for 'derived quantities and predictions' that calculates the 
#   	trait across a temperature gradient for each saved sample in the MCMC chain. This output 
#   	is what we plot and use later on to calculate S(T). In this approach The TPC models 
#		take longer to fit, but the code is a bit simpler. 


##### List of Models
	# 1. Truncated normally-distributed Briere - for individual biting rate data
	# 2. Gamma-distributed Quadratic - for individual lifespan data
	# 3. Negative Binomial-distributed Briere - for individual lifetime eggs data
	# 4. Negative Binomial-distributed Quadtratic - for individual lifetime eggs data
	# 5. Normally-distributed Quadratic - for group mean data
	# 6. Normally-distributed Quadratic - for group mean data - derived quantities always =< 1 (for traits that are proportions)
	# 7. Normally-distributed Briere - for group mean data


############## 1. Truncated normally-distributed Briere - for individual biting rate data

sink("briere_T.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()


############## 2. Gamma-distributed Quadratic - for individual lifespan data

sink("quad_gamma.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.ra ~ dunif(0, 100)

    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dgamma(sh[i], cf.ra)
    sh[i] <- cf.ra * trait.mu[i]
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])

    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()


############## 3. Negative Binomial-distributed Briere - for individual lifetime eggs data

sink("briere_negbin.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.r ~ dunif(0, 100)

    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dnegbin(p[i], cf.r)
    p[i] <- cf.r / (cf.r + trait.mu[i])
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])

    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()


############## 4. Negative Binomial-distributed Quadtratic - for individual lifetime eggs data

sink("quad_negbin.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.r ~ dunif(0, 100)

    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dnegbin(p[i], cf.r)
    p[i] <- cf.r / (cf.r + trait.mu[i])
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])

    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()


############## 5. Normally-distributed Quadratic - for group mean data

sink("quad.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()

############## 6. Normally-distributed Quadratic - for group mean data - derived quantities always =< 1 (for traits that are proportions)

sink("quadprob.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }

    } # close model
    ",fill=T)
sink()

############## 7. Normally-distributed Briere - for group mean data

sink("briere.txt")
cat("
    model{

    ## Priors
    cf.q ~ dunif(0, 10)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }

    } # close model
    ",fill=T)
sink()