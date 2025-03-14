# README 
Code and Software information for:
*'Mean daily temperatures predict the thermal limits of malaria transmission better than hourly rate summation'*, authored by:
- Marta S. Shocket (m.shocket@lancaster.ac.uk),
- Joey R. Bernhardt (joey.bernhardt@uoguelph.ca),
- Kerri L. Miazgowicz (NCD0@cdc.gov),
- Alyzeh Orakzai (alyzeh.orakzai@lmunet.edu),
- Van M. Savage (vsavage@ucla.edu),
- Richard J. Hall (rjhall@uga.edu),
- Sadie J. Ryan (sjryan@ufl.edu),
- Courtney C. Murdock (ccm256@cornell.edu).

## 1. System Requirements

Our custom data analysis pipeline should run on any relatively recent Mac, PC, or Linux operating system. No specialized hardware is required.

Please use a version of R ≥ v4.4.1 - (available at https://www.r-project.org/) and JAGS ≥ v4.3.2 - (available at https://mcmc-jags.sourceforge.io/).

Our pipeline uses the following packages for R:
- tidyverse v2.0.0
- r2jags v0.7.1.1
- mcmcplots v0.4.3
- progress v1.2.3
- cowplot v1.1.1
- patchwork v1.1.2
- gridExtra v2.3
- raster v3.6.26
- terra v1.7.78
- sf v1.0.17
- mapproj v1.2.11
- mapdata v2.3.1
- ggthemes v5.1.0
- rnaturalearth v1.0.1

## 2. Installation Guide

Follow the instructions for installing R and JAGS provided at the software download urls. 

The code and the mosquito life table experiment data are available here: https://github.com/mshocket/anopheles-rate-summation-release. You can clone the GitHub repository to your computer or download it as as zip file. 

World Clim 1.4 historical data for the mapping analyses are available here: https://worldclim.org/data/v1.4/worldclim14.html. Download the 5 minute resolution data for average temperature. 

## 3. Demo

We created a short demo that fits one thermal performance curve (TPC) to data from one life history trait (bite rate) from one fluctuation regime (lifespan data from constant temperature). It then applies rate summation to calculate a new TPC for two diurnal temperature ranges (DTR = 9 and 12 degrees Celcius). Finally, it plots the resulting TPCs as multiple curves on a single figure. 

To run the demo, open the file `R-scripts/AnalysisDemo.R` in R Studio or your R environment of choice and run the script.

You will need to set your working directory to the top-level anopheles-rate-summation-release directory.

## 4. Instructions for Use

TPC fitting and suitability modelling analysis pipeline

The top-level directory contains a R Project File `anopheles-rate-summation.Rproj`. If you are using R Studio as your R environment, begin by opening this file, and it will automatically set your working directory to that location. Otherwise you will need to set it manually at the beginning of your session.

The `R-Scripts` subdirectory contains a `working-versions-code` subdirectory with 7 R script files.

The files starting with 00 do not need to be run directly to reproduce the analysis.
- The file `00_JAGSModels.R` can be used to generate the required JAGS model text files in the R working directory (where they must be located for R/JAGS to find them); however, the text file versions of the models are already present in the  top level directory.
- The file `00_RSProjectFunctions.R` contains the functions that do most of the analysis - they are loaded and called in the other R scripts (see below).

The files starting with 01 - 05 must be run directly to generate the analysis. Open then in RStudio or your R environment of choice and run them in numerical sequence to reproduce the analysis. All of the necessary data are already included in the `data` and `processed-data` subdirectories.
- The file `01_TraitTPCFitting.R` fits and saves TPCS for each mosquito and malaria life history trait.
- The file `02_PlottingTraitTPCs.R` generates Figure 2 in the main text and SI Figure S1.
- The file `03_RSandSuitabilityCalcs.R` applies rate summation to the mosquito and malaria life history traits, calculates the thermal suitability for transmission, and generates Figures 3 and 4 in the main text.
- The file `04_SensitivityAnalysis.R` performs sensitivity and uncertainty analyses on the suitability models and generates SI Figures S2, S3, and S4. 

Mapping analysis pipeline

The mapping analysis is contained in a file called `05_RSMaps.R`. 

Before running it, you must first download the WorldClim data from the url provided in section 2 and save it in the `processed-data` subdirectory. (This file is not hosted on the GitHub repository because of its large size.)