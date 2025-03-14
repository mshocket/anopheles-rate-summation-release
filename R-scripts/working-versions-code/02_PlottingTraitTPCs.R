## Rate Summation Project
## Script for plotting trait data and TPC fits
## Written by Marta Shocket & Joey Bernhardt in 2022-2024

##	Table of contents:
##
##	1. Set-up workspace
##	2. Load and process data for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##	3. Plot panels for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##	4. Manuscript Figure 2
##	5. Load, process, and plot data for other traits (bc, EIP50, pEA, MDR, gamma)
##	6. Manuscript Figure S1


##########
###### 1. Set-up workspace
##########

### Load libraries
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)

##### Load functions
source("R-scripts/working-versions-code/00_RSProjectFunctions.R")

##### Make figures output folder
dir.create("figures", showWarnings = FALSE)


##########
###### 2. Load and process data for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##########

####################### Bite rate (a)

### load data - all fluctuation treatments
predictions_bite_rate_constant_summary <- read.csv("data-processed/predictions_bite_rate_constant_summary.csv")
params_bite_rate_constant_summary <- read.csv("data-processed/params_bite_rate_constant_summary.csv")
data_bite_rate_constant_summary <- read.csv("data-processed/data_bite_rate_constant_summary.csv")

predictions_bite_rate_dtr9_summary <- read.csv("data-processed/predictions_bite_rate_dtr9_summary.csv")
params_bite_rate_dtr9_summary <- read.csv("data-processed/params_bite_rate_dtr9_summary.csv")
data_bite_rate_dtr9_summary <- read.csv("data-processed/data_bite_rate_dtr9_summary.csv")

predictions_bite_rate_dtr12_summary <- read.csv("data-processed/predictions_bite_rate_dtr12_summary.csv")
params_bite_rate_dtr12_summary <- read.csv("data-processed/params_bite_rate_dtr12_summary.csv")
data_bite_rate_dtr12_summary <- read.csv("data-processed/data_bite_rate_dtr12_summary.csv")

### edit DTR9 treatment labels so they're alphabetically before DTR12 - ggplot2 assigns categorical variables colors in alphabetical order
predictions_bite_rate_dtr9_summary$treatment <- "bite_rate_dtr09"
params_bite_rate_dtr9_summary$treatment <- "bite_rate_dtr09"
data_bite_rate_dtr9_summary$treatment <- "bite_rate_dtr09"

### combine all fluctuation treatments
all_bite_rate_predictions <- bind_rows(predictions_bite_rate_constant_summary, predictions_bite_rate_dtr9_summary, predictions_bite_rate_dtr12_summary)
all_bite_rate_data <- bind_rows(data_bite_rate_constant_summary, data_bite_rate_dtr9_summary, data_bite_rate_dtr12_summary)

### combine and merge parameters with different names 
all_bite_rate_params <- bind_rows(params_bite_rate_constant_summary, params_bite_rate_dtr9_summary, params_bite_rate_dtr12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

####################### Lifespan (lf)

### load data - all fluctuation treatments
predictions_lifespan_constant_summary <- read.csv("data-processed/predictions_lifespan_constant_summary.csv")
params_lifespan_constant_summary <- read.csv("data-processed/params_lifespan_constant_summary.csv")
data_lifespan_constant_summary <- read.csv("data-processed/data_lifespan_constant_summary.csv")

predictions_lifespan_dtr9_summary <- read.csv("data-processed/predictions_lifespan_dtr9_summary.csv")
params_lifespan_dtr9_summary <- read.csv("data-processed/params_lifespan_dtr9_summary.csv")
data_lifespan_dtr9_summary <- read.csv("data-processed/data_lifespan_dtr9_summary.csv")

predictions_lifespan_dtr12_summary <- read.csv("data-processed/predictions_lifespan_dtr12_summary.csv")
params_lifespan_dtr12_summary <- read.csv("data-processed/params_lifespan_dtr12_summary.csv")
data_lifespan_dtr12_summary <- read.csv("data-processed/data_lifespan_dtr12_summary.csv")

### edit DTR9 treatment labels so they're alphabetically before DTR12 - ggplot2 assigns categorical variables colors in alphabetical order
predictions_lifespan_dtr9_summary$treatment <- "lifespan_dtr09"
params_lifespan_dtr9_summary$treatment <- "lifespan_dtr09"
data_lifespan_dtr9_summary$treatment <- "lifespan_dtr09"

### combine all fluctuation treatments
all_lifespan_predictions <- bind_rows(predictions_lifespan_constant_summary, predictions_lifespan_dtr9_summary, predictions_lifespan_dtr12_summary)
all_lifespan_data <- bind_rows(data_lifespan_constant_summary, data_lifespan_dtr9_summary, data_lifespan_dtr12_summary)

### combine and merge different names
all_lifespan_params <- bind_rows(params_lifespan_constant_summary, params_lifespan_dtr9_summary, params_lifespan_dtr12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

####################### Lifetime eggs (B)

### load data - all fluctuation treatments
predictions_eggs_constant_summary <- read.csv("data-processed/predictions_eggs_constant_summary.csv")
params_eggs_constant_summary <- read.csv("data-processed/params_eggs_constant_summary.csv")
data_eggs_constant_summary <- read.csv("data-processed/data_eggs_constant_summary.csv")

predictions_eggs_dtr9_summary <- read.csv("data-processed/predictions_eggs_dtr9_summary.csv")
params_eggs_dtr9_summary <- read.csv("data-processed/params_eggs_dtr9_summary.csv")
data_eggs_dtr9_summary <- read.csv("data-processed/data_eggs_dtr9_summary.csv")

predictions_eggs_dtr12_summary <- read.csv("data-processed/predictions_eggs_dtr12_summary.csv")
params_eggs_dtr12_summary <- read.csv("data-processed/params_eggs_dtr12_summary.csv")
data_eggs_dtr12_summary <- read.csv("data-processed/data_eggs_dtr12_summary.csv")

### edit DTR9 treatment labels so they're alphabetically before DTR12 - ggplot2 assigns categorical variables colors in alphabetical order
predictions_eggs_dtr9_summary$treatment <- "eggs_dtr09"
params_eggs_dtr9_summary$treatment <- "eggs_dtr09"
data_eggs_dtr9_summary$treatment <- "eggs_dtr09"

#### combine all fluctuation treatments
all_eggs_predictions <- bind_rows(predictions_eggs_constant_summary, predictions_eggs_dtr9_summary, predictions_eggs_dtr12_summary)
all_eggs_data <- bind_rows(data_eggs_constant_summary, data_eggs_dtr9_summary, data_eggs_dtr12_summary)

### combine and merge different names 
all_eggs_params <- bind_rows(params_eggs_constant_summary, params_eggs_dtr9_summary, params_eggs_dtr12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))


##########
###### 3. Plot panels for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##########

### plot for bite rate - xlim() gives warning because not all rows are used in the plot
bite_rate_plot <- all_bite_rate_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_bite_rate_data, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_bite_rate_data, size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "Empirical DTR 9", "Empirical DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("Constant", "Empirical DTR 9", "Empirical DTR 12")) +
	ylab(parse(text = "Bite~rate~(day^-1)~-~a")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = c(0.2, 0.8), legend.title=element_blank(), legend.text = element_text(size=11)) + 
	annotate("text", x = 0, y = 0.53, label = "A", size = 5)

params_plot_bite_rate <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_bite_rate_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "D", size = 5)

bite_rate_plot_all <- bite_rate_plot / params_plot_bite_rate + plot_layout(heights = c(3, 0.5))

### plot for lifespan
lifespan_plot <- all_lifespan_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_lifespan_data, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_lifespan_data, size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	ylab("Lifespan (days) - lf") + xlab("Temperature (°C)") + xlim(0,40) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 51, label = "B", size = 5)

params_plot_lifespan <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_lifespan_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) + ylim(0,40) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "E", size = 5)

lifespan_plot_all <- lifespan_plot / params_plot_lifespan + plot_layout(heights = c(3, 0.5))

### plot for eggs
eggs_plot <- all_eggs_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_eggs_data, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_eggs_data, size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	ylab("Lifetime eggs - B") + xlab("Temperature (°C)") +xlim(0, 40) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 450, label = "C", size = 5)

params_plot_eggs <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_eggs_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "F", size = 5)

eggs_plot_all <- eggs_plot / params_plot_eggs + plot_layout(heights = c(3, 0.5))


##########
###### 4. Manuscript Figure 2
##########

Fig2_empTPCs <- wrap_plots(bite_rate_plot_all, lifespan_plot_all, eggs_plot_all)
ggsave('figures/Fig2_trait_TPCs_JAGS.pdf', Fig2_empTPCs, width = 15, height = 6)


##########
###### 5. Load and process data for other traits (bc, EIP50, pEA, MDR, gamma)
##########

####################### Gamma (proportion surviving EIP)

### load data - separate fluctuation treatments
predictions_gamma_constant_summary <- read.csv("data-processed/predictions_gamma_constant_summary.csv")
params_gamma_constant_summary <- read.csv("data-processed/params_gamma_constant_summary.csv")
data_gamma_constant_summary <- read.csv("data-processed/data_gamma_constant_summary.csv")

predictions_gamma_dtr9_summary <- read.csv("data-processed/predictions_gamma_dtr9_summary.csv")
params_gamma_dtr9_summary <- read.csv("data-processed/params_gamma_dtr9_summary.csv")
data_gamma_dtr9_summary <- read.csv("data-processed/data_gamma_dtr9_summary.csv")

predictions_gamma_dtr12_summary <- read.csv("data-processed/predictions_gamma_dtr12_summary.csv")
params_gamma_dtr12_summary <- read.csv("data-processed/params_gamma_dtr12_summary.csv")
data_gamma_dtr12_summary <- read.csv("data-processed/data_gamma_dtr12_summary.csv")

### edit DTR9 treatment labels so they're alphabetically before DTR12 - ggplot2 assigns categorical variables colors in alphabetical order
predictions_gamma_dtr9_summary$treatment <- "bite_rate_dtr09"
params_gamma_dtr9_summary$treatment <- "bite_rate_dtr09"
data_gamma_dtr9_summary$treatment <- "bite_rate_dtr09"

#### combine all fluctuation treatments
all_gamma_predictions <- bind_rows(predictions_gamma_constant_summary, predictions_gamma_dtr9_summary, predictions_gamma_dtr12_summary)
all_gamma_data <- bind_rows(data_gamma_constant_summary, data_gamma_dtr9_summary, data_gamma_dtr12_summary)

### plot for gamma
gamma_plot <- all_gamma_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_gamma_data, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_gamma_data, size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "Empirical DTR 9", "Empirical DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("Constant", "Empirical DTR 9", "Empirical DTR 12")) +
	ylab("Gamma (proportion)") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = c(0.15, 0.8), legend.title=element_blank(), legend.text = element_text(size=11)) + 
	annotate("text", x = 0, y = 1, label = "A", size = 5)

params_plot_gamma <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_lifespan_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) + ylim(0,40) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "F", size = 5)

gamma_plot_all <- gamma_plot / params_plot_gamma + plot_layout(heights = c(3, 0.5))


####################### Vector competence (bc)

### load data
predictions_bc_summary <- read.csv("data-processed/predictions_bc_summary.csv")
params_bc_summary <- read.csv("data-processed/params_bc_summary.csv")
data_bc_summary <- read.csv("data-processed/data_bc_summary.csv")

### plot for bc
bc_plot <- predictions_bc_summary %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = data_bc_summary, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = data_bc_summary, size = 0.8) +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	ylab("Vector competence (proportion) - bc") + xlab("Temperature (°C)") + xlim(0, 45) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 0.65, label = "B", size = 5)

### filter and rename parameters
all_params_bc <- params_bc_summary %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### parameter plot for bc
params_plot_bc <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_bc, position=position_dodge(width=1)) +
	ylim(0, 45) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "G", size = 5)

bc_plot_all <- bc_plot / params_plot_bc + plot_layout(heights = c(3, 0.5))


####################### Extrinsic Incubation Period 50% (EIP50) / Pathogen development rate (PDR)

### load data
predictions_eip50_summary <- read.csv("data-processed/predictions_eip50_summary.csv")
params_eip50_summary <- read.csv("data-processed/params_eip50_summary.csv")
data_eip50_summary <- read.csv("data-processed/data_eip50_summary.csv")

### plot for eip50/PDR
eip50_plot <- predictions_eip50_summary %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = data_eip50_summary, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = data_eip50_summary, size = 0.8) +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	ylab(parse(text = "Pathogen~dev.~rate~(day^-1)~-~PDR")) + xlab("Temperature (°C)") + xlim(0, 45) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 0.16, label = "C", size = 5)

### filter and rename parameters
all_params_eip50 <- params_eip50_summary %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### parameter plot for eip50/PDR
params_plot_eip50 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_eip50, position=position_dodge(width=1)) +
	ylim(0, 45) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "H", size = 5)

eip50_plot_all <- eip50_plot / params_plot_eip50 + plot_layout(heights = c(3, 0.5))


####################### Juvenile survival (pEA)

### Load data
predictions_pea_summary <- read.csv("data-processed/predictions_pea_summary.csv")
params_pea_summary <- read.csv("data-processed/params_pea_summary.csv")
data_pea_summary <- read.csv("data-processed/data_pea_summary.csv")

### plot for pea
pea_plot <- predictions_pea_summary %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = data_pea_summary, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = data_pea_summary, size = 0.8) +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	ylab("Juvenile survival - pEA") + xlab("Temperature (°C)") + xlim(0, 45) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 1, label = "D", size = 5)

### filter and rename parameters
all_params_pea <- params_pea_summary %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### parameter plot for pEA
params_plot_pea <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_pea, position=position_dodge(width=1)) +
	ylim(0, 45) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "I", size = 5)

pea_plot_all <- pea_plot / params_plot_pea + plot_layout(heights = c(3, 0.5))


####################### Mosquito development rate (MDR)

### load data
predictions_mdr_summary <- read.csv("data-processed/predictions_mdr_summary.csv")
params_mdr_summary <- read.csv("data-processed/params_mdr_summary.csv")
data_mdr_summary <- read.csv("data-processed/data_mdr_summary.csv")

### plot for mdr
mdr_plot <- predictions_mdr_summary %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = data_mdr_summary, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = data_mdr_summary, size = 0.8) +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	ylab(parse(text = "Mosquito~development~rate~(day^-1)~-~MDR")) + xlab("Temperature (°C)") + xlim(0, 45) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 0.16, label = "E", size = 5)

### filter and rename parameters
all_params_mdr <-params_mdr_summary %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### parameter plot for MDR
params_plot_mdr <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_mdr, position=position_dodge(width=1)) +
	ylim(0, 45) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c_constant) +
	scale_fill_manual(values = ct_constant) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "J", size = 5)

mdr_plot_all <- mdr_plot / params_plot_mdr + plot_layout(heights = c(3, 0.5))


##########
######6. Manuscript Figure S1 
##########

other_trait_plots <- wrap_plots(gamma_plot_all, bc_plot_all, eip50_plot_all, pea_plot_all, mdr_plot_all, 
							nrow = 3, ncol = 2)
ggsave('figures/FigS1_other_traits.pdf', other_trait_plots, width = 13, height = 16)
ggsave('figures/FigS1_other_traits.png', other_trait_plots, width = 13, height = 16)

