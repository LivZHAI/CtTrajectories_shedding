library(tidyverse)
library(lazyeval)
library(rstan) 
library(shinystan) 
library(purrr)
library(data.table)
options(mc.cores=parallel::detectCores())
source('code/utils.R')
source("code/set_global_pars.R")

# Store the number of people we've kept: 
n_indiv <- length(unique(refined_data$PersonIDClean))

# transform refined_data into qualified Stan-required dataframe
# Define a pared-down dataset for passing to Stan: 
indiv_data <- refined_data %>% 
	select(PersonIDClean, TestDateIndex, CtT1) %>%
	rename(id=PersonIDClean) %>% 
	rename(t=TestDateIndex) %>%
	rename(y=CtT1) %>%
  trim_negatives(global_pars)

# Trim to 6 b117 and 6 non-b117 (comment for a full run):
# indiv_data <- indiv_data %>% 
# 	split(.$b117) %>%
# 	map(~ group_by(., id)) %>% 
# 	map(~ sample_n_groups(., 10)) %>% 
# 	bind_rows() %>%
# 	rename(PersonID=id) %>%
# 	select(-id_clean) %>%
# 	clean_person_id %>%
# 	rename(id=PersonID) %>% 
# 	rename(id_clean=PersonIDClean) %>%
# 	select(id, id_clean, t, y, b117) %>% 
# 	ungroup() 
# n_indiv <- length(unique(indiv_data$id))


# Useful dataframe for mapping official ID to Stan ID:
#id_map <- indiv_data %>% 
#	group_by(id) %>%
#	summarise(id_clean=first(id_clean)) %>% 
#	select(id, id_clean) %>%
#	mutate(id_clean=as.character(id_clean))

# Useful dataframe for mapping official ID to symptoms
#b117_map <- indiv_data %>% 
#	group_by(id) %>%
#	summarise(b117=first(b117))

prior_pars <- list(
	tpsd=2,
	dpmean_prior=global_pars[["lod"]]/2,
	dpsd_prior=global_pars[["lod"]]/6,
	wpmax=14,
	wpmean_prior=14/2,
	wpsd_prior=14/6,
	wrmax=30,
	wrmean_prior=30/2,
	wrsd_prior=30/6,
	sigma_max=10,
	sigma_prior_scale=5,
	lambda=0.01,
	fpmean=1/log(10)  # so that 90% of mass is <1 and 99% is <2
	)