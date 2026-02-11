library(lazyeval)
library(rstan) 
library(shinystan) 
library(data.table)
options(mc.cores=parallel::detectCores())
source('code/utils.R')
source("code/set_global_pars.R")

# the number of subjects
n_indiv <- length(unique(refined_data$PersonID))

# transform refined_data into qualified Stan-required dataframe
# Define a pared-down dataset for passing to Stan: 
trim_data <- refined_data %>% 
	select(PersonID, TestDate, 3) %>%
	rename(id=PersonID) %>% 
	rename(t=TestDate) %>%
	rename(y=3) %>%
  when(
    names(.)[3] == "CtValue" ~ trim_negatives(., global_pars),
    ~ .
  )

full_data <- refined_data %>% 
  select(PersonID, TestDate, 3) %>%
  rename(id=PersonID) %>% 
  rename(t=TestDate) %>%
  rename(y=3)


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