# read fit_posteriors.stan
ct_model <- stan_model("code/fit_posteriors.stan") 

fit_startq <- Sys.time()

set.seed(2)
ct_fit <- sampling(ct_model, 
	data=list(
		N=nrow(trim_data), 
		n_id=length(unique(trim_data$id)),
		lod=global_pars[["lod"]], 
		# add value_type
		value_type = ifelse((names(refined_data)[3] == 'CtValue'), 1, 0),
		id=trim_data$id,
		t=trim_data$t, 
		y=trim_data$y, 
		tpsd=as.list(prior_pars)$tpsd,
		dpmean_prior=as.list(prior_pars)$dpmean_prior,
		dpsd_prior=as.list(prior_pars)$dpsd_prior,
		wpmax=as.list(prior_pars)$wpmax,
		wpmean_prior=as.list(prior_pars)$wpmean_prior,
		wpsd_prior=as.list(prior_pars)$wpsd_prior,
		wrmax=as.list(prior_pars)$wrmax,
		wrmean_prior=as.list(prior_pars)$wrmean_prior,
		wrsd_prior=as.list(prior_pars)$wrsd_prior,
		sigma_max=as.list(prior_pars)$sigma_max,
		sigma_prior_scale=as.list(prior_pars)$sigma_prior_scale,
		lambda=as.list(prior_pars)$lambda,
		fpmean=as.list(prior_pars)$fpmean), 
	# random sampling parameters
	  control = list(adapt_delta=0.80,stepsize = 0.05, max_treedepth=20), 
	iter=2000, chains=4)

fit_endq <- Sys.time()
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))

# launch_shinystan_nonblocking(ct_fit)

params <- rstan::extract(ct_fit)
indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wr"), n_indiv) 
shared_params_df <- make_shared_params_df(params, c("dpmean","wpmean","wrmean","dpsd","wpsd","wrsd")) 

params_df <- indiv_params_df %>% 
	left_join(shared_params_df, by="iteration") %>% 
	select(-iteration) 

# write.csv(shared_params_df,file='output/shared_params_df.csv',row.names=FALSE)
# write.csv(indiv_params_df,file='output/indiv_params_df.csv',row.names=FALSE)