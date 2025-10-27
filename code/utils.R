library(tidyverse) 
library(ggplot2)
library(yaml)
library(here)


# import .yaml file
raw_data <- yaml::read_yaml("data/team2020clinical.yaml")

# transform .yaml file to dataframe
participant_data <- map_df(seq_along(raw_data$participants), function(i) {
  participant <- raw_data$participants[[i]]
  measurements <- map_df(participant$measurements, function(meas) {
    ct_value <- ifelse(meas$value == "negative", 40, as.numeric(meas$value))
    data.frame(
      time = meas$time,
      value = ct_value,
      stringsAsFactors = FALSE
    )
  })
  measurements |>
    mutate(
      participant_id = i
    ) |>
    select(participant_id, everything())
})

# standardize dataframe form (recalculate time, change column name)
refined_data <- participant_data |>
  rename(PersonIDClean = participant_id,
         TestDateIndex = time,
         CtT1 = value
         ) |>
  group_by(PersonIDClean) |>
  # create new testDate based on min ct-value date
  mutate(TestDateIndex = TestDateIndex - TestDateIndex[which.min(CtT1)]) |>
  ungroup()
  

################################
# trim redundant negative ct results
# keep the point if any item is satisfied：
# ispositive==1: current point is pos
# ispositive_lag==1: the previous one point is pos
# ispositive_lag2==1: the previous two points are pos
# ispositive_lead==1: the next one point is pos
# ispositive_lead2==1: the next two points are pos
trim_negatives <- function(full_data, global_pars){
	# lod <- 40
	out <- full_data %>% 
		split(.$id) %>% 
		map(~ arrange(., t)) %>% 
		map(~ mutate(., rowindex=1:n())) %>% 
		map(~ mutate(., ispositive=case_when(y<global_pars[["lod"]] ~ 1, TRUE~0))) %>% 
		map(~ mutate(., ispositive_lag=lag(ispositive))) %>%
		map(~ mutate(., ispositive_lag2=lag(ispositive,2))) %>%
		map(~ mutate(., ispositive_lead=lead(ispositive))) %>%
		map(~ mutate(., ispositive_lead2=lead(ispositive,2))) %>%
		map(~ filter(., ispositive==1 | ispositive_lag==1 | ispositive_lag2==1 | ispositive_lead==1 | ispositive_lead2==1)) %>%
		map(~ select(., id, t, y)) %>%
		bind_rows()
	return(out)
}

launch_shinystan_nonblocking <- function(fit) {
  library(future)
  plan(multisession)
  future(
    launch_shinystan(fit) 
  )
}

# randomly select samples from grouped data
sample_n_groups <- function(grouped_df, size, replace=FALSE, weight=NULL){
	# From https://cmdlinetips.com/2019/07/how-to-randomly-select-groups-in-r-with-dplyr/
	grp_var <- grouped_df %>% 
	    groups %>%
	    unlist %>% 
	    as.character

	random_grp <- grouped_df %>% 
	    summarise() %>% 
	    slice_sample(n=size, replace=replace, weight_by=weight)# %>% 
	    # mutate(unique_id = 1:n())

	out <- grouped_df %>% 
	    right_join(random_grp, by=grp_var)

	return(out)
}

################################  
# functions for extracted parameters from Stan
# For generating names for a matrix-turned-data frame: 
makenames <- function(parname, n_indiv){
	unlist(lapply(1:n_indiv, function(x) paste0(parname, "_", x)))
}

# For parsing the 'parameters' output from Stan: 
# the following functions are used in fit_posteriors.R
parseparam <- function(extracted_params, parname, n_indiv){
	as_tibble(setNames(as.data.frame(extracted_params[[parname]]), makenames(parname,n_indiv)))
}

make_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value") %>%
		select(-iteration) 
}

make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value")
}

make_shared_params_df <- function(extracted_params, parnames){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) 
		as_tibble(setNames(as.data.frame(extracted_params[[x]]),x))
		), cbind) %>%
		as_tibble() %>%
		mutate(iteration=1:n())
}

################################   
# indiv_data applies function of 'trim_negatves', to ignore meaningless negative ct values

# indiv_data is defined in fit_posteriors_preamble.R
# 100 trajectories for each person shown in separate plots
plot_ct_fit_SARS_indivdata <- function(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=100){
  with(as.list(global_pars),{
    params_df %>% 
      group_by(id) %>%
      sample_n(ntraces) %>% 
      ungroup() %>%
      ggplot() + 
      # Plot traces:
      geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
      geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
      geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
      geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
      # Plot data:
      geom_point(data=indiv_data, aes(x=t, y=y), size=0.5, color="blue") + 
      theme_minimal() + 
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
      labs(x="Time since min Ct (days)", y="Ct") + 
      scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
      facet_wrap(~id)
  })
}

plot_ct_fit_SARS_fulldata <- function(params_df, global_pars, ct_dat_refined, ctalpha=0.01, ntraces=100){
	with(as.list(global_pars),{
	params_df %>% 
		group_by(id) %>%
		sample_n(ntraces) %>% 
		ungroup() %>%
		ggplot() + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
	    geom_point(data=indiv_data, aes(x=t, y=y), size=0.5, color="blue") + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
			labs(x="Time since min Ct (days)", y="Ct") + 
			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
			facet_wrap(~id)
			})
}

grid_off <- list(theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))

y_ticks_off <- list(theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()))

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
	out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
	return(out) 
}

srise <- function(x, dp, wp){
	out <- dp*(1+x/wp)
	return(out)
}
sfall <- function(x, dp, wr){
	out <- dp*(1-x/wr)
	return(out)
}

################################
# get estimate of posterior parameters by calculating mean; generate ct trajectory by using estimates.
make_sample_trajectory <- function(shared_params_df, global_pars, siglevel=0.9, ge=FALSE){
	# For asymptomatic:
	with(as.list(global_pars),{

	wp_mean_W <- mean(shared_params_df$wpmeanW)
	wp_lwr_W <- quantile(shared_params_df$wpmeanW,(1-siglevel)/2)
	wp_upr_W <- quantile(shared_params_df$wpmeanW,1-(1-siglevel)/2)

	wr_mean_W <- mean(shared_params_df$wrmeanW)
	wr_lwr_W <- quantile(shared_params_df$wrmeanW,(1-siglevel)/2)
	wr_upr_W <- quantile(shared_params_df$wrmeanW,1-(1-siglevel)/2)

	dp_mean_W <- mean(shared_params_df$dpmeanW)
	dp_lwr_W <- quantile(shared_params_df$dpmeanW,(1-siglevel)/2)
	dp_upr_W <- quantile(shared_params_df$dpmeanW,1-(1-siglevel)/2)

	xvals_proliferation_W <- seq(from=-wp_upr_W, 0, length.out=500)
	xvals_clearance_W <- seq(from=0, wr_upr_W, length.out=500)

	yvals_upr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),(1-siglevel)/2)))
	yvals_upr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),1-(1-siglevel)/2)))
	yvals_lwr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),(1-siglevel)/2)))


	if(ge==FALSE){
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=yvals_lwr_proliferation_W, 
					yvals_upr=yvals_upr_proliferation_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="blue") + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=lod,yend=lod-dp_mean_W),col="blue") + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=yvals_lwr_clearance_W, 
					yvals_upr=yvals_upr_clearance_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="blue") + 
			geom_segment(aes(x=0,xend=wr_mean_W,y=lod-dp_mean_W,yend=lod),col="blue") + 
			coord_cartesian(ylim=c(40,15), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Ct") + 
			scale_y_reverse() + 
			theme(text=element_text(size=18))
	} else {
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=(yvals_lwr_proliferation_W), 
					yvals_upr=(yvals_upr_proliferation_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="blue") + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_W)),col="blue") + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=(yvals_lwr_clearance_W), 
					yvals_upr=(yvals_upr_clearance_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="blue") + 
			geom_segment(aes(x=0,xend=wr_mean_W,y=10^convert_Ct_logGEML(lod-dp_mean_W),yend=10^convert_Ct_logGEML(lod)),col="blue") + 
			coord_cartesian(ylim=c(10^convert_Ct_logGEML(40),10^convert_Ct_logGEML(15)), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="RNA copies per ml") + 
			scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
			theme(text=element_text(size=18))
	}
	return(out)
	})
}

################################
# Generate a data frame that contains the pdf of the prior normal distribution within the specified range to 
# facilitate the visualization of this prior distribution.

make_normal_prior_df <- function(mean, sd, pmin, pmax, step){
	xmin <- qnorm(pmin, mean=mean, sd=sd)
	xmax <- qnorm(pmax, mean=mean, sd=sd)
	xvals <- seq(from=xmin, to=xmax, by=step)
	out <- as_tibble(data.frame(
		x=xvals, 
		density=dnorm(xvals, mean=mean, sd=sd)))
	return(out)
}
