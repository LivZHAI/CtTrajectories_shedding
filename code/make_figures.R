# Sampled trajectories:
fig_ct_fit <- plot_ct_fit_SARS_trimdata(params_df, global_pars, trim_data, ctalpha=0.01, ntraces=100)

# Sampled trajectories:
fig_ct_fit_fulldata <- plot_ct_fit_SARS_fulldata(params_df, global_pars, full_data, ctalpha=0.01, ntraces=100)

meanvalsindiv <- params_df %>% 
	group_by(id) %>% 
	summarise(tp=mean(tp), dp=mean(dp), wp=mean(wp), wr=mean(wr))

jitterfactor <- 0.01
fig_dpmean <- shared_params_df %>% 
	select(dpmean) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=global_pars[["lod"]]-value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(10,NA)) + 
		scale_color_manual(values=c("blue")) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=0)) + 
		labs(x="Mean peak Ct", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_dpmean_withprior <- shared_params_df %>% 
	select(dpmean) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=global_pars[["lod"]]-value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(10,NA)) + 
		scale_color_manual(values=c("blue")) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=global_pars[["lod"]]-dp, y=0))  + 
		geom_line(data=make_normal_prior_df(mean=prior_pars$dpmean_prior, sd=prior_pars$dpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
		labs(x="Mean peak Ct", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

# transform ct value to viral load
fig_gemlmean <- shared_params_df %>% 
	select(dpmean) %>% 
	pivot_longer(everything()) %>% 
	mutate(value=10^convert_Ct_logGEML(global_pars[["lod"]]-value)) %>%
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x)), breaks=c(10^6, 10^7, 10^8, 10^9, 10^10)) + 
		# scale_x_continuous(limits=c(10,NA)) + 
		scale_color_manual(values=c("blue")) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=10^convert_Ct_logGEML(global_pars[["lod"]]-dp), y=0))  + 
		labs(x="Mean peak RNA copies per ml", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off


fig_wpmean <- shared_params_df %>% 
	select(wpmean) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("blue")) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=wp, y=0))  + 
		labs(x="Mean proliferation stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_wpmean_withprior <- shared_params_df %>% 
	select( wpmean) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("blue")) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=wp, y=0))  + 
		geom_line(data=make_normal_prior_df(mean=prior_pars$wpmean_prior, sd=prior_pars$wpsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
		labs(x="Mean proliferation stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off


fig_wrmean <- shared_params_df %>% 
	select(wrmean) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("blue")) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=wr, y=0))  + 
		labs(x="Mean clearance stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off


fig_wrmean_withprior <- shared_params_df %>% 
	select(wrmean) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("blue")) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=wr, y=0))  + 
		geom_line(data=make_normal_prior_df(mean=prior_pars$wrmean_prior, sd=prior_pars$wrsd_prior, pmin=0.001, pmax=0.999, step=0.1), aes(x=x, y=density), col="black", linetype="dashed") + 
		labs(x="Mean clearance stage duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_infdurmean <- shared_params_df %>% 
	mutate(infdurmean=wpmean+wrmean) %>%
	select(infdurmean) %>% 
	pivot_longer(everything()) %>% 
	ggplot(aes(x=value)) + 
		# geom_histogram(aes(y=..density.., fill=name), alpha=0.2, position="identity", bins=50) + 
		geom_density(aes(col=name, fill=name), alpha=0.2, adjust=2) + 
		scale_x_continuous(limits=c(0,NA)) + 
		scale_color_manual(values=c("blue"),) + 
		scale_fill_manual(values=c("blue")) + 
		geom_point(data=meanvalsindiv, aes(x=wp+wr, y=0))  + 
		labs(x="Mean acute infection duration (days)", y="Density") + 
		theme_minimal() + 
		# theme(legend.title=element_blank(), text=element_text(size=18)) + 
		theme(legend.position="none", text=element_text(size=18)) + 
		y_ticks_off + 
		grid_off

fig_ct_trajectory_inference <- make_sample_trajectory(shared_params_df, global_pars)
fig_ct_trajectory_inference_geml <- make_sample_trajectory(shared_params_df, global_pars, ge=TRUE)

# piecewise function of viral trajectory
fig_raw_trajectories <- with(as.list(global_pars),{
meanvalsindiv %>% 
	ggplot() + 
		geom_segment(aes(x=-Inf, y=lod, xend=-wp, yend=lod), alpha=0.5, color = "blue") + 
		geom_segment(aes(x=-wp, y=lod, xend=0, yend=lod-dp), alpha=0.5, color = "blue") + 
		geom_segment(aes(x=0, y=lod-dp, xend=wr, yend=lod), alpha=0.5, color = "blue") + 
		geom_segment(aes(x=wr, y=lod, xend=Inf, yend=lod), alpha=0.5, color = "blue") + 
		theme_minimal() + 
		theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), text=element_text(size=18)) + 
		labs(x="Time since min Ct (days)", y="Ct") + 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml"))
})


