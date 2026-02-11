# calculate quantile of posterior estimates
# output includes three cols，which are mean，90low，90up
dist_summary <- shared_params_df %>% 
  summarise(
    # Peak Ct values
    peak.ct.mean=mean(global_pars[["lod"]]-dpmean),
    peak.ct.lwr95=quantile(global_pars[["lod"]]-dpmean,0.05),
    peak.ct.upr95=quantile(global_pars[["lod"]]-dpmean,0.95),
    
    # Peak GEML (log RNA copies/ml)
    peak.geml.mean=mean(convert_Ct_logGEML(global_pars[["lod"]]-dpmean)),
    peak.geml.lwr95=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmean),0.05),
    peak.geml.upr95=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmean),0.95),
    
    # Proliferation time
    proliferation.time.mean=mean(wpmean),
    proliferation.time.lwr95=quantile(wpmean,0.05),
    proliferation.time.upr95=quantile(wpmean,0.95),
    
    # Clearance time
    clearance.time.mean=mean(wrmean),
    clearance.time.lwr95=quantile(wrmean,0.05),
    clearance.time.upr95=quantile(wrmean,0.95),
    
    # Total duration
    total.duration.mean=mean(wpmean+wrmean),
    total.duration.lwr95=quantile(wpmean+wrmean,0.05),
    total.duration.upr95=quantile(wpmean+wrmean,0.95),
  ) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("parameter", "statistic"), sep = "\\.(?=[^\\.]+$)") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  arrange(parameter)
