scenario_df <- data.frame(scenario=c("one_winner","null"),
n_per_arm = rep(125,2),control_surv_rate = rep(0.7,2), hr1=c(0.7,1), hr2=c(1.3,1), 
  admin_cens_time = rep(14,2), drop_out=rep(0.1,2), stringsAsFactors = FALSE)

write.csv(scenario_df, file="scenario_df.csv",row.names = FALSE)