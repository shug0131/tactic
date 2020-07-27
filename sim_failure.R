#function to provide survival data
# goes up to 14 days - discretised to integer days.
# set HR in multiple arms
# final response rate in control arm. 

## @knitr sim_failure


sim_failure <- function(n_per_arm, control_surv_rate, hr, site_raneff_sd=0.1, last_time=14){
  require(magrittr)
  #set parameters
  n_arms <- length(hr)+1
  rx <- rep(1:n_arms, rep(n_per_arm,n_arms)) %>% 
    factor(labels=c("Cx",paste0("E",1:(n_arms-1))))
  #lazy replication
  site <- rep(1, length(rx))
  site[1:length(site)] <- 1:10
  site_re <- rnorm(n=10,mean=0, sd=site_raneff_sd)
  #assume exponential distribution
  hazard_cx <-   -log(control_surv_rate)/last_time
  hazard <- hazard_cx*c(1,hr)
  time_exact <- rexp(length(rx),rate=hazard[rx %>% as.numeric]*exp(site_re[site]))
  time <- pmin(ceiling(time_exact), last_time)
  status <- 1*(time_exact<last_time)
  data.frame(rx, time, status, time_exact,site)
}

add_censor <- function(data, admin_cens_time, drop_out, drop_out_time=max(admin_cens_time)){
  haz <-  -log(1-drop_out)/drop_out_time
  n <- nrow(data)
  cens <- rexp(n, rate=haz)
  cens <- pmin(cens, admin_cens_time)
  time_exact <- pmin(data$time_exact, cens)
  status <- ifelse( cens<data$time_exact, 0, data$status)
  time <-  ceiling(time_exact)
  data$time <- time
  data$status <- status
  data$cens <- cens  
  data
}

#test
if(FALSE){
  df <- sim_failure(n_per_arm = 1000,control_surv_rate = 0.7, hr=c(0.7,1.3)) %>% 
    add_censor(admin_cens_time = 14, drop_out=0.1)
  survfit(Surv(time, status)~rx,data=df) %>% plot(mark.time=TRUE,lty=1:3)
  coxph(Surv(time, status)~rx,data=df) %>% summary
#check the censoring.
  survfit(Surv(time,1-status)~1, data=df) %>% plot
  abline(h=0.9)
}