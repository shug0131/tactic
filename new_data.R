## @knitr new_data


#Generate data for posterior_predict(newdata)

#Bootstrap from the original data the desired nunmber of times

#n_new is defined globally
if(!exists("n_new")){warning("setting n_new to 200"); n_new <- 200}

df_new <- df[sample(1:nrow(df), n_new, replace=TRUE),]
#summary(df_new)
#force the randomisation argument to be an equal split seperate from teh original data
n <- nrow(df_new)
n_arms <- nlevels(df_new$rx)
n_per_arm <- n%/%n_arms
spare <- n%%n_arms
df_new$rx <- factor(c(rep(1:n_arms, n_per_arm),rep(1,spare)), labels=levels(df$rx))
df_new$id <- 1:nrow(df_new)
#convert to long format with a row for each day
long <- expand.grid(time=sort(fail_times), id=df_new$id)
df_new_long <- left_join(long, df_new %>% select(-time), by="id")


# FUnction to add back in a simulated status endppoint
# and convert to standard survival data : either long counting, or short time/status

add_sim_status <- function(df, posterior,...){
  # ... give the covariates to keep.
  ##  future_analysis doesn't pass this on directly at the moment
  
  require(dplyr)
  require(magrittr)
  # this works one row at a time for posterior
  df$status <- posterior
 df_fail_time <- df %>%  group_by(id) %>% 
    filter(status>0) %>% summarise( fail_time=min(time), .groups="drop")
 # Say you had a  discrete survival distribution to simulate from
 # geenrate a U(0,1) and then invert the CDF
 # grouping together the poisson outcomes >=1, does this in effect
 # a value >1 means the precise failure occured earlier, but within that interval.
 
  
  df_long<- left_join(df, df_fail_time, by=c("id"="id")) %>% 
    mutate( fail_time=ifelse(is.na(fail_time),Inf, fail_time)) %>% 
    filter( time<=fail_time)
  
 
    covars <- c(quos(id,rx), quos(...))
  
  df_short <- df_long %>% group_by(!!!covars) %>% 
    summarise(time=max(time), status=min(1,max(status)), .groups="drop")
  list("long"=df_long, "short"=df_short)
}
#test
if(FALSE){
  ans <- add_sim_status(df_new_long , pred_data[6,],"rx")
  plot(survfit(Surv(time, status)~rx, data=ans$short))
  plot(survfit(Surv(time, status)~rx, data=df))
  df %>% group_by(rx) %>% summarise(event=mean(status))
  ans$short %>% group_by(rx) %>% summarise(event=mean(status))
  mean(df$status)
  mean(ans$short$status)
}




### Now add back in censoring ...


# Generate a censoring distribution survfit(Surv(time, 1- status)~1).  
censor_dist <- survfit(Surv(time, 1-status)~1, data=df) 
time <- censor_dist$time
S <- censor_dist$surv
p <- -diff(c(1,S))
index <- which(p!=0)
censor_dist <- data.frame(time=time[index], prob=p[index])
# But add in an admin censor at the end if needed - might vary between patients if we have calendar censoring.
admin_cens <- 15
# make this the last day+1
# the rule will be that censoring wins for ties
# so this allows people to fail on day 14, rather than being censored - if they get day 15
early <- subset(censor_dist, time <admin_cens)
prob_cens <- sum(early$prob)
censor_dist <- rbind( early,
       c(admin_cens, 1-prob_cens)
)

add_censor_sample <- function(data, censor_dist){
  n <- nrow(data)
  data$censor_time <- sample(censor_dist$time, n, replace=TRUE, prob=censor_dist$prob)
  data <- mutate(data, 
                 time=pmin(time, censor_time),
                 status=ifelse(censor_time<=time, 0, status)
                 )
  data
}
#test
if(FALSE){
  survfit(Surv(time, status)~rx, data=ans$short) %>% plot(mark.time=TRUE)
  ans2 <- add_censor_sample(ans$short, censor_dist)
  survfit(Surv(time, status)~rx, data=ans2) %>% plot(mark.time=TRUE)
}




