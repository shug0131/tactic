command_line_args = commandArgs(trailingOnly=TRUE)
#### command line arguments
if (length(command_line_args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

scenario_df <- read.csv(file="scenario_df.csv", stringsAsFactors = FALSE)
#print(scenario_df)
scenario_arg <- ifelse(is.na(command_line_args[1]),"one_winner",command_line_args[1])
# number of cores to use internally for parallelism
parallel_arg <- ifelse(is.na(command_line_args[2]),0,command_line_args[2])

#########



.libPaths("library")
library(survival)
library(magrittr)
library(dplyr)

library(doParallel)
library(foreach)

if( 1 < parallel_arg){ 
  ncores <- min(parallel::detectCores(), parallel_arg)
  `%my_do%` <- `%dopar%`
} else{
  ncores <- 1
  `%my_do%` <- `%do%`
}

cl<-makeCluster(ncores)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("library"))
getDoParWorkers()
source("sim_failure.R")
inputs <- subset(scenario_df,scenario==scenario_arg) %>% as.list()
inputs$hr <- with(inputs, c(hr1,hr2))

index <- which(names(inputs) %in% names(formals(sim_failure)))
df <- do.call(sim_failure, inputs[index]) #ignore the warning
index <- which(names(inputs) %in% names(formals(add_censor)))
df <- do.call(add_censor, c(list(df), inputs[index]))
rm(index)
#survfit(Surv(time, status)~rx, data=df) %>% plot(lty=1:3)
#coxph(Surv(time, status)~rx, data=df) %>% summary()
#df <- sim_failure(n_per_arm = 125,control_surv_rate = 0.7, hr=c(0.7,1.3)) %>% 
#  add_censor(admin_cens_time = 14, drop_out=0.1)

fail_times <- df %>% filter(status==1) %>% select(time) %>% unique() %>% unlist
df_long <- survSplit(Surv(time, status)~rx+site, df, cut=fail_times)  
library(rstanarm)
#library(brms)
number_chains=2
options(mc.cores=ncores)
#library(bayesplot)
source("priors.R")
prior_list <- list("vague"=vague, "optimistic"=optimistic, "pessimistic"=pessimistic)
effects <- c("rxE1","rxE2")

system.time(
fits <- foreach(i = 1:length(prior_list),
                .inorder = TRUE,
                .packages=c("rstanarm","survival"),
                .export=c("df_long")
) %my_do% {
  stan_glmer(status~-1+strata(time)+rx+(1|site), data=df_long, family=poisson,
             chains=number_chains, iter=3000/number_chains,QR=FALSE, thin=2,
             prior =prior_list[[i]])
}
)


names(fits) <- names(prior_list)

#library("bayesplot")
#plot(fits[["vague"]],pars=effects, plotfun="areas")

#save it 
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

file_name <- paste0("/home/sjb277/rds/hpc-work/",scenario_arg,"/fits_image_",task_id,".rds")
saveRDS(fits, file=file_name)

# Calculate the the posterior probabilities of various regions


posterior_probabilities <- function(fit, effects=c("rxE1","rxE2")){
 
  posterior <- as.array(fit)
 
  prob_interval <- function(X, interval){
    Y <- interval[1]< X & X <interval[2]
    mean(Y)
  }
  efficacy <- apply(posterior, 3, FUN=prob_interval, interval=c(-Inf, 0))[effects]
  moderate_efficacy <- apply(posterior, 3, FUN=prob_interval, interval=c(-Inf, log(0.8)))[effects]
  similarity <- apply(posterior, 3, FUN=prob_interval, interval= log(c(0.8,1.25)))[effects]
  harm <- apply(posterior, 3, FUN=prob_interval, interval= c(0,Inf))[effects]
  X <- rbind(efficacy, moderate_efficacy, similarity, harm) %>% as.data.frame
  X <- cbind(HR=c("<1","<0.8","0.8-1.25",">1"),X)
  return(X)
}

prob_list <- lapply(fits, posterior_probabilities)

file_name <- paste0("/home/sjb277/rds/hpc-work/",scenario_arg,"/posterior_probs_",task_id,".rds")
saveRDS(prob_list, file=file_name)

#basic coxph model
cox_fit <- coxph(Surv(time,status)~rx+frailty.gaussian(site), data=df) %>% summary
file_name <- paste0("/home/sjb277/rds/hpc-work/",scenario_arg,"/coxfit_",task_id,".rds")
saveRDS(cox_fit, file=file_name)



# Now do the power calculations 


n_new <- 1500
source("new_data.R")
# this gets us output : df_new_long, censor_dist,  several functions add_*()
source("future_analysis.R")

system.time(
  pred_data_list <- lapply(fits, posterior_predict, newdata = df_new_long)
)

arms <- combinations(levels(df$rx)[-1])
arms <- lapply(arms, function(x){c("Cx",x)})
ref_ss <- c(458,687,938,1407)

library(abind)

abind3 <- function(...){abind(..., along=3)}

print("starting power loops")

system.time(
  power_list <- foreach(pred_index = 1:length(pred_data_list)) %:% 
    foreach( arm_index = 1:length(arms)) %:%
    foreach(row =1:nrow(pred_data_list[[1]]),
            .inorder = FALSE, .combine=abind3, .final=power,
            .packages = c("survival","magrittr"),
            .export=c("df","df_new_long","ref_ss")
    ) %my_do% {
      future_analysis(prediction=pred_data_list[[pred_index]][row,],
                      newdata=df_new_long, olddata=df,
                      arms=arms[[arm_index]], n_total=ref_ss
      )
    }
)




power_df <- cbind(
  Reduce(rbind, unlist(power_list, recursive = FALSE)),
  expand.grid(n_total=ref_ss,arms=sapply(arms, paste, collapse=", "),prior=names(fits))
)

#library(ggplot2)
#ggplot(power_list_df, aes(x=n_total, y=power12, group=prior, colour=prior))+
#  geom_point()+geom_line()+facet_grid(arms~.)

file_name <- paste0("/home/sjb277/rds/hpc-work/",scenario_arg,"/power_",task_id,".rds")
saveRDS(power_df, file=file_name)


stopCluster(cl)





