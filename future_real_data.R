
##
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




# need a symbolic link in the tactic directory


system("ln -s ~/rds/hpc-work/ ~/tactic/results")


###  Pull out old data into standard survival data set
## This might be saved as a data from in future versions of code
## But I didn't do this D'Oh, so grab it from the Rstanarm object

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

model_index <- task_id


filename <- paste0("results/",scenario_arg,"/fits_image_", model_index,".rds")
fit <- readRDS(file=filename)
library(tidyr)
library(dplyr)
df <- fit$optimistic$data
df$id <- 0
id <- 0
for( row in 1:nrow(df)){
  if(df[row,"time"]==1){id <- id+1}
  df[row,"id"] <- id
}
#pull out the last row for each id
df <- df %>% group_by(id) %>% slice(n()) %>% ungroup %>% as.data.frame()


## Get the parameters needed to generate more data
## Again coudl be saved, but easier to just repeat the same code from the interim data/analysis

## BAD CODE
#scenario_arg <- "null"

scenario_df <- read.csv(file="scenario_df.csv", stringsAsFactors = FALSE)
inputs <- subset(scenario_df,scenario==scenario_arg) %>% as.list()
inputs$hr <- with(inputs, c(hr1,hr2))

# HARD CODING
inputs$n_per_arm <- 469-125  # The max sample size minus the data already generated
 
## Generate future data for 2/3 arms at each of the 2 future analysis points.
source("sim_failure.R")
index <- which(names(inputs) %in% names(formals(sim_failure)))
df_new <- do.call(sim_failure, inputs[index]) #ignore the warning
index <- which(names(inputs) %in% names(formals(add_censor)))
df_new <- do.call(add_censor, c(list(df_new), inputs[index]))
rm(index)

## Produce coxph analysis at each 4 combinations

# blocks of treatments, size 344, then site is 1:10 replicated lots of times.
# work out which rows are in teh second interim, but not the final


#HARD COding
n2 <- 229-125
# c(  1:n2, 344+1:n2, 2*344+1:n2)
index2 <-rep(1:n2, 3)+rep(0:2*344, rep(n2,3))

df_interim2 <- rbind(
  df[,c("rx","site","time","status")],
  df_new[index2, c("rx","site","time","status")]
)
##  cox fit,  and triangular test version. 
library(survival)
cox_interim2 <- coxph(Surv(time, status)~rx+ frailty.gaussian(site), data=df_interim2) %>% summary

df_final <- rbind(
  df[,c("rx","site","time","status")],
  df_new[, c("rx","site","time","status")]
)
cox_final <- coxph(Surv(time, status)~rx+ frailty.gaussian(site), data=df_final) %>% summary

filename <- paste0("results/null/coxfit_", model_index,".rds")
cox_interim1 <- readRDS(file=filename)

## Run Bayes model at each 4 combinations

# INterim 2 bayes analysis and set up for baye/parallel processing.

fail_times <- df_interim2 %>% filter(status==1) %>% select(time) %>% unique() %>% unlist

df_long_interim2 <- survSplit(Surv(time, status)~rx+site, df_interim2, cut=fail_times)  
library(rstanarm)
#library(brms)
number_chains=2
### CHECK FOR SLURM ARRAY ETC
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


#library(bayesplot)
source("priors.R")
prior_list <- list("vague"=vague, "optimistic"=optimistic, "pessimistic"=pessimistic)
effects <- c("rxE1","rxE2")

system.time(
  fits_interim2 <- foreach(i = 1:length(prior_list),
                  .inorder = TRUE,
                  .packages=c("rstanarm","survival"),
                  .export=c("df_long_interim2")
  ) %my_do% {
    stan_glmer(status~-1+strata(time)+rx+(1|site), data=df_long_interim2, family=poisson,
               chains=number_chains, iter=3000/number_chains,QR=FALSE, thin=2,
               prior =prior_list[[i]])
  }
)


names(fits_interim2) <- names(prior_list)

# Final bayes analaysis

fail_times <- df_final %>% filter(status==1) %>% select(time) %>% unique() %>% unlist
df_long_final<- survSplit(Surv(time, status)~rx+site, df_final, cut=fail_times)  
system.time(
  fits_final <- foreach(i = 1:length(prior_list),
                           .inorder = TRUE,
                           .packages=c("rstanarm","survival"),
                           .export=c("df_long_final")
  ) %my_do% {
    stan_glmer(status~-1+strata(time)+rx+(1|site), data=df_long_final, family=poisson,
               chains=number_chains, iter=3000/number_chains,QR=FALSE, thin=2,
               prior =prior_list[[i]])
  }
)
names(fits_final) <- names(prior_list)



## Put the estimates, SE, and hypothesis test into a data.frame.
# frequentest unadjusted  and triangular test version, and 

#source("group_sequential_design.R")
# but don't need to do this repeatedly. The set of boundaries are saved

print("Line 140")

boundary_list <- readRDS("boundary_list.rds")

#1-sided p-values, negative log hr is good. 
inference <- function(fit, interim,arm, boundary_list){
  # only look at one arm at a time
  estimates <- coef(fit)[arm,]
  coef <- estimates["coef"]
  se <- estimates["se(coef)"]
  z <- -coef/se #negative log hr is good. 
  pvalue <- 1-pnorm(z)
  decisions <- sapply(boundary_list, function(x,z,interim){
    if(z< x[interim,1]){output <- "accept"}
    if(x[interim,1]<= z & z <= x[interim,2]){output <- "continue"}
    if( x[interim,2]< z){ output <- "reject"}
    if( interim==3 & output=="continue"){output <- "accept"}
    output
  },z=z, interim=interim)
  data.frame(arm=arm,interim=interim, "coef"=coef, "se"=se, z=z, pvalue=pvalue,
        pocock=decisions["pocock"],
        triangular=decisions["triangular"],
        obf=decisions["obf"], stringsAsFactors = FALSE
        )
}

inference_df <- expand.grid(arm=c("rxE1","rxE2"), interim=1:3, stringsAsFactors = FALSE)
inference_df %<>% mutate( fit=c("cox_interim1","cox_interim2","cox_final")[interim])

print("Line 169")

nrow <- nrow(inference_df)
results_df <- data.frame(arm=character(nrow),
                         interim=integer(nrow),
                         coef=numeric(nrow),
                         se=numeric(nrow),
                         z=numeric(nrow),
                         pvalue=numeric(nrow),
                         pocock=character(nrow),
                         triangular=character(nrow),
                         obf=character(nrow), stringsAsFactors = FALSE
                         )
for( row in 1:nrow(inference_df)){
  x <- as.list(inference_df[row,])
  results_df[row,] <- do.call("inference", list(fit=as.name(x$fit),interim=x$interim, arm=as.character(x$arm),
        boundary_list=boundary_list) )
}

print("Line 188")
# bayesian p-values.
bayes_p <- function(fit){
  posterior <- as.array(fit)
  effects=c("rxE1","rxE2")
  apply(posterior, 3, FUN=function(x){mean(x<0)})[effects]
}

bayes_int1 <- lapply(fit, bayes_p)
bayes_int2 <- lapply(fits_interim2, bayes_p)
bayes_final <- lapply(fits_final, bayes_p)

bayes_results <- cbind( bind_rows(c(bayes_int1,bayes_int2, bayes_final)),
       prior=rep(names(bayes_final),3),
       interim=rep(1:3,rep(3,3))
)

print("Line 205")

bayes_results %<>% gather(key = "arm", value = "prob", -prior, -interim ) %>% 
  spread(key=prior, value=prob)
results_df %>% merge(bayes_results, by=c("interim","arm"))

#get bayes estimates and se 

bayes_est <- function(fit){
  data.frame(arm=effects,
  bayes_coef=fixef(fit$vague)[effects],
  bayes_se=se(fit$vague)[effects],
  stringsAsFactors = FALSE
  )
}

print("Line 221")

bayes_est_df <- lapply(list(interim1=fit, interim2= fits_interim2, final=fits_final), bayes_est) %>% 
  Reduce(rbind,.)
bayes_est_df$interim <- rep(1:3,rep(2,3))

output <- results_df %>% merge(bayes_results, by=c("interim","arm")) %>% 
  merge(bayes_est_df, by=c("interim","arm"))

output$model_index <- model_index
output$scenario <- scenario_arg 

# coverage 
log_hr <- filter(scenario_df, scenario==scenario_arg) %>% 
  select(hr1,hr2) %>% unlist %>% log

print("Line 237")

output %<>% mutate( 
  freq_coverage=coef - 1.96 * se < log_hr & log_hr < coef + 1.96 * se,
  bayes_coverage=bayes_coef - 1.96 * bayes_se < log_hr & log_hr < bayes_coef + 1.96 * bayes_se,
  freq_bias = coef-log_hr,
  bayes_bias= bayes_coef - log_hr
)

filename=paste0("results/",scenario_arg,"/inference_results_",model_index,".rds")
saveRDS(output, file=filename)
print(scenario_arg)
stopCluster(cl)
