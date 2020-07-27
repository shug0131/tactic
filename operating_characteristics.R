command_line_args = commandArgs(trailingOnly=TRUE)
#### command line arguments
if (length(command_line_args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

scenario_arg <- ifelse(is.na(command_line_args[1]),"null",command_line_args[1])


# need a symbolic link in the tactic directory
system("ln -s ~/rds/hpc-work/ ~/tactic/results")


output <- list.files(path=paste0("./results/", scenario_arg), pattern="inference*")

df <- readRDS(paste0("results/",scenario_arg,"/inference_results_100.rds"))

results_df <- cbind(df[0,], data_set=integer(0))
for( file in output){
  index <- gsub(pattern = "inference_results_(\\d+).rds","\\1", file)
  df <- readRDS(paste0("results/",scenario_arg,"/",file))
  df$data_set <- as.integer(index)
  results_df <- rbind(results_df, df)
}

#pick out first row that has either reject or accept within each arm and dataset

library(dplyr)
library(tidyr)
library(magrittr)


oc_group_seq <- function(family){
family_var <- enquo(family)
df <- results_df %>% group_by(data_set, arm) %>% 
  filter( !!family_var != "continue") %>% 
  filter(row_number()==1)

df$interim <- factor(df$interim, levels=1:3, labels=paste0("stage_",1:3))

df2 <- df %>% group_by(arm) %>% 
  summarise(power=mean( !!family_var=="reject"),
            bias=mean(freq_bias),
            mse=sqrt(mean(freq_bias^2)),
            coverage=mean(freq_coverage),
            bayes_mse=sqrt(mean(bayes_bias^2)),
            # Note the order is important - bad code
            bayes_bias=mean(bayes_bias),
            bayes_coverage=mean(bayes_coverage),
            rule=as_label(family_var)
            ) #%>%  print

n <- xtabs(~arm, data=df)
stage <- sweep(xtabs(~arm+interim, data=df),1,n,"/") %>% as.array()
#print(stage)
#cbind(df2, stage)
stage %<>% as.data.frame %>% 
  spread(key="interim", value="Freq")
left_join(df2, stage, by="arm")
}



triangular <- oc_group_seq(triangular)
pocock <- oc_group_seq(pocock)
obf <- oc_group_seq(obf)



oc_bayes <- function(prior, threshold=0.95){
  prior_var <- enquo(prior)
  df <- results_df %>% group_by(data_set, arm) %>% 
    filter( !!prior_var > threshold | interim==3) %>% 
    filter(row_number()==1)
  df$interim <- factor(df$interim, levels=1:3, labels=paste0("stage_",1:3))
  
  
  df2 <- df %>% group_by(arm) %>% 
    summarise(power=mean( !!prior_var> threshold),
              bias=mean(freq_bias),
              mse=sqrt(mean(freq_bias^2)),
              coverage=mean(freq_coverage),
              bayes_mse=sqrt(mean(bayes_bias^2)),
              # Note the order is important - bad code
              bayes_bias=mean(bayes_bias),
              bayes_coverage=mean(bayes_coverage),
              rule=as_label(prior_var)
    ) #%>%  print
  
  n <- xtabs(~arm, data=df)
  stage <- sweep(xtabs(~arm+interim, data=df),1,n,"/")
  stage %<>% as.data.frame %>% 
    spread(key="interim", value="Freq")
  left_join(df2, stage, by="arm")
}

vague <- oc_bayes(vague)
pessimistic <- oc_bayes(pessimistic)
optimistic <- oc_bayes(optimistic)

oc_df <- Reduce(rbind, list(triangular, pocock, obf, vague, pessimistic, optimistic))

oc_df$scenario=scenario_arg

print(oc_df)
saveRDS(oc_df, file=paste0("./results/",scenario_arg,"/oc_results.rds"))




            
