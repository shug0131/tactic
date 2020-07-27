library(magrittr)
library(tidyr)
library(ggplot2)
library(dplyr)

file_list <- list.files("results/null")
index <- gsub("coxfit_(\\d+)\\.rds","\\1", file_list) %>% as.numeric 
index <- index[complete.cases(index)]

pval <- array(0, dim=c(length(index),2), dimnames = list(NULL, c("E1","E2")))
for( i in 1:length(index)){
  mod <- readRDS(file=paste0("results/null/coxfit_",index[i],".rds"))
  pval[i,] <- mod$coefficients[1:2, "p"]
}

pval %>% data.frame %>% gather(key="arm",value="pvalue") %>% 
  ggplot(aes(x=arm,y=pvalue))+geom_boxplot()

plot(pval)

boxplot(pval[,"E2"])

## Posterior Probabilities



file_list <- list.files("results/null")
index <- gsub("posterior_probs_(\\d+)\\.rds","\\1", file_list) %>% as.numeric 
index <- index[complete.cases(index)]

prob_df <- function(probs, index){
cbind(Reduce(rbind,probs),
prior=rep(names(probs),rep(nrow(probs[[1]]),length(probs))),
range=rep(row.names(probs[[1]]), length(probs)),
simulation=index
)
}





probs <- readRDS(file="results/null/posterior_probs_NA.rds")
all_prob <- prob_df(probs,NA)
all_prob <- all_prob[0,]
for(i in 1:length(index)){
  probs <- readRDS(file=paste0("results/null/posterior_probs_",index[i],".rds"))
  all_prob <- rbind(all_prob,  prob_df(probs, index[i]))
}
summary(all_prob)


#index <- with(all_prob, prior=="vague" & range=="similarity")
all_prob %>% gather(key="arm", value="prob",rxE1,rxE2) %>% 
  ggplot(aes(x=range, y=prob))+geom_boxplot()+facet_grid(prior~arm)


## Power

#???
