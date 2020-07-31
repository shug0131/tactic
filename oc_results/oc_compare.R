scenario_df <- read.csv(file="../scenario_df.csv", stringsAsFactors = FALSE)
for(scenario in scenario_df$scenario){
  df <- readRDS(paste0(scenario,".rds"))
  assign(scenario, df)
}

obj_list <- lapply(scenario_df$scenario, as.name)
all <- Reduce(rbind, lapply(obj_list,eval))

rm(list=scenario_df$scenario)
rm(obj_list, scenario)

.libPaths("../library")
library(tidyr)
library(dplyr)
all <- left_join(all, scenario_df, by="scenario")

library(ggplot2)
library(magrittr)
all %<>% mutate(hr=ifelse( arm=="rxE1",hr1,hr2))

all_long <- gather(all,key="metric",value="value", power, bias,mse, coverage, bayes_mse, bayes_bias, bayes_coverage)


all_long %>% ggplot(aes(y=value,x=hr, group=rule, colour=rule))+geom_line()+
  facet_grid(arm~metric)
