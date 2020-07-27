library(survival)
library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)

# Simulate data, or function to do so
# 14 day semi-discrete failure time
# Rates based on power calculations
# small degree of random censoring before day 14
source("sim_failure.R")

df <- sim_failure(n_per_arm = 125,control_surv_rate = 0.7, hr=c(0.7,1.3)) %>% 
  add_censor(admin_cens_time = 14, drop_out=0.1)

#add in site , Any other covars?
#df$site <- rep(1:25,nrow(df)/25) %>% factor


# surv_split the data

fail_times <- df %>% filter(status==1) %>% select(time) %>% unique() %>% unlist
df_long <- survSplit(Surv(time, status)~rx+site, df, cut=fail_times)  


# Frequentist equivalent inference check.
glmer(status~-1+strata(time)+rx +(1|site), family=poisson, data=df_long) %>% summary
fit_cox <- coxph(Surv(time,status)~rx+frailty.gaussian(site), data=df, ties = "breslow") 
fit_cox %>% summary
survfit(Surv(time,status)~rx, data=df) %>% plot(lty=1:3)
legend(0,0.6,lty=1:3, legend=levels(df$rx))


#fit poisson glm bayesian 

library(rstanarm)
#library(brms)
options(mc.cores = parallel::detectCores())
library(bayesplot)
theme_set(bayesplot::theme_default())
#library(brms)

# work on HR priors. vague, strong postivie & negative

source("priors.R")

prior_list <- list("vague"=vague, "optimistic"=optimistic, "pessimistic"=pessimistic)

effects <- c("rxE1","rxE2")

if(FALSE){
  fits <- lapply(prior_list,
                 function(my_prior){
                   stan_glm(status~-1+strata(time)+rx, data=df_long, family=poisson,
                            chains=2, iter=1000,QR=FALSE,
                            prior =my_prior)
                 }
  )




#fit %>% summary

# extract inference, 
## HR, probabilities of ranges,
#prior_summary(fit)

posterior_vs_prior(fits[["optimistic"]],pars=effects)
posterior_vs_prior(fits$pessimistic,pars=effects)
posterior_vs_prior(fits$vague,pars=effects)

lapply(fits, posterior_interval, pars=effects)
sapply(fits, coef)[effects,]
}

# focus on the vage prior model now
fit <- stan_glmer(status~-1+strata(time)+rx+(1|site), data=df_long, family=poisson,
                chains=2, iter=1000,QR=FALSE,
                prior =vague)
fit %>% summary
#transform back to HR scale?
plot(fit, pars=effects, plotfun = "areas", prob = 0.5)#, prob_outer = 0.9)
#transform back to HR scale?
plot(fit, pars=effects, plotfun = "areas", transformations="exp")+
  labs(x="Hazard Ratio")+scale_y_discrete(labels=c("E1 vs Control","E2 vs Control"))
plot(fit, pars=effects, "hist")
plot(fit, pars=effects, "violin")
plot(fit, pars=effects, "dens_overlay")
plot(fit, pars=effects, "trace")
plot(fit, pars=effects, "acf_bar")

posterior <- as.array(fit)
#dimnames(posterior)

prob_interval <- function(X, interval){
  Y <- interval[1]< X & X <interval[2]
  mean(Y)
}

apply(posterior, 3, FUN=prob_interval, interval=c(-Inf, log(0.8)))[effects]
apply(posterior, 3, FUN=prob_interval, interval= log(c(0.8,1.1)))[effects]
apply(posterior, 3, FUN=prob_interval, interval= c(log(1.1),Inf))[effects]





# posterior_predictive simulation up to large N, all arms
# bootstrap covariates
# random treatment
# full follow-up
n_new <- 1500
source("new_data.R")
# this gets us output : df_new_long, censor_dist,  several functions add_*()
pred_data <- posterior_predict(fit, newdata = df_new_long)
#convert back to times, rather than long counting process format, takign the first non-zero status
source("future_analysis.R")


# Now loop across the 3 combinations (Check there are 3 arms...)
# and a few choices of SS
#### START
# now loop


# all possible combination of arms, with Cx allways taken forwards
arms <- combinations(levels(df$rx)[-1])
arms <- lapply(arms, function(x){c("Cx",x)})


ref_ss <- c(458,687,938,1407)




power_list <- lapply(arms, power, n_total=ref_ss)
#power_list
#ref_ss
arms_char <- sapply(arms, paste, collapse=", ")
power_df_adj <- data.frame(power=unlist(power_list),
                       ss=rep(ref_ss, length(arms)),
                       arms=rep(arms_char, rep(length(ref_ss),length(arms)))
                       )

ggplot(power_df_adj, aes(x=ss,y=power, group=arms, colour=arms))+
  geom_line()+geom_point()

##END

## Next step woudl be to consider decision rules based on these predictive powers
# then simulate actual data following the rules
# loop around teh strategy

## NEED to start timing and optimising code at these steps.
