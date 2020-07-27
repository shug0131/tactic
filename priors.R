# define a selection of priors for the log-hazard coefficients
# use normal distribution
# find hyper-parameters that match quantiles that are specified

## @knitr priors

inverse_normal <- function(quantiles, probs=c(0.05,0.95)){
  A <- cbind(1,qnorm(probs))
  ans <- solve(A, quantiles)
  names(ans) <- c("location","scale")
  ans
}

#test
# inverse_normal(c(-1,2))
# pnorm(c(-1,2),mean=0.5, sd=0.9119)

#control hazard function

# set of cumhaz increments per day for days 1- 14
# if exponential with 14-day survival S0, then log-haz= log(-ln(1-S0)/14)

(quantile_log_haz <- log(-log(c(0.95,0.05))/14))
(base_log_haz_pars <- inverse_normal(quantile_log_haz))

(log_hr <- inverse_normal(log(c(0.5,1/0.5))))

#number of events in a 2 arm trial
4/log_hr["scale"]^2 # 23


vague <- normal(location=rep( c(base_log_haz_pars["location"], log_hr["location"]), c(14,2)),
                scale=rep( c(base_log_haz_pars["scale"], log_hr["scale"]), c(14,2)),
                autoscale=FALSE
                )

#optimistic
(log_hr <- inverse_normal(log(c(0.7,1)), probs=c(0.5,0.95)))

#number of events in a 2 arm trial
4/log_hr["scale"]^2 # 85

optimistic <- normal(location=rep( c(base_log_haz_pars["location"], log_hr["location"]), c(14,2)),
                scale=rep( c(base_log_haz_pars["scale"], log_hr["scale"]), c(14,2)),
                autoscale=FALSE
)


#pessimistic
(log_hr <- inverse_normal(log(c(0.7,1)), probs=c(0.05,0.5)))

#number of events in a 2 arm trial
4/log_hr["scale"]^2 # 85

pessimistic <- normal(location=rep( c(base_log_haz_pars["location"], log_hr["location"]), c(14,2)),
                    scale=rep( c(base_log_haz_pars["scale"], log_hr["scale"]), c(14,2)),
                    autoscale=FALSE
)

