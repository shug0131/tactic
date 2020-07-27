setwd("~/Documents/Work/tactic")
.libPaths("~/Documents/Work/tactic/library")
library(mvtnorm)
library(magrittr)

# Design a triangular test - overall optimal
#  Fix the interim analysis points
# The geometric shape on the Z-scale is fixed - a triangle
# Find the one parmeter ( two y-intercepts +/- the same) to give significnace overall


n <- c(125,229,469)
cor <- outer(n,n, function(x,y){sqrt(pmin(x,y)/pmax(x,y))})

bounds_triangle <- function(a,n){
  lower=-a*(1 - 3*(n/max(n)))
  upper=a*(1 + (n/max(n)))
  cbind(lower, upper)
}

bounds_pocock <- function(a,n){
  cbind(rep(-a,length(n)), rep(a, length(n)))
}

bounds_obf <- function(a,n){
  upper <- a*sqrt(max(n)/n)
  cbind(-upper, upper)
}


bounds <- function(a,n, family=c("triangular","pocock","obf")){
  family=match.arg(family)
  switch(family,
         triangular=bounds_triangle(a,n),
         pocock=bounds_pocock(a,n),
         obf=bounds_obf(a,n)
         )
}

bounds(1,n, "obf")


power <- function(a,n, theta=0, target=0, family="triangular"){
  require(mvtnorm)
  cor <- outer(n,n, function(x,y){sqrt(pmin(x,y)/pmax(x,y))})
  bound <- bounds(a,n, family=family)
  lower <- bound[,1]
  upper <- bound[,2]
  reject <- rep(0,length(n))
  accept <- rep(0,length(n))
  reject[1] <- 1-pnorm(upper[1], mean=theta*sqrt(n[1]))
  accept[1] <- pnorm(lower[1], mean=theta*sqrt(n[1]))
  for(analysis in 2:length(n)){
    reject[analysis] <- pmvnorm(
      lower=c(lower[1:(analysis-1)], upper[analysis]),
      upper=c(upper[1:(analysis-1)], Inf),
      mean=theta*sqrt(n[1:analysis]),
      corr= cor[1:analysis, 1:analysis]
    )
    accept[analysis] <- pmvnorm(
      lower=c(lower[1:(analysis-1)], -Inf),
      upper=c(upper[1:(analysis-1)], lower[analysis]),
      mean=theta*sqrt(n[1:analysis]),
      corr= cor[1:analysis, 1:analysis]
    )
  }
  output <- sum(reject) - target
  attr(output,"individual_reject") <- reject
  attr(output,"individual_accept") <- accept
  output
}

boundary_list <- list()

a <- uniroot(power, c(1,4),n=n, target=0.025, family="pocock")$root
#e_tri <- power(a, n, family="poc")
boundary_list <- c(boundary_list, list("pocock"=bounds(a,n, family="p")))
a <- uniroot(power, c(1,4),n=n, target=0.025, family="triangular")$root
boundary_list <- c(boundary_list, list("triangular"=bounds(a,n, family="triangular")))
a <- uniroot(power, c(1,4),n=n, target=0.025, family="obf")$root
boundary_list <- c(boundary_list, list("obf"=bounds(a,n, family="obf")))

saveRDS(boundary_list, file="boundary_list.rds")

e_tri
f_tri <- power(a, n, theta=abs(log(0.6))/sqrt(1/0.32+1/0.21))
f_tri

#type 1 and Type 2 errors
cbind(attr(e_tri,"individual_reject") %>% cumsum, attr(f_tri,"individual_accept") %>% cumsum)

plot(n, bounds(a,n)[,1])
points(n, bounds(a,n)[,2])
# reverse engineer teh pwoer calculation
#HR =0.6, with event rateof 0.2 amd 0.12

#events <- (1-c(0.68,0.79))*229 
#1-pnorm( qnorm(0.975) - abs(log(0.6))/sqrt(sum(1/events)), sd = 1)
#points(log(n), log(cumsum(attr(e,"individual_reject"))), col="red")
#plot(log(n), log(cumsum(attr(f,"individual_accept"))))

# Error spendign approach with alpha*(V)^0.5, and beta*(V)^0.5 for upper an lower error.

# So Trianlgle - fix bounds as 1 parameter, solve to give final error
# Error-Spend- fix error, and calculate bounds based on increment in error & type 2

cor
mu <- abs(log(0.6))/sqrt(1/0.32+1/0.21)*sqrt(n)
e <- 0.025*(n/max(n))^0.5
f <- 0.2*(n/max(n))^0.5
lower <- upper <- rep(0,3)
upper[1] <- qnorm(e[1], lower.tail = FALSE)
lower[1] <- qnorm(f[1], mean=mu[1], lower.tail = TRUE)

find_boundary_ee <- function(x, lower_bound, upper_bound, mu, cor, target=0, lower.tail=TRUE){
  require(mvtnorm)
  lower <- c(lower_bound, ifelse(lower.tail,-Inf, x))
  upper <- c(upper_bound, ifelse(lower.tail,x, Inf))
  pmvnorm( lower, upper, mean = mu, corr=cor)-target
}

for( analysis in 2:length(n)){
  upper[analysis] <- uniroot(find_boundary_ee, interval=c(upper[analysis-1],100),
                             lower_bound=lower[1:(analysis-1)],
                             upper_bound=upper[1:(analysis-1)],
                             mu=0,
                             cor=cor[1:analysis,1:analysis],
                             target=e[analysis]-e[analysis-1],
                             lower.tail=FALSE, extendInt = "yes"
  )$root
  lower[analysis] <- uniroot(find_boundary_ee, interval=c(lower[analysis-1], upper[analysis-1]),
                             lower_bound=lower[1:(analysis-1)],
                             upper_bound=upper[1:(analysis-1)],
                             mu=mu[1:analysis],
                             cor=cor[1:analysis,1:analysis],
                             target=f[analysis]-f[analysis-1],
                             lower.tail=TRUE, extendInt = "yes"
  )$root
                             
}



cbind(lower, upper)
# ignor the final lower bound - must set the same, ie. the upper for significnace.
cbind(e,f)
points( n, lower, col="red")
points( n, upper, col="red")
