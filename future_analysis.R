
## @knitr future_analysis


# For each row of pred_DataA EudraCT results user with this username and email address cannot be found
future_analysis <- function(prediction, newdata, olddata, 
                            arms=NULL, 
                            n_total=nrow(newdata)+nrow(olddata)
                            ){
  
  
  df_new <- add_sim_status(newdata, prediction,site)
  # add_censor_sample()
  df_new <- add_censor_sample(df_new$short, censor_dist)
  # combine with the original data
  common_names <- intersect(names(df_new), names(olddata))
  df_all <- rbind(olddata %>% as.data.frame %>%  select(all_of(common_names)),
                  df_new %>% as.data.frame %>%  select(all_of(common_names))
  )
  n_arms <- length(arms)
  if(!is.null(arms)){
    df_all <- subset(df_all, rx %in% arms)
    #n_arms <- nlevels(olddata$rx)
  }
  if( min(n_total)<nrow(olddata)*n_arms/nlevels(olddata$rx)){
    stop("you are throwing out the existing data!")
  }
  
  
  ans <- array(NA, c(length(n_total), 4))# auto-fill 2*(n_arms-1))) 
  for( row in 1:length(n_total)){
    df_subset <- df_all[1:n_total[row],]
    df_subset$rx <- droplevels(df_subset$rx)
    # consider which arms to use, how many new patients.
    # calculate coefs and SE using coxph - but leave open to stan_glm()...
    analysis <- coxph(Surv(time,status)~rx+frailty.gaussian(site), data=df_subset)
    # loop across rows and then calculate OCs: power, bias, SE, coverage...
    se <- analysis$var %>% diag %>% sqrt
    coef <- analysis$coef
    #coef1, se1, coef2, se2,   ## even if we repeat and 1=2 
    if(length(coef)==2){
    ans[row,] <- c(coef[1],se[1],coef[2],se[2])
    } else{
      ans[row,] <- rep(c(coef,se), 2)
    }
  }
 
  ans
}



combinations <- function(x){
  m <- length(x)
  n <- 2^m
  ans <- vector("list",n-1)
  for(i in 1:(n-1)){
    index <- binary(i)
    pad <- m-length(index)
    index <- c(as.logical(index),rep(FALSE,pad))
    ans[[i]] <- x[index]
  }
  ans
}

#combinations(c("A","B","C"))

binary <- function(x){
  if(x<2){
    return(x)
  }else{
    alpha <- x %% 2
    x <- x %/% 2
    return( c(alpha, Recall(x)))
  }
}

#as.logical(binary(4))


power <- function(X){
    z1 <- X[,1,]/X[,2,]
    z2 <-  X[,3,]/X[,4,]
    #1-side alpah/2 test
  result1 <- pnorm(z1)<0.025
  result2 <- pnorm(z2)<0.025
  result12 <- result1 | result2
  power1 <- apply(result1,1, mean)
  power2 <- apply(result2,1, mean)
  power12 <- apply(result12,1, mean)
  cbind(power1,power2, power12)
}
