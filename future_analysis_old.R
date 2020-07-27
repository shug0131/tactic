
## @knitr future_analysis


# For each row of pred_Data
future_analysis <- function(prediction, newdata, olddata, 
                            arms=NULL, 
                            n_total=nrow(newdata)+nrow(olddata),
                            parallel=TRUE){
  # add_sim_stats()
  require(doParallel)
  require(foreach)
  
  
  #prediction <- pred_data[1,]# FUN=future_analysis, 
  #newdata=df_new_long
  #o#lddata=df
  
  df_new <- add_sim_status(newdata, prediction,rx,site)
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
  
  if(!parallel){
  ans <- array(NA, c(length(n_total), 2*(n_arms-1))) 
  for( row in 1:length(n_total)){
    df_subset <- df_all[1:n_total[row],]
    df_subset$rx <- droplevels(df_subset$rx)
    # consider which arms to use, how many new patients.
    # calculate coefs and SE using coxph - but leave open to stan_glm()...
    analysis <- coxph(Surv(time,status)~rx+frailty.gaussian(site), data=df_subset)
    # loop across rows and then calculate OCs: power, bias, SE, coverage...
    se <- analysis$var %>% diag %>% sqrt
    ans[row,] <- c(analysis$coefficients, se)
  }
  }
  
  if(parallel){
  #cl<-makeCluster(2)
  #registerDoSNOW(cl)
  ans <- foreach( row = 1:length(n_total), .packages = c("survival","magrittr"),
                   .combine=rbind, .export=c("df_all","n_total"), .inorder=TRUE
                   )  %dopar% {
    df_subset <- df_all[1:n_total[row],]
    df_subset$rx <- droplevels(df_subset$rx)
    # consider which arms to use, how many new patients. 
    # calculate coefs and SE using coxph - but leave open to stan_glm()...
    analysis <- coxph(Surv(time,status)~rx+frailty.gaussian(site), data=df_subset)
    # loop across rows and then calculate OCs: power, bias, SE, coverage...
    se <- analysis$var %>% diag %>% sqrt
    c(analysis$coefficients, se)
  }
  #stopCluster(cl)
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


power <- function(arms,n_total){
  require(doParallel)
  require(foreach)
  # apply doesn't respect the array structure returned by future_analysis()
  #oc_data <- array(NA, dim = c(dim(pred_data)[1], length(n_total), 2*(length(arms)-1)))
  
 
#  for(row in 1:dim(oc_data)[1])
 
  oc_data_list <- foreach( row = 1:dim(pred_data)[1],  
                      .export=c("future_analysis","pred_data","df_new_long",
                                #"arms","n_total",
                                "df","add_sim_status","add_censor_sample",
                                "censor_dist"),
                      .inorder=FALSE,
                      .packages=c("survival","magrittr")) %dopar%
  {
   future_analysis(
     pred_data[row,],newdata=df_new_long, olddata=df,
     arms=arms, n_total=n_total
   )
  }
  #stopCluster(cl)
  #oc_data <- simplify2array(oc_data_list)
  # for(row in 1:dim(oc_data)[1]){
  #   oc_data[row, ,] <-  oc_list[[row]]
  #     
  #   # future_analysis(
  #   #   pred_data[row,],newdata=df_new_long, olddata=df, 
  #   #   arms=arms, n_total=n_total
  #   # )
  # }
  # 
  
  
  if( length(arms)==3){
    z1 <- sapply(oc_data_list,  FUN = function(x){x[,1]/x[,3]}, simplify="array") #oc_data[,1,]/oc_data[,3,]
    z2 <-  sapply(oc_data_list,  FUN = function(x){x[,2]/x[,4]}, simplify="array")#oc_data[,2,]/oc_data[,4,]
  } else{
    z1 <- z2 <- sapply(oc_data_list,  FUN = function(x){x[,1]/x[,2]}, simplify="array")#oc_data[,1,]/oc_data[,2,]
  }
  
  #1-side alpah/2 test
  any_result <- (pnorm(z1)<0.025) |(pnorm(z2)<0.025)
  power <- apply(any_result,1, mean)
  #apply(pnorm(z1)<0.025,2, mean)
  #apply(pnorm(z2)<0.025,2, mean)
  return(power)
}
