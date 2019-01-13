#' Create permutation for group vector
#' 
#' @param trt A binary vector indicating group.
#' @param delta_orginal A binary status vector for the original data, 0 for censored observations and 1 for events.
#' @param imputed_altern_delta_vec A binary imputed status vector for observations that will change group in the permuted data.
#' @return A permuted vector of trt such that there are at least two events in each group .
#' 
#' @details This is an inner function inside KONPSURV package, not for the users.
sample_trt <- function(trt,delta_orginal,imputed_altern_delta_vec)
{
  ptrt <- sample(trt)
  delta <- ifelse(ptrt==trt,delta_orginal,imputed_altern_delta_vec)
  while(sum(delta[ptrt==0])<2 | sum(delta[ptrt==1])<2) #making sure we have at least 2 events in each group
  {
    ptrt<-sample(trt)
    delta <- ifelse(ptrt==trt,delta_orginal,imputed_altern_delta_vec)
  }
  return(ptrt)
}


#' KONP Test for Equality of Two Distributions for Right Censored Data
#' 
#' @param time The follow up time.
#' @param delta A binary status vector, where 0 stands for censored observations and 1 stands for events.
#' @param trt Group vector, must contain only two unique values.
#' @param n.perm The number of permutations.
#' @param n.impu The number of imputations, for each imputation n.perm permutations will be executed.
#' 
#' @return Two statistics and their appropriate P-values are returned: \cr 
#' 
#' \code{pv_chisq} - returns the P-value based on the Chi-square statistic. \cr 
#' \code{pv_lr} - returns the P-value based on the likelihood ratio statistic. \cr 
#' \code{pv_cauchy} - returns the P-value based on the Cauchy statistic. \cr 
#' \code{chisq_test_stat} - returns the Chi-square test statistic. \cr 
#' \code{lr_test_stat} - returns the likelihood ratio test statistic.
#' 
#' @details This is an inner function inside KONPSURV package, not for the users.
#' The user should run the K-sample test and if K=2 this function will run.
#'  
konp_2_sample<-function(time,delta,trt,n.perm,n.impu = 1)
{
  
  #Input testing 
  if (all(unique(delta)!= c(0,1)) & all(unique(delta)!= c(1,0)) & all(unique(delta)!= 1) & all(unique(delta)!= 0))
  {
    stop ("ERROR - delta vecotr must contain 0 or 1 only\n")
  }
  
  if (length(unique(trt)) != 2)
  {
    stop ("ERROR - there should be exactly 2 treatment groups\n")
  }
  
  if (class(time) != "numeric" & class(time) != "integer")
  {
    stop ("ERROR - time class sould be numeric or integer\n")
  }
  
  if (length(time) != length(trt) | length(time) != length(delta))
  {
    stop ("ERROR - Vectors time, trt and delta must be in the same length\n")
  }
  
  if (sum(is.na(time))+ sum(is.na(delta)) +sum(is.na(trt))>0)
  {
    stop ("ERROR - time or delta or trt has NA's in the vector\n")
  }
  
  if (n.perm%%1 != 0 | n.perm<1)
  {
    stop ("ERROR - n.perm must be a natural number\n")
  }
  
  if (n.impu%%1 != 0 | n.impu<1)
  {
    stop ("ERROR - n.impu must be a natural number\n")
  }
  
  #Making the trt vector be equal to 0 and 1
  trt_ex <- rep(NA,length(trt))
  trt_unq <- unique(trt)
  trt_ex[trt==trt_unq[1]] <- 0
  trt_ex[trt==trt_unq[2]] <- 1
  trt <- trt_ex
  
  if (sum(delta[trt==0])<2 | sum(delta[trt==1])<2)
  {
    stop ("ERROR - Data must have at least two events in each groups in order to preform test\n")
  }
  
  
  if (min(time)<=0)
  {
    stop ("ERROR - the time vector has negative or zero values\n")
  }
  
  
  n<-length(time)
  
  fit <- survival::survfit(survival::Surv(time[trt==0], delta[trt==0]) ~ 1)
  s0 <- c(1,fit$surv)
  time0 <- c(0,fit$time)
  
  fit <- survival::survfit(survival::Surv(time[trt==1], delta[trt==1]) ~ 1)
  s1 <- c(1,fit$surv)
  time1 <- c(0,fit$time)
  
  M <- Inf # represents a really large number
  
  max_ev_0 <- max(time[trt==0 & delta==1])  #last event time in group 0
  max_ev_1 <- max(time[trt==1 &delta==1])  
  
  max_obs_0 <- max(time[trt==0])  #last event\censor time in group 0
  max_obs_1 <- max(time[trt==1])  
  
  max_0 <- ifelse(max_ev_0==max_obs_0,M,max_ev_0) #last time to estimate km in group 0
  max_1 <- ifelse(max_ev_1==max_obs_1,M,max_ev_1) #last time to estimate km in group 1
  
  
  tau <- min(max_0,max_1) # the last time in which we can estimate km in both groups
  
  test_stat_list <- hhgsurv_test_stat(s0 = s0,s1 = s1,time0 = time0,time1 = time1,time = time,delta = delta,
                            trt = trt,tau = tau)
  
  chisq_test_stat <- test_stat_list$chisq_stat
  
  lr_test_stat <- test_stat_list$lr_stat
  
  
  #prepering for time only imputation
  fit <- survival::survfit(survival::Surv(time, delta) ~ 1)
  prob.t <- -diff(c(1,fit$surv)) 
  values.t <- fit$time
  if(sum(prob.t)<1)
  {
    prob.t <- c(prob.t,1-sum(prob.t))
    values.t <- c(values.t,max(values.t)+1) #each observation that will get max(values.t)+1 will be censored for sure
  }
  
  
  #prepearing for imputation for censoring (and some of the time)
  
  cen <- 1 - delta
  fit <- survival::survfit(survival::Surv(time[trt==1], cen[trt==1]) ~ 1)
  prob.c1 <- -diff(c(1,fit$surv)) # probabilities for treatment group censorship
  time1 <- fit$time
  if (sum(prob.c1)==0) #this claim will be correct if and only if there are no censorship in the group
  {
    prob.c1 <- rep(1/2,2)
    time1 <- rep(Inf,2) #giving the censorship value a value that will be bigger then the time event for sure
  }else
  {
    if(sum(prob.c1)<1) # giving the correct value of being later then the last observation time
    {
      prob.c1 <- c(prob.c1,1-sum(prob.c1))
      time1 <- c(time1,max(time1))
    }
  }
  fit <- survival::survfit(survival::Surv(time[trt==0], cen[trt==0]) ~ 1)
  prob.c0 <- -diff(c(1,fit$surv)) # probabilities for control group censorship
  time0 <-fit$time
  if (sum(prob.c0)==0) #this claim will be correct if and only if there are no censorship in the group
  {
    prob.c0 <- rep(1/2,2)
    time0 <-  rep(Inf,2) #giving the censorship value a value that will be bigger then the time event for sure
  }else
  {
    if(sum(prob.c0)<1) # giving the correct value of being later then the last observation time
    {
      prob.c0 <- c(prob.c0,1-sum(prob.c0))
      time0 <- c(time0,max(time0))
    }
  }
  
  imputed_altern_time_mat<-matrix(data = NA,ncol = n.impu,nrow = n,byrow = F) #a time mat for each observation if the observation had changed a group
  imputed_altern_delta_mat<-matrix(data = NA,ncol = n.impu,nrow = n,byrow = F) #a delta mat for each observation if the observation had changed a group
  
  
  for (subj in 1:n) #running for trt==0
  {
    if (trt[subj]==0)
    {
      c <- sample(time1,n.impu,prob=prob.c1,replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      t <- rep(time[subj],n.impu)
      prob <- prob.t[values.t > time[subj]]
      values <- values.t[values.t > time[subj]]
      if ((length(prob)>0) & (sum(prob)>0) & delta[subj]==0)
      { t <- sample(c(values,values), n.impu, prob = c(prob,prob),replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      }
      imputed_altern_delta_mat[subj,] <- ifelse(t<=c,1,0)
      imputed_altern_time_mat[subj,] <- pmax(ifelse(t<=c,t,c),10^-10) #pmax is for not getting negative values
    }
    if (trt[subj]==1)
    {
      c <- sample(time0,n.impu,prob=prob.c0,replace = T) +stats::rnorm(n = n.impu,sd = 10^-4)
      t <-  rep(time[subj],n.impu)
      prob <- prob.t[values.t > time[subj]]
      values <- values.t[values.t > time[subj]]
      if ((length(prob)>0) & (sum(prob)>0) & delta[subj]==0)
      { t <-  sample(c(values,values), n.impu, prob = c(prob,prob),replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      }
      imputed_altern_delta_mat[subj,] <- ifelse(t<=c,1,0)
      imputed_altern_time_mat[subj,] <- pmax(ifelse(t<=c,t,c),10^-10) #pmax is for not getting negative values
    }
  }
  
  
  chisq_stat_perm <- c()
  lr_stat_perm <- c()
  
  pv_chisq_vec <- rep(NA,n.impu)
  pv_lr_vec <- rep(NA,n.impu)

  
  for (imp in 1:n.impu)
  {
    ptrt_mat <- replicate(n.perm,sample_trt(trt,delta,imputed_altern_delta_mat[,imp])) #creating the permutation for n.perm permutations in a matrix
    
    perm <- get_perm_stats(trt,ptrt_mat,time,delta,imputed_altern_time_mat[,imp],
                          imputed_altern_delta_mat[,imp],n_perm = n.perm)
    
    pv_chisq_vec[imp] <- (sum(chisq_test_stat<=perm$chisq_stat)+1)/(n.perm +1)
    pv_lr_vec[imp] <- (sum(lr_test_stat<=perm$lr_stat)+1)/(n.perm +1)
    

    chisq_stat_perm <- c(chisq_stat_perm,perm$chisq_stat)
    lr_stat_perm <- c(lr_stat_perm,perm$lr_stat)
    
  }
  
  pv_chisq <- (sum(chisq_test_stat<=chisq_stat_perm) + 1)/(n.impu*n.perm+1)
  pv_lr <- (sum(lr_test_stat<=lr_stat_perm) + 1)/(n.impu*n.perm+1)
  
  #Calculate Cauchy test statistic:
  #First calculate logrank P-value
  fit_lr <- survival::survdiff(survival::Surv(time, delta) ~ trt , rho=0)
  Pvalue_logrank <- 1 - stats::pchisq(fit_lr$chisq, 1)
  cau <- mean(c(tan((0.5-pv_chisq)*pi),
                tan((0.5-pv_lr)*pi), tan((0.5-Pvalue_logrank)*pi)))
  pv_cauchy <- 0.5-atan(cau)/pi
  
  return(list(pv_chisq=pv_chisq, pv_lr=pv_lr, pv_cauchy = pv_cauchy,
              chisq_test_stat=chisq_test_stat, lr_test_stat=lr_test_stat))
}




# konp 2 sample with imputation output ------------------------------------

#' KONP Test for Equality of Two Distributions for Right Censored Data
#' 
#' @param time The follow up time.
#' @param delta A binary status vector, where 0 stands for censored observations and 1 stands for events.
#' @param trt Group vector, must contain only two unique values.
#' @param n.perm The number of permutations.
#' @param n.impu The number of imputations, for each imputation n.perm permutations will be executed.
#' 
#' @return Two statistics and their appropriate P-values are returned: \cr 
#' 
#' \code{pv_chisq} - returns the P-value based on the Chi-square statistic. \cr 
#' \code{pv_lr} - returns the P-value based on the likelihood ratio statistic. \cr 
#' \code{pv_cauchy} - returns the P-value based on the Cauchy statistic. \cr 
#' \code{chisq_test_stat} - returns the Chi-square test statistic. \cr 
#' \code{lr_test_stat} - returns the likelihood ratio test statistic.
#' 
#' @details This is inner funtion inside KONP package made only for running simulation and not for users
konp_2_sample_impu_output<-function(time,delta,trt,n.perm,n.impu = 1)
{
  
  #Input testing 
  if (all(unique(delta)!= c(0,1)) & all(unique(delta)!= c(1,0)) & all(unique(delta)!= 1) & all(unique(delta)!= 0))
  {
    stop ("ERROR - delta vecotr must contain 0 or 1 only\n")
  }
  
  if (length(unique(trt)) != 2)
  {
    stop ("ERROR - there should be exactly 2 treatment groups\n")
  }
  
  if (class(time) != "numeric" & class(time) != "integer")
  {
    stop ("ERROR - time class sould be numeric or integer\n")
  }
  
  if (length(time) != length(trt) | length(time) != length(delta))
  {
    stop ("ERROR - Vectors time, trt and delta must be in the same length\n")
  }
  
  if (sum(is.na(time))+ sum(is.na(delta)) +sum(is.na(trt))>0)
  {
    stop ("ERROR - time or delta or trt has NA's in the vector\n")
  }
  
  if (n.perm%%1 != 0 | n.perm<1)
  {
    stop ("ERROR - n.perm must be a natural number\n")
  }
  
  if (n.impu%%1 != 0 | n.impu<1)
  {
    stop ("ERROR - n.impu must be a natural number\n")
  }
  
  #Making the trt vector be equal to 0 and 1
  trt_ex <- rep(NA,length(trt))
  trt_unq <- unique(trt)
  trt_ex[trt==trt_unq[1]] <- 0
  trt_ex[trt==trt_unq[2]] <- 1
  trt <- trt_ex
  
  if (sum(delta[trt==0])<2 | sum(delta[trt==1])<2)
  {
    stop ("ERROR - Data must have at least two events in each groups in order to preform test\n")
  }
  
  
  if (min(time)<=0)
  {
    stop ("ERROR - the time vector has negative or zero values\n")
  }
  
  
  n<-length(time)
  
  fit <- survival::survfit(survival::Surv(time[trt==0], delta[trt==0]) ~ 1)
  s0 <- c(1,fit$surv)
  time0 <- c(0,fit$time)
  
  fit <- survival::survfit(survival::Surv(time[trt==1], delta[trt==1]) ~ 1)
  s1 <- c(1,fit$surv)
  time1 <- c(0,fit$time)
  
  M <- Inf # represents a really large number
  
  max_ev_0 <- max(time[trt==0 & delta==1])  #last event time in group 0
  max_ev_1 <- max(time[trt==1 &delta==1])  
  
  max_obs_0 <- max(time[trt==0])  #last event\censor time in group 0
  max_obs_1 <- max(time[trt==1])  
  
  max_0 <- ifelse(max_ev_0==max_obs_0,M,max_ev_0) #last time to estimate km in group 0
  max_1 <- ifelse(max_ev_1==max_obs_1,M,max_ev_1) #last time to estimate km in group 1
  
  
  tau <- min(max_0,max_1) # the last time in which we can estimate km in both groups
  
  test_stat_list <- hhgsurv_test_stat(s0 = s0,s1 = s1,time0 = time0,time1 = time1,time = time,delta = delta,
                                      trt = trt,tau = tau)
  
  chisq_test_stat <- test_stat_list$chisq_stat
  
  lr_test_stat <- test_stat_list$lr_stat
  
  tab_usage <- test_stat_list$tab_usage
  
  
  #prepering for time only imputation
  fit <- survival::survfit(survival::Surv(time, delta) ~ 1)
  prob.t <- -diff(c(1,fit$surv)) 
  values.t <- fit$time
  if(sum(prob.t)<1)
  {
    prob.t <- c(prob.t,1-sum(prob.t))
    values.t <- c(values.t,max(values.t)+1) #each observation that will get max(values.t)+1 will be censored for sure
  }
  
  
  #prepearing for imputation for censoring (and some of the time)
  
  cen <- 1 - delta
  fit <- survival::survfit(survival::Surv(time[trt==1], cen[trt==1]) ~ 1)
  prob.c1 <- -diff(c(1,fit$surv)) # probabilities for treatment group censorship
  time1 <- fit$time
  if (sum(prob.c1)==0) #this claim will be correct if and only if there are no censorship in the group
  {
    prob.c1 <- rep(1/2,2)
    time1 <- rep(Inf,2) #giving the censorship value a value that will be bigger then the time event for sure
  }else
  {
    if(sum(prob.c1)<1) # giving the correct value of being later then the last observation time
    {
      prob.c1 <- c(prob.c1,1-sum(prob.c1))
      time1 <- c(time1,max(time1))
    }
  }
  fit <- survival::survfit(survival::Surv(time[trt==0], cen[trt==0]) ~ 1)
  prob.c0 <- -diff(c(1,fit$surv)) # probabilities for control group censorship
  time0 <-fit$time
  if (sum(prob.c0)==0) #this claim will be correct if and only if there are no censorship in the group
  {
    prob.c0 <- rep(1/2,2)
    time0 <-  rep(Inf,2) #giving the censorship value a value that will be bigger then the time event for sure
  }else
  {
    if(sum(prob.c0)<1) # giving the correct value of being later then the last observation time
    {
      prob.c0 <- c(prob.c0,1-sum(prob.c0))
      time0 <- c(time0,max(time0))
    }
  }
  
  
  
  imputed_altern_time_mat<-matrix(data = NA,ncol = n.impu,nrow = n,byrow = F) #a time mat for each observation if the observation had changed a group
  imputed_altern_delta_mat<-matrix(data = NA,ncol = n.impu,nrow = n,byrow = F) #a delta mat for each observation if the observation had changed a group
  
  
  for (subj in 1:n) #running for trt==0
  {
    if (trt[subj]==0)
    {
      c <- sample(time1,n.impu,prob=prob.c1,replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      t <- rep(time[subj],n.impu)
      prob <- prob.t[values.t > time[subj]]
      values <- values.t[values.t > time[subj]]
      if ((length(prob)>0) & (sum(prob)>0) & delta[subj]==0)
      { t <- sample(c(values,values), n.impu, prob = c(prob,prob),replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      }
      imputed_altern_delta_mat[subj,] <- ifelse(t<=c,1,0)
      imputed_altern_time_mat[subj,] <- pmax(ifelse(t<=c,t,c),10^-10) #pmax is for not getting negative values
    }
    if (trt[subj]==1)
    {
      c <- sample(time0,n.impu,prob=prob.c0,replace = T) +stats::rnorm(n = n.impu,sd = 10^-4)
      t <-  rep(time[subj],n.impu)
      prob <- prob.t[values.t > time[subj]]
      values <- values.t[values.t > time[subj]]
      if ((length(prob)>0) & (sum(prob)>0) & delta[subj]==0)
      { t <-  sample(c(values,values), n.impu, prob = c(prob,prob),replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      }
      imputed_altern_delta_mat[subj,] <- ifelse(t<=c,1,0)
      imputed_altern_time_mat[subj,] <- pmax(ifelse(t<=c,t,c),10^-10) #pmax is for not getting negative values
    }
  }
  
  chisq_stat_perm <- c()
  lr_stat_perm <- c()
  
  pv_chisq_vec <- rep(NA,n.impu)
  pv_lr_vec <- rep(NA,n.impu)
  
  tab_usage_perm_vec <- rep(NA,n.impu)
  
  
  for (imp in 1:n.impu)
  {
    ptrt_mat <- replicate(n.perm,sample_trt(trt,delta,imputed_altern_delta_mat[,imp])) #creating the permutation for n.perm permutations in a matrix
    
    perm <- get_perm_stats(trt,ptrt_mat,time,delta,imputed_altern_time_mat[,imp],
                           imputed_altern_delta_mat[,imp],n_perm = n.perm)
    
    pv_chisq_vec[imp] <- (sum(chisq_test_stat<=perm$chisq_stat)+1)/(n.perm +1)
    pv_lr_vec[imp] <- (sum(lr_test_stat<=perm$lr_stat)+1)/(n.perm +1)
    
    tab_usage_perm_vec[imp] <- mean(perm$tab_usage_perm)
    
    chisq_stat_perm <- c(chisq_stat_perm,perm$chisq_stat)
    lr_stat_perm <- c(lr_stat_perm,perm$lr_stat)
    
  }
  
  
  pv_chisq <- (sum(chisq_test_stat<=chisq_stat_perm) + 1)/(n.impu*n.perm+1)
  pv_lr <- (sum(lr_test_stat<=lr_stat_perm) + 1)/(n.impu*n.perm+1)
  
  #Calculate Cauchy test statistic:
  #First calculate logrank P-value
  fit_lr <- survival::survdiff(survival::Surv(time, delta) ~ trt , rho=0)
  Pvalue_logrank <- 1 - stats::pchisq(fit_lr$chisq, 1)
  cau <- mean(c(tan((0.5-pv_chisq)*pi),
                tan((0.5-pv_lr)*pi), tan((0.5-Pvalue_logrank)*pi)))
  pv_cauchy <- 0.5-atan(cau)/pi
  
  tab_usage_perm <- mean(tab_usage_perm_vec)
  
  
  return(list(pv_chisq=pv_chisq, pv_lr=pv_lr, pv_cauchy=pv_cauchy,
              chisq_test_stat=chisq_test_stat, lr_test_stat=lr_test_stat,
              tab_usage = tab_usage, tab_usage_perm = tab_usage_perm,
              imputed_altern_delta_mat=imputed_altern_delta_mat,
              imputed_altern_time_mat=imputed_altern_time_mat,
              ptrt_mat=ptrt_mat))
}
