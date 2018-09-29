#' Create permutation for group vector with K values
#' 
#' @param trt A binary vector indicating group, contains values from 1 to K
#' @param imputed_delta_matrix A matrix contianing the delta each obsevations will recieve for each group it belongs.
#' @return A permuted vector of trt such that there are at least two events in each group .
#' 
#' @details This is an inner function inside KONPSURV package, not for the users.
sample_trt_K_sample <- function(trt,imputed_delta_matrix)
{
  ptrt <- sample(trt)
  n <- length(trt)
  delta <- imputed_delta_matrix[cbind(1:n,ptrt)] # This is the delta vector for the specific permutation
  tf_vec <- sapply(unique(trt),function(k){sum(delta[ptrt==k])<2}) # A vector incdicating whether one of the groups has less then two events
  while(any(tf_vec)) #making sure all groups has et least two events
  {
    ptrt <- sample(trt)
    delta <- imputed_delta_matrix[cbind(1:n,ptrt)] # This is the delta vector for the specific permutation
    tf_vec <- sapply(unique(trt),function(k){sum(delta[ptrt==k])<2}) # A vector incdicating whether one of the groups has less then two events
  }
  return(ptrt)
}

#' KONP Tests for Equality of K Distributions for Right-Censored Data
#' 
#' @param time The follow-up time.
#' @param delta A binary status vector, where 0 stands for censored observation and 1 for events.
#' @param trt Group vector, must contain at least two unique values.
#' @param n.perm The number of permutations.
#' @param n.impu The number of imputations, for each imputation n.perm permutations will be executed.
#' 
#' @return Two test statistics and their respective P-values are returned: \cr 
#' 
#' \code{pv_chisq} - returns the P-value based on the Chi-square statistic. \cr 
#' \code{pv_lr} - returns the P-value based on the likelihood ratio statistic. \cr 
#' \code{chisq_test_stat} - returns the Chi-square test statistic. \cr 
#' \code{lr_test_stat} - returns the likelihood ratio test statistic.
#' @examples
#' ## Generate some data to preform the test
#' set.seed(1)
#' n <- 50
#' time <- rexp(n)
#' delta <- sample(c(0,1),n,TRUE)
#' trt <- c(rep(1,25),rep(2,25))
#' 
#' konp_test(time,delta,trt,n.perm=10^3)
#' 
#' @details The KONP test is a powerful non-parametric K-sample survival test, for equality of K distributions.
#'  It is consistent against any form of difference between the two groups,
#'  and has been shown to offer more power than alternative tests, in non-proportional hazards data.
#'  
konp_test<-function(time,delta,trt,n.perm,n.impu = 1)
{
  
  #Input testing 
  
  if (length(unique(trt)) == 2)
  {
    return(konp_2_sample(time,delta,trt,n.perm,n.impu )) #running the faster function when K=2
  }
  
  if (all(unique(delta)!= c(0,1)) & all(unique(delta)!= c(1,0)) & all(unique(delta)!= 1) & all(unique(delta)!= 0))
  {
    stop ("ERROR - delta vecotr must contain 0 or 1 only\n")
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
  
  #Making the trt vector be equal to 1 until K
  trt_ex <- rep(NA,length(trt))
  trt_unq <- unique(trt)
  K <- length(trt_unq) #Nuber of different groups
  new_group <- seq(1,K,1)
  curr_group <- 1
  for (old_group in trt_unq)
  {
    trt_ex[trt==old_group] <- curr_group
    curr_group <- curr_group+1
  }
  trt <- trt_ex
  
  
  #Making sure we have at least two events in each group
  for (group in 1:K)
  {
    if (sum(delta[trt==group])<2)
    {
      stop ("ERROR - Data must have at least two events in each groups in order to preform test\n")
    }
  }
  
  
  if (min(time)<=0)
  {
    stop ("ERROR - the time vector has negative or zero values\n")
  }
  
  n<-length(time)
  M <- Inf # represents a really large number
  
  test_prep <- sapply(1:K,function(k){ # this function will return all that is neccessery to run the test
    fit <- survival::survfit(survival::Surv(time[trt==k], delta[trt==k]) ~ 1)
    sk <- c(1,fit$surv)
    time_k <- c(0,fit$time)
    n_k=sum(trt==k)
    
    max_ev_k <- max(time[trt==k & delta==1])  #last event time in group k
    
    max_obs_k <- max(time[trt==k])  #last event\censor time in group k
    
    tau_k <- ifelse(max_ev_k==max_obs_k,M,max_ev_k) #last time to estimate km in group 0
    
    return(list(sk=sk,time_k=time_k,n_k=n_k,tau_k=tau_k))
  })
  
  tau_k <- unlist(test_prep["tau_k",])
  tau <- sort(tau_k,decreasing = T)[2]
  n_k <- unlist(test_prep["n_k",])
  
  a <- hhgsurv_test_stat_K_sample(s_group = test_prep["sk",],time_group = test_prep["time_k",],n_vec = n_k,
                                  time=time,delta=delta,trt=trt,tau_k = tau_k ,tau = tau)
  
  chisq_test_stat <- a$chisq_stat
  
  lr_test_stat <- a$lr_stat
  
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
  
  
  imputed_time_array <- array(data = NA,dim=c(n,K,n.impu)) #a time mat for each observation if the observation had changed a group
  imputed_delta_array <- array(data = NA,dim=c(n,K,n.impu)) #a delta mat for each observation if the observation had changed a group
  
  cen <- 1 - delta
  
  for (k in 1:K)
  {
    fit <- survival::survfit(survival::Surv(time[trt==k], cen[trt==k]) ~ 1)
    prob.ck <- -diff(c(1,fit$surv)) # probabilities for treatment group censorship
    time_k <- fit$time
    if (sum(prob.ck)==0) #this claim will be correct if and only if there are no censorship in the group
    {
      prob.ck <- rep(1/2,2)
      time_k <- rep(Inf,2) #giving the censorship value a value that will be bigger then the time event for sure
    }else
    {
      if(sum(prob.ck)<1) # giving the correct value of being later then the last observation time
      {
        prob.ck <- c(prob.ck,1-sum(prob.ck))
        time_k <- c(time_k,max(time_k))
      }
    }
    for (subj in (1:n)[trt!=k])
    {
      c <- sample(time_k,n.impu,prob=prob.ck,replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      t <- rep(time[subj],n.impu)
      prob <- prob.t[values.t > time[subj]]
      values <- values.t[values.t > time[subj]]
      if ((length(prob)>0) & (sum(prob)>0) & delta[subj]==0)
      { t <- sample(c(values,values), n.impu, prob = c(prob,prob),replace = T) + stats::rnorm(n = n.impu,sd = 10^-4)
      }
      imputed_delta_array[subj,k,] <- ifelse(t<=c,1,0)
      imputed_time_array[subj,k,] <-  pmax(ifelse(t<=c,t,c),10^-10) #pmax is for not getting negative values
    }
    imputed_delta_array[trt==k,k,] <- delta[trt==k] #giving the original data when the observation didn't changed group
    imputed_time_array[trt==k,k,] <- time[trt==k]   #giving the original data when the observation didn't changed group
  }
  
  chisq_stat_perm <- c()
  lr_stat_perm <- c()
  
  pv_chisq_vec <- rep(NA,n.impu)
  pv_lr_vec <- rep(NA,n.impu)

  
  for (imp in 1:n.impu)
  {
    ptrt_mat <- replicate(n.perm,sample_trt_K_sample(trt,imputed_delta_array[,,imp])) #creating the permutation for n.perm permutations in a matrix
    
    perm <- get_perm_stats_K_sample(ptrt_mat,imputed_time_array[,,imp],imputed_delta_array[,,imp],n_perm = n.perm,n_vec = n_k)
    
    pv_chisq_vec[imp] <- (sum(chisq_test_stat<=perm$chisq_stat)+1)/(n.perm +1)
    pv_lr_vec[imp] <- (sum(lr_test_stat<=perm$lr_stat)+1)/(n.perm +1)
    
    chisq_stat_perm <- c(chisq_stat_perm,perm$chisq_stat)
    lr_stat_perm <- c(lr_stat_perm,perm$lr_stat)
  }
  
  
  pv_chisq <- (sum(chisq_test_stat<=chisq_stat_perm) + 1)/(n.impu*n.perm+1)
  pv_lr <- (sum(lr_test_stat<=lr_stat_perm) + 1)/(n.impu*n.perm+1)
  
  
  return(list(pv_chisq=pv_chisq, pv_lr=pv_lr,
              chisq_test_stat=chisq_test_stat, lr_test_stat=lr_test_stat))
}
