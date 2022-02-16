# -------------------------------------------------------------------------------------------
##### This code can be used to implement the survival-CRM,                             ##### 
##### the informative survival-CRM and benchmark for Phase 1                           ##### 
##### dose-finding trials with right censored time-to-event endpoints;                 #####  
##### and to reproduce the result in the paper ``Dose-finding                          ##### 
##### design and benchmark for a right censored endpoint''                             ##### 
#####                                                                                  ##### 
##### by Anais Andrillon, Sylvie Chevret, Shing M Lee and Lucie Biard (2020)           ##### 
# -------------------------------------------------------------------------------------------

## Load the required packages 
install.packages("pacman")
pacman::p_load(dplyr, dfcrm, survival)

# -------------------------------------------------------------------------------------------
#' getprior_exp 
#' Returns a vector of initial guesses of toxicity probabilities associated to the doses 
#' for the exponential working model for a given set of indifference intervals
#' 
#' @param halfwidth      The desired halfwidth of the indifference intervals
#' @param target         The target toxicity probability of DLT
#' @param nu             The prior guess of MTD
#' @param nlevel         The number of test doses 
#' @param tstar          The end of toxicity observation window t*
#' @return               A vector of initial guesses of toxicity probabilities 
#' 
getprior_exp <- function (halfwidth, target, nu, nlevel, tstar) {
  dosescaled <- prior <- rep(NA, nlevel)
  b <- rep(NA, nlevel + 1)
  b[1] <- -Inf
  b[(nlevel + 1)] <- Inf
  
  dosescaled[nu] <- log(-log(1-target)/tstar)
  for (k in nu:2) {
    b[k] <- log(log(-log(1-(target + halfwidth))/tstar)/ dosescaled[k])      
    if (nu > 1) {
      dosescaled[k - 1] <- log(-log(1-(target - halfwidth))/tstar)/ exp(b[k])
    }
  }
  if (nu < nlevel) {
    for (k in nu:(nlevel - 1)) {
      b[k + 1] <- log(log(-log(1-(target - halfwidth))/tstar)/ dosescaled[k])  
      dosescaled[k + 1] <- log(-log(1-(target + halfwidth))/tstar)/ exp(b[k+1])
    }
  }
  return(1-exp(-exp(dosescaled)*tstar))
  
}


# -------------------------------------------------------------------------------------------
#' find_param
#' From True proba of DLT, find parameters of time-to-event distribution 
#' 
#' @param p             The true probability of DLT 
#' @param tstar         The end of toxicity observation window t*
#' @param hazard        The time-to-toxicity hazard. "decreasing" for early toxicities; "increasing" for late toxicities or "constant"
#' @return              Distribution parameters of the time-to-event distribution 
#' 
find_param <- function(p, tstar, hazard){
  if(hazard!="constant"){
    a <- ifelse(hazard=="increasing",3,0.5)
    b = exp(1/a * log(-tstar^a/log(1-p)))
    param <- cbind(rep(a, length(p)), b)
  } else {param <- (- log(1-p)/tstar)}#constant time-to-toxicity hazard
  return(param)
}



# -------------------------------------------------------------------------------------------
#' find_params
#' From True proba of DLT and discontinuation, find parameters of time-to-event distributions
#' 
#' @param p1           The true probability of DLT 
#' @param p2           The true probability of discontinuation 
#' @return             Distribution parameters of the time-to-event distribution 
#' 

find_params <- function(p1, p2,tstar){
  l1 <- -log(1-p1-p2)/ ((p2/p1+1)*tstar)
  l2 <- -log(1-p1-p2)-l1#old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2
  #old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2
  #old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2old (p1*l1)/p2
  return(list(l1=l1,
              l2=l2))
}


 
# -------------------------------------------------------------------------------------------
#' getdataset 
#' Generate complete dataset for a trial
#'
#' @param n             The sample size of the trial
#' @param scenario      The vector of the true toxicity probabilites associated with the doses
#' @param tstar         The end of toxicity observation window t*
#' @param hazard        The time-to-toxicity hazard. "decreasing" for early toxicities; "increasing" for late toxicities or "constant"
#' @return              A complete dataset about n patients for all dose levels
#' 
getdataset <- function(n, scenario, tstar, seed, hazard) {
  set.seed(seed)
  K <- length(scenario)
  Dataset <- cbind(Patient=rep(1:n, each=K),Ui= rep(runif(n), each=K),Time=rep(NA,n*K),Event=rep(NA,n*K),Cens_time=rep(NA,n*K),Dose=rep(1:K,n))  
  i = 1
  while(i<=n){
    if(hazard=="constant"){
      Dataset[(K*(i-1)+1):(K*i),'Time'] <-  qexp(Dataset[(K*(i-1)+1):(K*i),'Ui'], rate = find_param(scenario,tstar,hazard))
    } else{
      param <- find_param(scenario,tstar,hazard)
      Dataset[(K*(i-1)+1):(K*i),'Time'] <-  qweibull(Dataset[(K*(i-1)+1):(K*i),'Ui'],param[,1],param[,2])
    }
    Dataset[(K*(i-1)+1):(K*i),'Event'] <-  ifelse(Dataset[(K*(i-1)+1):(K*i),'Time'] >= tstar,0,1) #Event: censor indicator. 1 if DLT before t*
    Dataset[(K*(i-1)+1):(K*i),'Cens_time'] <- pmin( Dataset[(K*(i-1)+1):(K*i),'Time'] ,tstar )   #Cens_time=1 if no event,  Cens_time=Time if DLT before t*
    
    i=i+1
  }
  return(Dataset)
  
}

# -------------------------------------------------------------------------------------------
#' getdataset 
#' Generate complete dataset for a trial in a setting of competing risks 
#'
#' @param n             The sample size of the trial
#' @param scenario1     The vector of the true toxicity probabilites associated with the doses
#' @param scenario2     The vector of the true discontinuation probabilites associated with the doses
#' @param tstar         The end of toxicity observation window t*
#' @return              A complete dataset about n patients for all dose levels
#' 


getdataset_CR <- function(n, scenario1, scenario2, tstar, seed) {
  set.seed(seed)
  K <- length(scenario1)
  Dataset <- cbind(Patient=rep(1:n, each=K),Ui= rep(runif(n), each=K),Time=rep(NA,n*K),Event=rep(0,n*K) ,Cens_time=rep(NA,n*K),Dose=rep(1:K,n),I_Tox=rep(NA,n*K) )  
  i = 1
  while(i<=n){
    params <- find_params(p1=scenario1, p2=scenario2,tstar=tstar)   #Instantaneous hazard back-computed from scenarios cumulative incidences 
    Dataset[(K*(i-1)+1):(K*i),'Time'] <-  qexp(Dataset[(K*(i-1)+1):(K*i),'Ui'], rate = params$l1+params$l2 ) #Time to any event sampled from exponential distribution with rate (l1+l2)
    d=1
    while(d<=K){
      Dataset[((K*(i-1)+1):(K*i))[d],'I_Tox'] <- rbinom(1,1,params$l1[d]/(params$l1[d]+params$l2[d])) #Event case determined by a random drawn from a Bernouilli(l1/(l1+l2)) for toxicity
      d <- d+1  
    }
    Dataset[(K*(i-1)+1):(K*i),'Event'][Dataset[(K*(i-1)+1):(K*i),'Time'] <= tstar &  Dataset[(K*(i-1)+1):(K*i),'I_Tox']==1] <- 1   #Toxicity
    Dataset[(K*(i-1)+1):(K*i),'Event'][Dataset[(K*(i-1)+1):(K*i),'Time'] <= tstar &  Dataset[(K*(i-1)+1):(K*i),'I_Tox']==0] <- 2   #Discontinuation
    #Administrative censoring at tstar
    Dataset[(K*(i-1)+1):(K*i),'Cens_time'] <- pmin( Dataset[(K*(i-1)+1):(K*i),'Time'] ,tstar ) #Times=1 if no event,  Times=Time if DLT before t*
    
    i=i+1
  }
  return(Dataset)
  
}




# -------------------------------------------------------------------------------------------
#' likelihood_tox_exp 
#' Survival likelihood for toxicity
#' 
#' @param beta         The unknown model parameter
#' @param event        The binary vector of toxicity occurrence
#' @param dose_level   The dose levels allocated to patients
#' @param xref         The set of numerical labels for the doses 
#' @param time         The right-censored time-to-toxicity
#' @return             The likelihood value
#' 
likelihood_tox_exp <- function(beta, event, dose_level, xref,time){
  n <- length(event)
  res <- 1
  for (i in 1:n){
    res <- res * (( (exp(exp(beta)*xref[dose_level[i]])) * (exp(-exp(xref[dose_level[i]]*exp(beta))*time[i])) )^(I(event[i]==1)*1) * (exp(-exp(xref[dose_level[i]]*exp(beta))*time[i]))^(1-I(event[i]==1)*1) )
  }
  return(res)
}

denom_tox <- function(beta, event, dose_level, xref,time){ likelihood_tox_exp(beta, event, dose_level, xref,time)*dnorm(beta, mean=0, sd=sqrt(1.34)) }
num_tox <- function(beta, event, dose_level, xref,time){beta*denom_tox(beta, event, dose_level, xref,time)}



# -------------------------------------------------------------------------------------------
#' likelihood_disc_exp 
#' Survival likelihood for discontinuation
#' 
#' @param beta         The unknown model parameter
#' @param event        The binary vector of toxicity occurrence
#' @param dose_level   The dose levels allocated to patients
#' @param xref        The set of numerical labels for the doses 
#' @param time         The right-censored time-to-toxicity
#' @return             The likelihood value
#' 
likelihood_disc_exp <- function(beta, event, dose_level, xref,time){
  res <- 1
  for (i in 1:length(event)){
    res <- res * (( (exp(-exp(beta)*xref[dose_level[i]])) * (exp(-exp(-xref[dose_level[i]]*exp(beta))*time[i])) )^(I(event[i]==2)*1) * (exp(-exp(-xref[dose_level[i]]*exp(beta))*time[i]))^(1-I(event[i]==2)*1) )
  }
  return(res)
}

denom_disc <- function(beta, event, dose_level, xref,time){ likelihood_disc_exp(beta, event, dose_level, xref,time)*dnorm(beta, mean=0, sd=sqrt(1.34)) }
num_disc <- function(beta, event, dose_level, xref,time){beta*denom_disc(beta, event, dose_level, xref,time)}


# -------------------------------------------------------------------------------------------
#' survcrm_nextdose  
#' Survival-CRM Dose finding algorithm
#' 
#' @param event        The binary vector of toxicity occurrence
#' @param dose_level   The dose levels allocated to patients
#' @param xref         The set of numerical labels for the doses 
#' @param time         The right-censored time-to-toxicity
#' @param tstar         The end of toxicity observation window t*
#' @return             The likelihood value
#' 
  
survcrm_nextdose <- function(event, dose_level, xref,time,tstar,target){
  #Model parameter estimation
  constante <- integrate(denom_tox , -Inf, Inf, event=event, dose_level=dose_level, xref=xref,time=time,abs.tol = 0)$value
  beta_hat <-   integrate(num_tox, -Inf, Inf,  event=event, dose_level=dose_level, xref=xref,time=time,abs.tol = 0)$value / constante
  
  lambda <- exp(exp(beta_hat)*xref)            #Updated model parameter
  CI_temp <- 1-exp(- lambda * tstar)           #Cumulative incidence at t* updated

  #Next dose finding algorithm
  Diff_pestim_target <- abs(CI_temp - target)
  Doses_min <- which(Diff_pestim_target==min(Diff_pestim_target)  )
  next_dose <- ifelse(length(Doses_min)==1,Doses_min,sample(Doses_min,1)) 
  
  return(list(next_dose=next_dose,
              constante=constante,
              CI_temp=CI_temp))
}

# -------------------------------------------------------------------------------------------
#' isurvcrm_nextdose  
#' iSurvival-CRM Dose finding algorithm
#' 
#' @param event        The binary vector of toxicity occurrence
#' @param dose_level   The dose levels allocated to patients
#' @param xref1        The scaled doses for toxicity
#' @param xref2        The scaled doses for progression
#' @param time         The right-censored time-to-toxicity
#' @param tstar        The end of toxicity observation window t*
#' @return             The likelihood value
#' 

isurvcrm_nextdose <- function(event, dose_level, xref1,xref2,time,tstar,target){
  #Model parameter estimation
  constante1 <- integrate(denom_tox , -Inf, Inf, event=event, dose_level=dose_level, xref=xref1,time=time,abs.tol = 0)$value
  beta_hat <-   integrate(num_tox, -Inf, Inf,  event=event, dose_level=dose_level, xref=xref1,time=time,abs.tol = 0)$value / constante1
  
  #prog
  constante2 <- integrate(denom_disc , -Inf, Inf, event=event, dose_level=dose_level, xref=xref2,time=time,abs.tol = 0)$value
  beta_hat2 <-   integrate(num_disc,-Inf, Inf,  event=event, dose_level=dose_level, xref=xref2,time=time,abs.tol = 0)$value / constante2
  
  #Updated model parameters
  lambda1 <- exp(exp(beta_hat)*xref1) 
  lambda2 <-  exp(-exp(beta_hat2)*xref2) 
  
  #Cumulative incidences at t* updated
  CI_temp1 <- lambda1*(1-exp(-(lambda1+lambda2)*tstar))/(lambda1+lambda2) 
  CI_temp2 <- lambda2*(1-exp(-(lambda1+lambda2)*tstar))/(lambda1+lambda2)
  
  
  #Next dose finding algorithm
  diff <- abs(CI_temp1 - target)
  Doses_min <- which(diff==min(diff)  )
  next_dose <- ifelse(length(Doses_min)==1,Doses_min,sample(Doses_min,1)) 
  
  
  return(list(next_dose=next_dose,
              constante1=constante1,
              CI_temp1=CI_temp1,
              CI_temp2=CI_temp2))
}


# -------------------------------------------------------------------------------------------
#' titecrm_nextdose 
#' TITE-CRM Dose finding algorithm
#' 
#' @param event        The binary vector of toxicity occurrence
#' @param dose_level   The dose levels allocated to patients
#' @param xref         The set of numerical labels for the doses 
#' @param time         The right-censored time-to-toxicity
#' @param tstar         The end of toxicity observation window t*
#' @return             The likelihood value
#' 

titecrm_nextdose <- function(xref, Tox, weights,skeleton){
  #Model parameter estimation
  constante <- integrate(crmh, -Inf, Inf, xref, Tox, weights, sqrt(1.34), abs.tol = 0)[[1]]
  est <- integrate(crmht, -Inf, Inf, xref, Tox, weights, sqrt(1.34), abs.tol = 0)[[1]]/constante
  
  CI_temp <- skeleton^exp(est)                  #Probability of toxicity updated 
  
  #Next dose finding algorithm
  Diff_pestim_target <- abs(CI_temp - target)
  Doses_min <- which(Diff_pestim_target==min(Diff_pestim_target)  )
  next_dose <- ifelse(length(Doses_min)==1,Doses_min,sample(Doses_min,1)) 
  
  return(list(next_dose=next_dose,
              constante=constante,
              CI_temp=CI_temp))
}


# ----------------------------------------------------------------------------------------------
#' survcrm_sim :       
#' Generate simulation replicates of phase I trial using the proposal survival-CRM
#' 
#' @param N            Number of simulation  
#' @param scenario     Vector of the true toxicity probabilites associated with the doses
#' @param n            Sample size of the trial
#' @param target       Target toxicity probability of DLT
#' @param tstar        End of the toxicity observation window t*
#' @param dose_init    Initial dose d1
#' @param rate         Expected number of arrivals per observation window
#' @param accrual      Inter-patient arrival scheme: "fixed" or "poisson"
#' @param hazard       The time-to-toxicity hazard. "decreasing" for early toxicities; "increasing" for late toxicities or "constant"
#' @param nu           The prior guess of MTD
#' @param halfwidth    The desired halfwidth of the indifference intervals
#' @return             Simulation results of the survival-CRM
#' 
survcrm_sim <- function(N,scenario,n,target,tstar,dose_init,rate,accrual,hazard,nu,halfwidth){
  K <- length(scenario)
  diff <- as.numeric(as.character(abs(scenario-target))) 
  correct <- which(diff== min(diff))
  toxic <- which(scenario > scenario[correct])

  #Initial guesses of toxicity probabilities associated the doses
  skeleton <- getprior_exp(halfwidth=halfwidth, target, nu=nu, nlevel=K, tstar=tstar) 
  xref <-  log(-log(1-skeleton)/tstar)              #Set of numerical labels for the doses investigated in the trial
  b1 <- log(log(-log(1-target)/tstar)/ xref[dose_init]) #Beta value such that the probability of toxicity at dose level 1 is equal to the target


  CI_estim <- matrix(nrow=N,ncol=length(scenario))  #Cumulative incidence estimate at t* (mean over N simulations) 
  Dose_select <- rep(0,length(scenario))            #No time each dose is selected over N simulations
  OD <- rep(NA,N)                                   #No patients treated at a toxic dose during a trial
  No_MTD <- rep(NA,N)                               #No patients treated at the MTD during a trial 
  No_DLT <- rep(NA,N)                               #No patients who experienced a DLT during a trial
  No_Stop <- rep(NA,N)                              #No trials stopped for safety decision
  
  s <- 1
  while (s<=N) {
    Dataset <- getdataset(n, scenario,  tstar, seed=s ,hazard)  #Complete data generation for n patients
    if (accrual == "fixed") {                       #Patient accrual scheme
      next_arrival <- tstar/rate
    }   else if (accrual == "poisson") {
      next_arrival <- rexp(1, rate/tstar)
    }
    arrival <- u <- NULL                            #Arrival of the first patient
    i <- 1
    Dataset_temp <- Dataset[Dataset[,'Patient']==1 & Dataset[,'Dose']==dose_init,c("Patient","Time", "Event", "Cens_time", "Dose"), drop = FALSE]
    
    while(i <=(n-1)){
      arrival <- c(arrival, next_arrival)
      if (i > 1){ 
        Dataset_temp <- rbind(Dataset_temp,Dataset[Dataset[,'Patient']==i & Dataset[,'Dose']==next_dose , c("Patient","Time", "Event", "Cens_time", "Dose")]) #Sequential data built from the complete Dataset
      }
      if(Dataset_temp[i,"Time"]<=tstar){            #If patient i experienced a toxicity
        unew <- Dataset_temp[i,"Time"]    
      } else  {                                     #Patient i has no toxicity
        unew <- Inf}
      u <- c(u,unew)
      uvent <- u+ arrival 
      
      if (accrual == "fixed") {          
        next_arrival <- next_arrival + tstar/rate
      }    else if (accrual == "poisson") {
        next_arrival <- next_arrival + rexp(1, rate/tstar)
      } 
      
      Tox <- rep(0, length(Dataset_temp[,"Event"])) #Vector of toxicity occurrence
      Tox[uvent <= next_arrival] <- 1 
      
      Cens <- pmin(next_arrival, uvent) - arrival    #Censor indicator
      Time <- pmin(Cens, tstar)                     #Censored time-to-toxicity
      
      # dose-finding algorithm   
      res <- survcrm_nextdose(event=Tox, dose_level=Dataset_temp[,'Dose'], xref=xref,time=Time,tstar,target)
      
      next_dose <- res$next_dose
      #Avoid skipping doses in escalation 
      next_dose <- min(next_dose, (Dataset_temp[i,'Dose'] + 1)) 
      
      # safety stopping rule
      if (integrate(denom_tox, -Inf, b1 ,event=Tox, dose_level=Dataset_temp[,'Dose'], xref=xref,time=Time,abs.tol = 0)$value / res$constante >=0.95) {
        break
      }
       i <- i+1
       
    } 
    
    #Completion of the trial with the last patient 
    if(i==n){
      Dataset_temp <- rbind(Dataset_temp,Dataset[Dataset[,'Patient']==i & Dataset[,'Dose']==next_dose , c("Patient","Time", "Event", "Cens_time", "Dose")])
      Tox <- Dataset_temp[,"Event"]
      Time <- Dataset_temp[,"Cens_time"]
      
      res <- survcrm_nextdose(event=Tox, dose_level=Dataset_temp[,'Dose'], xref=xref,time=Time,tstar,target)
      next_dose <- res$next_dose      
      
      CI_estim[s,] <- res$CI_temp                    #Cumulative incidence estimate at t* for trial s 
      Dose_select[next_dose] <- Dose_select[next_dose]+1  
      No_MTD[s] <- sum(Dataset_temp[,"Dose"] %in% correct)
      OD[s] <- sum(Dataset_temp[,"Dose"] %in% toxic)
      No_DLT[s] <- sum(Dataset_temp[,"Event"])
    }
    else(No_Stop[s] <- i)
    
    s <- s+1
  }
  No_trials <- sum(is.na(No_Stop))                    #Number of terminated trials 

  PS <- (Dose_select/No_trials)*100                   #Percent of selection of each dose over No_trials trials
  
  
  return(list(Scenario=round(scenario,2),
              PS=round(PS,0),
              PCS = round(sum(PS[correct]),0),
              Accuracy_index = round(1-length(scenario)*(sum(abs(scenario-target)*(Dose_select/No_trials))/sum(abs(scenario-target))),2), 
              POS = round(sum(PS[toxic]), 0),
              R_Bais_MTD = round(mean((colMeans(CI_estim,na.rm = TRUE)[correct]-scenario[correct])/scenario[correct] ) ,3), 
              OD=round(mean(OD,na.rm = TRUE),2),
              No_DLT = round(mean(No_DLT,na.rm = TRUE),2),
              No_MTD=round(mean(No_MTD,na.rm = TRUE),2),
              P_stop = (N-No_trials)/N*100
  ))
  
}     



# ----------------------------------------------------------------------------------------------
#' isurvcrm_sim :       
#' Generate simulation replicates of phase I trial using the proposal iSurvival-CRM
#' 
#' @param N            Number of simulation  
#' @param scenario1    Vector of the true toxicity probabilites associated with the doses
#' @param scenario2    Vector of the true discontinuation probabilites associated with the doses
#' @param n            Sample size of the trial
#' @param target       Target toxicity probability of DLT
#' @param tstar        End of the toxicity observation window t*
#' @param dose_init    Initial dose d1
#' @param rate         Expected number of arrivals per observation window
#' @param accrual      Inter-patient arrival scheme: "fixed" or "poisson"
#' @param skeleton1    Initial guesses of toxicity probabilities 
#' @param skeleton2    Initial guesses of progression probabilities 
#' @return             Simulation results of the survival-CRM
#' 

isurvcrm_sim <- function(N,scenario1,scenario2,n,target,tstar,dose_init,rate,accrual,skeleton1,skeleton2){
  K <- length(scenario1)
  diff <- as.numeric(as.character(abs(scenario1-target))) 
  correct <- which(diff== min(diff))
  toxic <- which(scenario1 > scenario1[correct])
  
  #Set of scaled doses 
  xref1 <-  log(-log(1-skeleton1-skeleton2)/(1+(skeleton2/skeleton1))*tstar) 
  xref2 <- -log(skeleton2*exp(xref1)/skeleton1)
  
  b1 <- log(log(-log(1-target)/tstar)/ xref1[dose_init])     #Beta value such that the probability of toxicity at dose level 1 is equal to the target
  
  CI_estim1 <- matrix(nrow=N,ncol=length(scenario1))  #Cumulative incidence of toxicity estimate at t* (mean over N simulations) 
  CI_estim2 <- matrix(nrow=N,ncol=length(scenario1))  #Cumulative incidence of discontinuation estimate at t* (mean over N simulations) 
  
  Dose_select <- rep(0,length(scenario1))           #No time each dose is selected over N simulations
  OD <- rep(NA,N)                                   #No patients treated at a toxic dose during a trial
  No_MTD <- rep(NA,N)                               #No patients treated at the MTD during a trial 
  No_DLT <- rep(NA,N)                               #No patients who experienced a DLT during a trial
  No_Stop <- rep(NA,N)                              #No trials stopped for safety decision
  
  
  s <- 1
  while (s<=N) {
    Base <- getdataset_CR(n, scenario1, scenario2, tstar,  seed=s)#Complete data generation for n patients
    #Patient accrual scheme
    if (accrual == "fixed") {
      next.arrival <- tstar/rate
    }   else if (accrual == "poisson") {
      next.arrival <- rexp(1, rate/tstar)
    }
    arrival <- u <- NULL #Arrival first patient
    
    i <- 1
    #Dose_init is given to the first patient
    Dataset_temp <- Base[Base[,'Patient'] ==1 & Base[,'Dose']==dose_init,c("Patient","Time", "Event", "Cens_time", "Dose"), drop = FALSE]
    
    
    while(i <=(n-1)){
      arrival <- c(arrival, next.arrival)
      if (i > 1){#Base_tmp (sequential data), built from the complete data (Base)
        Dataset_temp <- rbind(Dataset_temp,Base[Base[,'Patient']==i & Base[,'Dose']==next_dose , c("Patient","Time", "Event", "Cens_time", "Dose")])
      }
      if(Dataset_temp[i,"Time"]<=tstar){#if patient i experienced an event 
        unew <- Dataset_temp[i,"Time"]#we recover the time patient i made the event
      } else  {#patient i has no event
        unew <- Inf}
      u <- c(u,unew)
      uvent <- u+ arrival 
      
      if (accrual == "fixed") {
        next.arrival <- next.arrival + tstar/rate
      }    else if (accrual == "poisson") {
        next.arrival <- next.arrival + rexp(1, rate/tstar)
      } 
      Event <- rep(0, length(Dataset_temp[,"Event"])) #Vector of event occurrence
      Event[uvent <= next.arrival] <- Dataset_temp[uvent <= next.arrival,"Event"] 
      
      Cens <- pmin(next.arrival, uvent) - arrival #Censor indicator
      Time <- pmin(Cens, tstar) #Censored time-to-toxicity
      
      # dose-finding algorithm   
      res <- isurvcrm_nextdose(event=Event, dose_level=Dataset_temp[,'Dose'], xref1=xref1,xref2=xref2,time=Time,tstar,target)
      #Avoid skipping doses in escalation 
      next_dose <- min(res$next_dose, (Dataset_temp[i,'Dose'] + 1))
      
      # safety stopping rule
      if (integrate(denom_tox, -Inf, b1 ,event=Event, dose_level=Dataset_temp[,'Dose'], xref=xref1,time=Time,abs.tol = 0)$value / res$constante1 >=0.95) {
        break
      }
      
      i <- i+1
      
    } 
    #Completion of the trial with the last patient 
    if(i==n){
      
      Dataset_temp <- rbind(Dataset_temp,Base[Base[,'Patient']==i & Base[,'Dose']==next_dose , c("Patient","Time", "Event", "Cens_time", "Dose")])
      Event <- Dataset_temp[,"Event"]
      Time <- Dataset_temp[,"Cens_time"]
      
      #Model
      res <- isurvcrm_nextdose(event=Event, dose_level=Dataset_temp[,'Dose'], xref1=xref1,xref2=xref2,time=Time,tstar,target)
      
      CI_estim1[s,] <- res$CI_temp1 ##Cumulative incidence estimate at t* (for s)
      CI_estim2[s,] <- res$CI_temp2
      Dose_select[res$next_dose] <-  Dose_select[res$next_dose]+1 #dose selected for s
      No_MTD[s] <- sum(Dataset_temp[,"Dose"]%in%correct)
      OD[s] <- sum(Dataset_temp[,"Dose"]%in%toxic)
      No_DLT[s] <- sum(Dataset_temp[,"Event"])
    }
    else(No_Stop[s] <- i)
    
    s <- s+1
    print(s)
  }
  #nb of trials done = nb of trials that did not stop 
  No_trials <- sum(is.na(No_Stop))
  PS <- (Dose_select/No_trials)*100
  
  return(list(Scenario1=scenario1,
              Scenario2=scenario2,
              CI_estim1=colMeans(CI_estim1,na.rm = TRUE),
              CI_estim2=colMeans(CI_estim2,na.rm = TRUE),
              PS=PS,
              correct = correct,
              PCS = round(sum(PS[correct]),0),
              Accuracy_Index = round(1-length(scenario1)*(sum(abs(scenario1-target)*(Dose_select/No_trials))/sum(abs(scenario1-target))),2), 
              POS = round(sum(PS[toxic]), 0),
              R_Bais_MTD = round(mean((colMeans(CI_estim1,na.rm = TRUE)[correct]-scenario1[correct])/scenario1[correct] ) ,3), 
              OD=round(mean(OD,na.rm = TRUE),2),
              Nb_DLT = round(mean(No_DLT,na.rm = TRUE),2),
              Nb_MTD=round(mean(No_MTD,na.rm = TRUE),2),
              P_stop = (N-No_trials)/N*100,
              No_Stop = mean(No_Stop,na.rm = TRUE)
              
  ))
}     


# ----------------------------------------------------------------------------------------------
#' bmk_sim :       
#' Generate simulation replicates of phase I trial using the proposal benchmark
#' 
#' @param N            Number of simulation  
#' @param scenario     Vector of the true toxicity probabilites associated with the doses
#' @param n            Sample size of the trial
#' @param target       Target toxicity probability of DLT
#' @param tstar        End of the toxicity observation window t*
#' @param hazard       The time-to-toxicity hazard. "decreasing" for early toxicities; "increasing" for late toxicities or "constant"
#' @return             Simulation results of the survival-CRM
bmk_sim <- function(N,scenario,n,target,tstar,hazard){
  K <- length(scenario)
  diff <- as.numeric(as.character(abs(scenario-target))) 
  correct <- which(diff== min(diff))
  toxic <- which(scenario > scenario[correct])
  
  CI_estim <- matrix(nrow=N,ncol=length(scenario))  #Cumulative incidence estimate at t* (mean over N simulations) 
  Dose_select <- rep(0,length(scenario))            #No time each dose is selected over N simulations
  j = 1
  while(j<=N){
    set.seed(j)
    Ui=runif(n,0,1)
    Time <- matrix(nrow=n,ncol=length(scenario)) 
    Event <- matrix(nrow=n,ncol=length(scenario)) 
    Cens_time <- matrix(nrow=n,ncol=length(scenario))  
    i = 1
    while(i<=n){
      if(hazard=="constant"){
        Time[i,] <- qexp(Ui[i], rate = find_param(scenario,tstar,hazard) )
      } else{
        param <- find_param(scenario,tstar,hazard)
        Time[i,] <-  qweibull(Ui[i],param[,1],param[,2])
      }
      Event[i,] <-  ifelse(Time[i,] >= tstar,0,1) ##Indicateur de censure (0), évenement (1)
      i=i+1
    }
    d=1
    while(d<=length(scenario)){
      Cens_time[,d] <- pmin(Time[,d],tstar)
      #Si on a que des évenements qui arrivent après la censure, soit aucune tox pdt la fenetre d'observation. 
      #On fixe l'incidence cumulée estimée à 0. 
      if(dim(table(Event[,d]))==2) {
        km <- summary(survfit(Surv(Cens_time[,d],Event[,d])~ 1 ))
        CI_estim[j,d] <-  c(IC= 1-km$surv[length(km$surv)])
        }
      
    else{CI_estim[j,d] <- 0}
      d=d+1
    }
    #Next dose finding algorithm
    Diff_pestim_target <- abs(CI_estim[j,] - target)
    Doses_min <- which(Diff_pestim_target==min(Diff_pestim_target)  )
    MTD <- ifelse(length(Doses_min)==1,Doses_min,sample(Doses_min,1)) 
    Dose_select[MTD] <- Dose_select[MTD]+1
    
    
    j=j+1
    
  }


PS <- (Dose_select/N)*100               #Percent of selection of each dose over No_trials trials
  
  
  
  
  return(list(Scenario=round(scenario,2),
              PS=round(PS,0),
              PCS = round(sum(PS[correct]),0),
              Accuracy_index = round(1-length(scenario)*(sum(abs(scenario-target)*(Dose_select/N))/sum(abs(scenario-target))),2), 
              POS = round(sum(PS[toxic]), 0),
              R_Bais_MTD = round(mean((colMeans(CI_estim,na.rm = TRUE)[correct]-scenario[correct])/scenario[correct] ) ,3)
              
  ))
}











 
# ----------------------------------------------------------------------------------------------
#' titecrm_sim :       
#' Generate simulation replicates of phase I trial using the TITE-CRM
#' 
#' @param N            Number of simulation  
#' @param scenario     Vector of the true toxicity probabilites associated with the doses
#' @param n            Sample size of the trial
#' @param target       Target toxicity probability of DLT
#' @param tstar        End of the toxicity observation window t*
#' @param dose_init    Initial dose d1
#' @param rate         Expected number of arrivals per observation window
#' @param accrual      Inter-patient arrival scheme: "fixed" or "poisson"
#' @param hazard       The time-to-toxicity hazard. "decreasing" for early toxicities; "increasing" for late toxicities or "constant"
#' @return             Simulation results of the survival-CRM
titecrm_sim <- function(N,scenario,n,target,tstar,dose_init,rate,accrual,hazard){
  K <- length(scenario)
  diff <- as.numeric(as.character(abs(scenario-target)))
  correct <- which(diff== min(diff))
  toxic <- which(scenario > scenario[correct])
  

  #Initial guesses of toxicity probabilities associated the doses
  skeleton <-  getprior(halfwidth=0.07, target, nu=3, nlevel=K, model = "empiric", intcpt = 3)
  xref <-  log(-log(1-skeleton)/tstar)              #Set of numerical labels for the doses investigated in the trial
  b1 <- log(log(target)/log(skeleton[dose_init]))   #Beta value such that the probability of toxicity at dose level 1 is equal to the target

    
  CI_estim <- matrix(nrow=N,ncol=length(scenario))  #Cumulative incidence estimate at t* (mean over N simulations) 
  Dose_select <- rep(0,length(scenario))            #No time each dose is selected over N simulations
  OD <- rep(NA,N)                                   #No patients treated at a toxic dose during a trial
  No_MTD <- rep(NA,N)                               #No patients treated at the MTD during a trial 
  No_DLT <- rep(NA,N)                               #No patients who experienced a DLT during a trial
  No_Stop <- rep(NA,N)                              #No trials stopped for safety decision
 
  s <- 1
  
  while (s<=N) {
    Dataset <- getdataset(n, scenario,  tstar, seed=s ,hazard)  #Complete data generation for n patients
    if (accrual == "fixed") {                       #Patient accrual scheme
      next_arrival <- tstar/rate
    }   else if (accrual == "poisson") {
      next_arrival <- rexp(1, rate/tstar)
    }
    arrival <- u <- NULL                            #Arrival of the first patient
    i <- 1
    Dataset_temp <- Dataset[Dataset[,'Patient']==1 & Dataset[,'Dose']==dose_init,c("Patient","Time", "Event", "Cens_time", "Dose"), drop = FALSE]
    
    while(i <=(n-1)){
      arrival <- c(arrival, next_arrival)
      if (i > 1){ 
        Dataset_temp <- rbind(Dataset_temp,Dataset[Dataset[,'Patient']==i & Dataset[,'Dose']==next_dose , c("Patient","Time", "Event", "Cens_time", "Dose")]) #Sequential data built from the complete Dataset
      }
      if(Dataset_temp[i,"Time"]<=tstar){            #If patient i experienced a toxicity
        unew <- Dataset_temp[i,"Time"]    
      } else  {                                     #Patient i has no toxicity
        unew <- Inf}
      u <- c(u,unew)
      uvent <- u+ arrival 
      
      if (accrual == "fixed") {          
        next_arrival <- next_arrival + tstar/rate
      }    else if (accrual == "poisson") {
        next_arrival <- next_arrival + rexp(1, rate/tstar)
      } 

      
      Tox <- rep(0, length(Dataset_temp[,"Event"])) #Vector of toxicity occurrence
      Tox[uvent <= next_arrival] <- 1 
      
      Cens <- pmin(next_arrival, uvent) - arrival    #Censor indicator
      time <- pmin(Cens, tstar)                   
      
      
      weights <- pmin(time/tstar, 1)
      weights[Tox == 1] <- 1
      
      xref <- skeleton[Dataset_temp[,'Dose']]
      
      #dose-findind algorithm
      res <- titecrm_nextdose(xref, Tox, weights,skeleton)
      next_dose <- res$next_dose   
      
      #Avoid skipping doses in escalation 
      next_dose <- min(next_dose, (Dataset_temp[i,'Dose'] + 1))  
      
      
      # safety stopping rule
      if (integrate(crmh, -Inf, b1 ,xref, Tox, weights, sqrt(1.34),abs.tol = 0)$value / res$constante >= 0.95) {
        break
      }
 
      i <- i+1
      
    } 
    #Completion of the trial with the last patient 
    if(i==n){
      Dataset_temp <- rbind(Dataset_temp,Dataset[Dataset[,'Patient']==i & Dataset[,'Dose']==next_dose , c("Patient","Time", "Event", "Cens_time", "Dose")])
      Tox <- Dataset_temp[,"Event"]
      weights <- rep(1,n)
      xref <- skeleton[Dataset_temp[,'Dose']]
      
      #Dose-finding algorithm
      res <- titecrm_nextdose(xref, Tox, weights,skeleton)
      next_dose <- res$next_dose

      CI_estim[s,] <- res$CI_temp                    #Cumulative incidence estimate at t* for trial s 
      Dose_select[next_dose] <-  Dose_select[next_dose]+1  
      No_MTD[s] <- sum(Dataset_temp[,"Dose"] %in% correct)
      OD[s] <- sum(Dataset_temp[,"Dose"] %in% toxic)
      No_DLT[s] <- sum(Dataset_temp[,"Event"])
    }
    else(No_Stop[s] <- i)
    
    s <- s+1
  }
       
  No_trials <- sum(is.na(No_Stop))                    #Number of terminated trials 
  PS <- (Dose_select/No_trials)*100                   #Percent of selection of each dose over No_trials trials
  
  return(list(Scenario=round(scenario,2),
              PS=round(PS,0),
              PCS = round(sum(PS[correct]),0),
              Accuracy_index = round(1-length(scenario)*(sum(abs(scenario-target)*(Dose_select/No_trials))/sum(abs(scenario-target))),2), 
              POS = round(sum(PS[toxic]), 0),
              R_Bais_MTD = round(mean((colMeans(CI_estim,na.rm = TRUE)[correct]-scenario[correct])/scenario[correct] ) ,3), 
              OD=round(mean(OD,na.rm = TRUE),2),
              No_DLT = round(mean(No_DLT,na.rm = TRUE),2),
              No_MTD=round(mean(No_MTD,na.rm = TRUE),2),
              P_stop = (N-No_trials)/N*100
  ))
  
}     



# -------------------------------------------------------------------------------------------
#' generate_scen 
#' Generate toxicity probability for simulation scenarios according the expected OR between two successive probabilities
#' 
#' @param p_init        The targeted probability of toxicity for the MTD
#' @param OR            The expected OR between two successive probabilities
#' @param no            The number of probability values to be generated  
#' @param place         The values of expected probablity relative to p_init "inf" or "sup"
#' @return              The complete simulation scenarios
#' 

generate_scen <- function(p_init, OR,no,place){
  if(place=="sup")  {
    res <- p_init
    i <- 1
    while (i<=(no-1)) {
      res <- c(res,(OR*res[i]/(1-res[i])/(1+OR*res[i]/(1-res[i]))
      ) )
      i <- i+1
    }
  }
  if(place=="inf")  {
    res <- c()
    res[(no+1)] <- p_init
    for(i in (no+1):1){
      res[i-1] <-  (res[i]/((1-res[i])*OR))/(1+res[i]/((1-res[i])*OR))
    }
    res <- res[1:no]
  }
  return(res)
}



# -------------------------------------------------------------------------------------------
#' print_res 
#' Print simulation results 
#' #' 
#' @param sim_obj simulation object 
#
print_res <- function(sim_obj){
  if  (unique(table(names(unlist(sim_obj))))>1){
    for (i in 1:length(sim_obj)){
      cat("Scenario", i, sep=" ",  "\n");
      cat("Simulation scenarios:          ", sim_obj[i][[1]]$Scenario, sep="  ",  "\n");
      cat("Selection percentage:          ", sim_obj[i][[1]]$PS, sep="  ",  "\n");
      cat("Correct selection percentage:  ", sim_obj[i][[1]]$PCS,  sep="  ",  "\n");
      cat("Accuracy Index:                ", sim_obj[i][[1]]$Accuracy_index,  sep="  ",  "\n");
      cat("Overdose selection percentage: ", sim_obj[i][[1]]$POS, sep="  ",  "\n");
      cat("Relative bias at the true MTD: ", sim_obj[i][[1]]$R_Bais_MTD, sep="  ",  "\n");
      cat("Average overdose number:       ", sim_obj[i][[1]]$OD, sep="  ",  "\n");
      cat("Number of observed DLT:        ", sim_obj[i][[1]]$No_DLT, sep="  ",  "\n");
      cat("Number of treated with MTD:    ", sim_obj[i][[1]]$No_MTD, sep="  ",  "\n");
      cat("Stopped trials percentage:     ", sim_obj[i][[1]]$P_stop, sep="  ",  "\n");
      cat(" ",  sep="  ",  "\n")
    }
  } else{
    
    cat("Simulation scenarios:          ", sim_obj$Scenario, sep="  ",  "\n");
    cat("Selection percentage:          ", sim_obj$PS, sep="  ",  "\n");
    cat("Correct selection percentage:  ", sim_obj$PCS,  sep="  ",  "\n");
    cat("Accuracy Index:                ", sim_obj$Accuracy_index,  sep="  ",  "\n");
    cat("Overdose selection percentage: ", sim_obj$POS, sep="  ",  "\n");
    cat("Relative bias at the true MTD: ", sim_obj$R_Bais_MTD, sep="  ",  "\n");
    cat("Average overdose number:       ", sim_obj$OD, sep="  ",  "\n");
    cat("Number of observed DLT:        ", sim_obj$No_DLT, sep="  ",  "\n");
    cat("Number of treated with MTD:    ", sim_obj$No_MTD, sep="  ",  "\n");
    cat("Stopped trials percentage:     ", sim_obj$P_stop, sep="  ",  "\n");
    
  }
}

