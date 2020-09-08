# -------------------------------------------------------------------------------------------
##### This code can be used to implement the survival-CRM,                             ##### 
##### the informative survival-CRM and benchmark for Phase 1                           ##### 
##### dose-finding trials with right censored time-to-event endpoints;                 #####  
##### and to reproduce the the result in the paper ``Dose-finding                      ##### 
##### design and benchmark for a right censored endpoint''                             ##### 
#####                                                                                  ##### 
##### by Andrillon, Chevret, Lee and Biard (2020)                                      ##### 
# -------------------------------------------------------------------------------------------


library(dplyr) #function between



# -------------------------------------------------------------------------------------------
#' getprior.exp 
#' Returns a vector of initial guesses of toxicity probabilities associated the doses 
#' for the exponential working model for a given set of indifference intervals. 
#' 
#' @param halfwidth      The desired halfwidth of the indifference intervals
#' @param target         The target toxicity probability of DLT
#' @param nu             The prior guess of MTD
#' @param nlevel         The number of test doses 
#' @param tstar          The end of toxicity observation window t*
#' @return               A vector of initial guesses of toxicity probabilities 
#' 
getprior.exp <- function (halfwidth, target, nu, nlevel, tstar) {
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
  Dataset <- cbind(Patient=rep(1:n, each=K),Ui= rep(runif(n), each=K),Time=rep(NA,n*K),Event=rep(NA,n*K),Cens.time=rep(NA,n*K),Dose=rep(1:K,n))  
  i = 1
  while(i<=n){
    if(hazard=="constant"){
      Dataset[(K*(i-1)+1):(K*i),'Time'] <-  qexp(Dataset[(K*(i-1)+1):(K*i),'Ui'], rate = find.param(scenario,tstar,hazard))
    } else{
      param <- find.param(scenario,tstar,hazard)
      Dataset[(K*(i-1)+1):(K*i),'Time'] <-  qweibull(Dataset[(K*(i-1)+1):(K*i),'Ui'],param[,1],param[,2])
    }
    Dataset[(K*(i-1)+1):(K*i),'Event'] <-  ifelse(Dataset[(K*(i-1)+1):(K*i),'Time'] >= tstar,0,1) #Event: censor indicator. 1 if DLT before t*
    Dataset[(K*(i-1)+1):(K*i),'Cens.time'] <- pmin( Dataset[(K*(i-1)+1):(K*i),'Time'] ,tstar )   #Cens.time=1 if no event,  Cens.time=Time if DLT before t*
    
    i=i+1
  }
  return(Dataset)
  
}



# -------------------------------------------------------------------------------------------
#' find.param
#' From True proba of DLT, find parameters of time-to-event distribution 
#' 
#' @param p             The true probability of DLT 
#' @param tstar         The end of toxicity observation window t*
#' @param hazard        The time-to-toxicity hazard. "decreasing" for early toxicities; "increasing" for late toxicities or "constant"
#' @return              Distribution parameters of the time-to-event distribution 
#' 
find.param <- function(p, tstar, hazard){
  if(hazard!="constant"){
   a <- ifelse(hazard=="increasing",3,0.5)
   b = exp(1/a * log(-tstar^a/log(1-p)))
   param <- cbind(rep(a, length(p)), b)
  } else {param <- (- log(1-p)/tstar)}
  return(param)
}
  
 

# -------------------------------------------------------------------------------------------
#' likelihood.tox.exp 
#' Survival likelihood
#' 
#' @param beta         The unknown model parameter
#' @param event        The binary vector of toxicity occurrence
#' @param dose.level   The dose levels allocated to patients
#' @param xref         The set of numerical labels for the doses 
#' @param time         The right-censored time-to-toxicity
#' @return             The likelihood value
#' 
likelihood.tox.exp <- function(beta, event, dose.level, xref,time){
  res <- 1
  for (i in 1:length(event)){
    res <- res * (( (exp(exp(beta)*xref[dose.level[i]])) * (exp(-exp(xref[dose.level[i]]*exp(beta))*time[i])) )^event[i] * (exp(-exp(xref[dose.level[i]]*exp(beta))*time[i]))^(1-event[i]) )
  }
  return(res)
}

denom.tox <- function(beta, event, dose.level, xref,time){ likelihood.tox.exp(beta, event, dose.level, xref,time)*dnorm(beta, mean=0, sd=sqrt(1.34)) }
num.tox <- function(beta, event, dose.level, xref,time){beta*denom.tox(beta, event, dose.level, xref,time)}


# ----------------------------------------------------------------------------------------------
#' survcrm.sim :       
#' Generate simulation replicates of phase I trial using the proposal survival-CRM
#' 
#' @param N            Number of simulation  
#' @param scenario     Vector of the true toxicity probabilites associated with the doses
#' @param n            Sample size of the trial
#' @param target       Target toxicity probability of DLT
#' @param tstar        End of the toxicity observation window t*
#' @param dose.init    Initial dose d1
#' @param rate         Expected number of arrivals per observation window
#' @param accrual      Inter-patient arrival scheme: "fixed" or "poisson"
#' @param hazard       The time-to-toxicity hazard. "decreasing" for early toxicities; "increasing" for late toxicities or "constant"
#' @return             Simulation results of the survival-CRM
#' 
survcrm.sim <- function(N,scenario,n,target,tstar,dose.init,rate,accrual,hazard){
  K <- length(scenario)
  diff <- as.numeric(as.character(abs(scenario-target))) 
  correct <- which(diff== min(diff))
  toxic <- which(scenario > unique(scenario[correct]))
  
  prior <- getprior.exp(halfwidth=0.05, target, nu=3, nlevel=K, tstar)
  
  xref <-  log(-log(1-prior)/tstar) 
  b1 <- log(log(-log(1-target))/xref[dose.init])#-0.7499267
  
  CI.estim <- matrix(nrow=N,ncol=length(scenario))  #Cumulative incidence estimate at t* (mean over N simulations) 
  Dose.select <- rep(0,length(scenario))            #Number of time each dose is selected over N simulations
  OD <- rep(NA,N)                                   #Number of patients treated at a toxic dose during a trial
  No.MTD <- rep(NA,N)                               #Number of patients treated at the MTD during a trial 
  No.DLT <- rep(NA,N)                               #Number of patients who experienced a DLT during a trial
  duration <-  rep(NA,N)                            #Trial duration
  No.Stop <- rep(NA,N)                              #Number of trials stopped for safety decision
  
  s <- 1
  while (s<=N) {
    Dataset <- getdataset(n, scenario,  tstar, seed=s ,hazard)  #Complete data generation for n patients
    if (accrual == "fixed") {                       #Patient accrual scheme
      next.arrival <- tstar/rate
    }   else if (accrual == "poisson") {
      next.arrival <- rexp(1, rate/tstar)
    }
    arrival <- u <- NULL                            #Arrival of the first patient
    i <- 1
    Dataset.temp <- Dataset[Dataset[,'Patient']==1 & Dataset[,'Dose']==dose.init,c("Patient","Time", "Event", "Cens.time", "Dose"), drop = FALSE]
    
    while(i <=(n-1)){
      arrival <- c(arrival, next.arrival)
      if (i > 1){ 
        Dataset.temp <- rbind(Dataset.temp,Dataset[Dataset[,'Patient']==i & Dataset[,'Dose']==next.dose , c("Patient","Time", "Event", "Cens.time", "Dose")]) #Sequential data built from the complete Dataset
      }
      if(Dataset.temp[i,"Time"]<=tstar){            #If patient i experienced a toxicity
        unew <- Dataset.temp[i,"Time"]    
      } else  {                                     #Patient i has no toxicity
        unew <- Inf}
      u <- c(u,unew)
      utox <- u+ arrival 
      
      if (accrual == "fixed") {          
        next.arrival <- next.arrival + tstar/rate
      }    else if (accrual == "poisson") {
        next.arrival <- next.arrival + rexp(1, rate/tstar)
      } 
      
      Tox <- rep(0, length(Dataset.temp[,"Event"])) #Vector of toxicity occurrence
      Tox[utox <= next.arrival] <- 1 
      
      Cens <- pmin(next.arrival, utox) - arrival    #Censor indicator
      time <- pmin(Cens, tstar)                   
      
      #Model parameter estimation
      constante <- integrate(denom.tox , -Inf, Inf, event=Tox, dose.level=Dataset.temp[,'Dose'], xref=xref,time=time,abs.tol = 0)$value
      beta_hat <-   integrate(num.tox, -Inf, Inf,  event=Tox, dose.level=Dataset.temp[,'Dose'], xref=xref,time=time,abs.tol = 0)$value / constante
      
      # safety stopping rule
      if (integrate(denom.tox, -Inf, b1 ,event=Tox, dose.level=Dataset.temp[,'Dose'], xref=xref,time=time,abs.tol = 0)$value / constante >=0.95) {
        break
      }
      
      lambda <- exp(exp(beta_hat)*xref)            #Updated model parameter
      CI.temp <- 1-exp(- lambda * tstar)           #Cumulative incidence at t* updated
      
      #Next dose finding algorithm
      Diff_pestim_target <- abs(CI.temp - target)
      Doses_min <- which(Diff_pestim_target==min(Diff_pestim_target)  )
      next.dose <- ifelse(length(Doses_min)==1,Doses_min,sample(Doses_min,1)) 
      
      
      #Calculation of the trial duration 
      if(i==1){ duration.temp <- time} else{
        duration.temp <- duration.temp+time[i]
      }
      
      #Avoid skipping doses in escalation 
      next.dose <- min(next.dose, (Dataset.temp[i,'Dose'] + 1))  
       
      
      i <- i+1
      
    } 
    
    #Completion of the trial with the last patient 
    if(i==n){
      Dataset.temp <- rbind(Dataset.temp,Dataset[Dataset[,'Patient']==i & Dataset[,'Dose']==next.dose , c("Patient","Time", "Event", "Cens.time", "Dose")])
      Tox <- Dataset.temp[,"Event"]
      time <- Dataset.temp[,"Cens.time"]
      
      #Model parameter estimation
      constante <- integrate(denom.tox , -Inf, Inf, event=Tox, dose.level=Dataset.temp[,'Dose'], xref=xref,time=time,abs.tol = 0)$value
      beta_hat <-   integrate(num.tox, -Inf, Inf,  event=Tox, dose.level=Dataset.temp[,'Dose'], xref=xref,time=time,abs.tol = 0)$value / constante
      
      lambda <- exp(exp(beta_hat)*xref)            #Updated model parameter
      CI.temp <- 1-exp(- lambda * tstar)           #Cumulative incidence at t* updated
      
      #Dose level Choice for trial s
      dose.estim <- which(abs(CI.temp - target)==min(abs(CI.temp - target)))
      next.dose <- ifelse(length(dose.estim)==1,dose.estim,sample(dose.estim,1)) 
      
      duration[s] <- duration.temp+time[i]
      
      CI.estim[s,] <- CI.temp                      #Cumulative incidence estimate at t* for trial s 
      Dose.select[next.dose] <-  Dose.select[next.dose]+1  
      No.MTD[s] <- sum(Dataset.temp[,"Dose"] %in% correct)
      OD[s] <- sum(Dataset.temp[,"Dose"] %in% toxic)
      No.DLT[s] <- sum(Dataset.temp[,"Event"])
    }
    else(No.Stop[s] <- i)
    
    s <- s+1
  }
  No.trials <- sum(is.na(No.Stop))                    #Number of terminated trials 

  PS <- (Dose.select/No.trials)*100                   #Percent of selection of each dose over No.trials trials
  
  
  return(list(Scenario=round(scenario,2),
              PS=round(PS,0),
              PCS = round(sum(PS[correct]),0),
              Accuracy.index = round(1-length(scenario)*(sum(abs(scenario-target)*(Dose.select/No.trials))/sum(abs(scenario-target))),2), 
              POS = round(sum(PS[toxic]), 0),
              R.Bais.MTD = round(mean((colMeans(CI.estim,na.rm = TRUE)[correct]-scenario[correct])/scenario[correct] ) ,3), 
              OD=round(mean(OD,na.rm = TRUE),2),
              No.DLT = round(mean(No.DLT,na.rm = TRUE),2),
              No.MTD=round(mean(No.MTD,na.rm = TRUE),2),
              P.stop = (N-No.trials)/N*100,
              Duration=round(mean(duration,na.rm = TRUE),2)
  ))
  
}     


# -------------------------------------------------------------------------------------------
#' print.res 
#' Print simulation results 
#' #' 
#' @param sim.obj simulation object 
#
print.res <- function(sim.obj){
  cat("Simulation scenarios:          ", sim.obj$Scenario, sep="  ",  "\n");
  cat("Selection percentage:          ", sim.obj$PS, sep="  ",  "\n");
  cat("Correct selection percentage:  ", sim.obj$PCS,  sep="  ",  "\n");
  cat("Accuracy Index:                ", sim.obj$Accuracy.index,  sep="  ",  "\n");
  cat("Overdose selection percentage: ",  sim.obj$POS, sep="  ",  "\n");
  cat("Relative bias at the true MTD: ",  sim.obj$R.Bais.MTD, sep="  ",  "\n");
  cat("Average overdose number:       ",  sim.obj$OD, sep="  ",  "\n");
  cat("Number of observed DLT:        ",  sim.obj$No.DLT, sep="  ",  "\n");
  cat("Number of treated with MTD:    ", sim.obj$No.MTD, sep="  ",  "\n");
  cat("Stopped trials percentage:     ",  sim.obj$P.stop, sep="  ",  "\n");
  cat("Duration of the trial:         ", sim.obj$Duration,sep="  ",  "\n")
}
 


# -------------------------------------------------------------------------------------------
#' generate.scen 
#' Generate toxicity probability for simulation scenarios with the MTD shifted from dose 1 to dose 5
#' 
#' @param p.init        The targeted probability of toxicity for the MTD
#' @param OR            The expected OR between  two successive probabilities
#' @param no            The number of probability values to be generated  
#' @param place         The values of expected probablity relative to p.init "inf" or "sup"
#' @return              The complete simulation scenarios
#' 

generate.scen <- function(p.init, OR,no,place){
  if(place=="sup")  {
    res <- p.init
    i <- 1
    while (i<=(no-1)) {
      res <- c(res,(OR*res[i]/(1-res[i])/(1+OR*res[i]/(1-res[i]))
      ) )
      i <- i+1
    }
  }
  if(place=="inf")  {
    res <- c()
    res[(no+1)] <- p.init
    for(i in (no+1):1){
      res[i-1] <-  (res[i]/((1-res[i])*OR))/(1+res[i]/((1-res[i])*OR))
    }
    res <- res[1:no]
  }
  return(res)
}


