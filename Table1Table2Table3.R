# -------------------------------------------------------------------------------------------
##### This code can be used to reproduce simulation results for scenario 1 to 12       ##### 
##### with the survival-CRM, isurvival-CRM, TITE-CRM and benchmark given in            ##### 
#####                                                                                  ##### 
##### Table 1, 2, and 3                                                                ##### 
#####                                                                                  ##### 
##### of the paper ``Dose-finding design and benchmark for a right censored endpoint'' ##### 
#####                                                                                  ##### 
##### by Anais Andrillon, Sylvie Chevret, Shing M Lee and Lucie Biard (2020)           ##### 
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
##### Input
# 
# Number of Simulations
N <- 10000
# Sample Size
n <-25
#Target toxicity probability of DLT
target <- 0.25
#End of the toxicity observation window t*
tstar <- 1
#Initial dose d1
dose_init <- 1
#Expected number of arrivals per observation window
rate <- 4
#Inter-patient arrival scheme: "fixed" or "poisson"
accrual <- "fixed"
#The prior guess of MTD.
nu <- 3 
#The desired halfwidth of the indifference intervals
halfwidth <- 0.05
# -------------------------------------------------------------------------------------------
##### Simulation scenarios for toxicity
# 
#Cumulative incidence of DLT at time t* by dose.
#Sc1 to Sc5: scenarios with OR=2 ; MTD shifted from dose 1 to dose 5
OR <- 2
no <- 5
S1 <- generate_scen(target, OR,no,place="sup")
S2 <-c(generate_scen(target, OR,no=1,place="inf"),generate_scen(target, OR,no=4,place="sup"))
S3 <-c(generate_scen(target, OR,no=2,place="inf"),generate_scen(target, OR,no=3,place="sup"))
S4 <-c(generate_scen(target, OR,no=3,place="inf"),generate_scen(target, OR,no=2,place="sup"))
S5 <-c(generate_scen(target, OR,no=4,place="inf"),generate_scen(target, OR,no=1,place="sup"))
#Sc6: steeper dose-toxicity curve around the MTD
S6 = c(0.01,0.25,0.55,0.78,0.95) 
#Sc7: similar scenario to Sc3, with a less steep
S7 <- c( 0.13, 0.18, 0.23, 0.28, 0.33)
#Sc7: similar scenario to Sc1, with a less steep
S8 <- c(seq(0.30,0.5,length.out = 5))
#Sc9: plateau
S9 <- c(0.01,0.13,rep(0.25,3))
scenarios <- list(S1,S2,S3,S4,S5,S6,S7,S8,S9)



#-----------------------------------------------------------------------------------------
##### Simulation scenarios for progression
#Sc10: high risks of discontinuation; cumulative incidence at t* decreasing from 0.5 to 0.4
S10 <- seq(0.5,0.4,length.out = 5)
#Sc11: moderate risks of discontinuation; cumulative incidence at t* decreasing from 0.3 to 0.1
S11 <- seq(0.3,0.1,length.out = 5)
#Sc12: cumulative incidence ranges from 0.55 to 0.3
S12 <- seq(0.55,0.3,length.out = 5)
scenarios.prog <- list(S10, S11, S12)



# -------------------------------------------------------------------------------------------
##### Running the survival-CRM using parallel processing over several scenarios
#
#Example with simulation scenario 1
#Running 1,000 trials with survival-CRM design for Sc1 with constant hazards of toxicity over time 
#print_res(survcrm_sim(N=1000,S1,n,target,tstar,dose_init,rate,accrual,hazard="constant",nu,halfwidth))
  

#Running N=10,000 replicates of the survival-CRM for all scenarios with constant hazards of toxicity over time
library(parallel)

survcrm_result_constant <- mclapply(scenarios, survcrm_sim,
                                    N=N,n=n,target=target,tstar=tstar,dose_init=dose_init, rate=rate,accrual=accrual,hazard="constant",nu=nu,halfwidth=halfwidth,
                                    mc.cores = 9)

print.res(survcrm_result_constant)

# could take several hours
# mclapply function enable parallel processing. By default, mclapply will use all cores available to it. 
# If you have enough cores available, I advise you to use as many cores as the number of scenarios you run
# Setting mc.cores to 1 disables parallel processing. 


# Running the survival-CRM for all scenarios with increasing hazards of toxicity over time, i.e. late onset toxicities
survcrm.result.increasing <- mclapply(scenarios, survcrm_sim,
                                      N=N,n=n,target=target,tstar=tstar,dose_init=dose_init, rate=rate,accrual=accrual,hazard="increasing",nu=nu,halfwidth=halfwidth,
                                      mc.cores = 9)
print.res(survcrm.result.increasing)

# Running the survival-CRM for all scenarios with decreasing hazards of toxicity over time, i.e. early onset toxicities
survcrm.result.decreasing <- mclapply(scenarios, survcrm_sim,
                                      N=N,n=n,target=target,tstar=tstar,dose_init=dose_init, rate=rate,accrual=accrual,hazard="decreasing",nu=nu,halfwidth=halfwidth,
                                      mc.cores = 9)

print.res(survcrm.result.decreasing)


# -------------------------------------------------------------------------------------------
##### Running the TITE-CRM
#Example with simulation scenario 1
#Running 1,000 trials with TITE-CRM design for Sc1 with constant hazards of toxicity over time 
# print_res(titecrm_sim(N=1000,S1,n,target,tstar,dose_init,rate,accrual,hazard="constant"))

#Running N=10,000 replicates of the TITE-CRM for all scenarios with constant hazards of toxicity over time
titecrm.result.constant <- mclapply(scenarios, titecrm_sim,
                                    N=N,n=n,target=target,tstar=tstar,dose_init=dose_init, rate=rate,accrual=accrual,hazard="constant",
                                    mc.cores = 9)
print.res(titecrm.result.constant)



# Running the TITE-CRM for all scenarios with increasing hazards of toxicity over time, i.e. late onset toxicities
titecrm.result.increasing <- mclapply(scenarios, titecrm_sim,
                                      N=N,n=n,target=target,tstar=tstar,dose_init=dose_init, rate=rate,accrual=accrual,hazard="increasing",
                                      mc.cores = 9)
print.res(titecrm.result.increasing)

# Running the TITE-CRM for all scenarios with decreasing hazards of toxicity over time, i.e. early onset toxicities
titercm.result.decreasing <- mclapply(scenarios, titecrm_sim,
                                      N=N,n=n,target=target,tstar=tstar,dose_init=dose_init, rate=rate,accrual=accrual,hazard="decreasing",
                                      mc.cores = 9)

print.res(titercm.result.decreasing)


# -------------------------------------------------------------------------------------------
##### Running the continuous benchmark
#Example with simulation scenario 1
#Running 1,000 trials with benchmark for Sc1 with constant hazards of toxicity over time 
print_res(bmk_sim(N=1000,S1,n,target,tstar,hazard="constant"))

#Running N=10,000 replicates of the benchmark for all scenarios 
bmk.result <- mclapply(scenarios, bmk_sim,
                                    N=N,n=n,target=target,tstar=tstar,hazard="constant",
                                    mc.cores = 1)
print.res(bmk.result)


# -------------------------------------------------------------------------------------------
##### Running the iSurvival-CRM 
#Initial guesses of toxicity probabilities 
skeleton1 <- c(0.05, 0.10, 0.25, 0.35, 0.50)
#Initial guesses of progression probabilities 
skeleton2 <- c(0.50, 0.45, 0.40, 0.35, 0.35)
 

#Example with simulation scenario 1
#Running 1,000 trials with isurvival-CRM design for Sc10 with constant hazards of toxicity over time 
print_res(isurvcrm_sim(N=1000,scenario1=S3,scenario2=S10,n,target,tstar,dose_init,rate,accrual,skeleton1=skeleton1,skeleton2=skeleton2))


# Running the survival-CRM for all scenarios with decreasing hazards of toxicity over time, i.e. early onset toxicities
isurvcrm.result <- mclapply(scenarios.prog, isurvcrm_sim,
                            scenario1=S3,N=N,n=n,target=target,tstar=tstar,dose_init=dose_init, rate=rate,accrual=accrual,skeleton1=skeleton1,skeleton2=skeleton2,
                                    mc.cores = 3)
 

print.res(isurvcrm.result)





