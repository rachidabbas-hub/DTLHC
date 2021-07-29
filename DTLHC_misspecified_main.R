# A two-stage drop-the-losers design for time-to-event outcome using a historical control arm
# by R. Abbas, J. Wason, S. Michiels, and G. Le Teuff (2021)
# Phase II clinical trials are simulated according to a drop-the-losers and a fixed designs both using historical control arm.

# the following code allow to reproduce results of the simulation study for two-stage drop-the-losers design using a historical control arm
# in the misspecified case as shown in tables 4- 5 of the paper.

# source the functions
# setwd()
source(file ="DTLHC_misspecified_functions.r" ,echo = FALSE) # here data are generated according to a gamma

#Fit historical control data to a gamma distribution
# Call:
#   flexsurvreg(formula = Surv(time, death) ~ 1, data = hdata, dist = "gengamma")
# 
# Estimates: 
#   est       L95%      U95%      se      
# mu      0.20509  -0.06522   0.47540   0.13791
# sigma   1.23705   1.11343   1.37440   0.06645
# Q      -0.40050  -0.80898   0.00799   0.20841
# 
# N = 265,  Events: 189,  Censored: 76
# Total time at risk: 481.3826
# Log-likelihood = -345.3897, df = 3
# AIC = 696.7793
q=-0.4005
s=1.23705
# we will simulate historical data according to a gamma distribution with above parameters estimated from the historical cohort

# simulations parameters

# seed: seed of the random number generator 
# wei_shape: shape parameter of the weibull model fitted on historical data
# hazard.censoring: experimental arms hazards (according to a specified scenario)
# n1.perarm: number of patients included per arm in stage 1
# n2: number of patients included in stage 2
# requiredfwer: the nominal family-wise error rate
# niterations: number of iterations
# K: the number of experimental arms
# recruitmentrate: number of patients enrolled to the trial per unit of time
# followuptime: required follow-up time from the inclusion of the last patient to the final analysis      
# hazard.censoring: censoring hazard (exponential distribution)
#hazard.null: event hazard parameter under the null hypothesis


# simulation plan
# specify the treatment effects
mu.null  = uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.40)$root
mu.small = uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.50)$root # small effect
mu.medium= uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.58)$root # medium effect
mu.large = uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.65)$root # large effect
# specify censoring
mu.censoring=uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=exp(-0.5))$root 
# specify parameters from historical data: exponential function (misspecification)
hazard.null=-log(0.4)/2
# scenarios
s1 = c( mu.null,	mu.null,   mu.null) #1: No better than HC
s2 = c( mu.small,	mu.null,   mu.null) #2: 1 better than HC small  effect
s3 = c( mu.medium,	mu.null,   mu.null) #3: 1 better than HC medium effect
s4 = c( mu.large,	mu.null,   mu.null) #4: 1 better than HC large effect
s5 = c( mu.small,	mu.small,  mu.null) #5: 2 better than HC small effects
s6 = c( mu.medium,	mu.small,  mu.null) #6: 2 better than HC (small/medium) effects
s7 = c( mu.medium,	mu.medium, mu.null) #7: 2 better than HC medium effects
s8 = c( mu.large,	mu.medium, mu.null) #8: 2 better than HC (medium/large) effects
s9 = c( mu.large,	mu.large,  mu.null) #9: 2 better than HC large effects
scenarios <- list(s1,s2,s3,s4,s5,s6,s7,s8,s9)


# Reproduce results under scenario 3 for the drop-thelosers design with 30 patients per arm at first stage and 30 patients at second stage
sim <- droptheloser(seed=2404 ,mu.exp=s3,n1.perarm=30,n2=30, requiredfwer=0.1, K=3,mu.null=mu.null,hazard.null=hazard.null, mu.censoring = mu.censoring, followuptime=2, recruitmentrate=50)
sim$disjunctive.power
sim$theoric.d.power
