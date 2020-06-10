# A two-stage drop-the-losers design for time-to-event outcome using a historical control arm
# by R. Abbas, J. Wason, S. Michiels, and G. Le Teuff (2020)
# Phase II clinical trials are simulated according to a drop-the-losers and a fixed designs both using historical control arm.

# the following code allow to reproduce results of the simulation study for the fixed design with or without correction for multiplicity
# as shown in tables 6 -7 of the paper.

# source the functions
# setwd()
source(file ="MANI_function.r" ,echo = FALSE)

# simulations parameters

# seed: seed of the random number generator 
# wei_shape: shape parameter of the weibull model fitted on historical data
# hazard.censoring: experimental arms hazards (according to a specified scenario)
# n.total: number of patients in total
# alpha: the nominal family-wise error rate
# niterations: number of iterations
# K: the number of experimental arms
# recruitmentrate: number of patients enrolled to the trial per unit of time
# followuptime: required follow-up time from the inclusion of the last patient to the final analysis      
# hazard.censoring: censoring hazard (exponential distribution)


# simulation plan
# specify the treatment effects
hazard.null  = -log(0.4)  / 2
hazard.small = -log(0.50) / 2  # small  effect
hazard.medium= -log(0.58) / 2  # medium effect
hazard.large = -log(0.65) / 2  # large  effect
# scenarios
s1	=	c( hazard.null,	  hazard.null,	 hazard.null) #1: No better than HC
s2	=	c( hazard.small,	hazard.null,	 hazard.null) #2: 1 better than HC small  effect
s3	=	c( hazard.medium,	hazard.null,	 hazard.null) #3: 1 better than HC medium effect
s4	=	c( hazard.large,	hazard.null,	 hazard.null) #4: 1 better than HC large effect
s5	=	c( hazard.small,	hazard.small,	 hazard.null) #5: 2 better than HC small effects
s6	=	c( hazard.medium,	hazard.small,	 hazard.null) #6: 2 better than HC (small/medium) effects
s7	=	c( hazard.medium,	hazard.medium, hazard.null) #7: 2 better than HC medium effects
s8	=	c( hazard.large,	hazard.medium, hazard.null) #8: 2 better than HC (medium/large) effects
s9	=	c( hazard.large,	hazard.large,  hazard.null) #9: 2 better than HC large effects
scenarios <- list(s1,s2,s3,s4,s5,s6,s7,s8,s9)

# Simulate a fixed design with 40 patients under the scenario 7
sim <- mani(seed=2306,wei_shape=0.9509, hazard.exp=scenarios[[7]], n.total=40, alpha=0.1, niterations=10000, K=3, recruitmentrate = 50, followuptime = 2, hazard.censoring = 0.5, hazard.null = hazard.null)
sim_hommel <- mani_hommel(seed=2306,wei_shape=0.9509, hazard.exp=scenarios[[7]], n.total=40, alpha=0.1, niterations=10000, K=3, recruitmentrate = 50, followuptime = 2, hazard.censoring = 0.5, hazard.null = hazard.null)
sim$disjunctive.power
sim_hommel$disjunctive.power
