## simulation: drop the loser design with historical control
## here we assess the impact of misspecified functionnal form of historical data
source(file ="DTLHC-gamma_functions.r" ,echo = FALSE) # here data are generated according to a gamma

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
# we will simulate historical data according to a gamma distribution with above parameters
# it is the parameter estimated from the historical cohort

# simulations parameters
K=3   # number of exprimental arms

mu.censoring=uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=exp(-0.5))$root # hazard of censoring 
followuptime=2
interimanalysisfollowup=0    # time between inclusion of the last patient in stage j and the jth interim analysis
recruitmentrate=50
nit=10000 #number of iterations

# specify the treatment effects
mu.null  = uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.40)$root
mu.small = uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.50)$root # small effect
mu.medium= uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.58)$root # medium effect
mu.large = uniroot(f = find.parameter,lower = -2,upper=10,q=-0.4005,s=1.23705,t=2,S=0.65)$root # large effect

# (mis-)specify parameters from historical data: exponential function
hazard.null=-log(0.4)/2

## simulation plan
# Nominal fwer
requiredfwer=0.1

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

# sample size per arm
sample.sizes <- cbind(
  c( 10,15,20,25,30,40,50),
  c( 10,15,20,25,30,40,50)
)

# run simulations
set.seed(2404)
sim <- tab <-  list()
tab.sel1 <- tab.sel2 <- matrix(0,nrow(sample.sizes),7) 

for (ss in (1:nrow(sample.sizes))){
for (s in 1:length(scenarios)) {
sim[[s]] = droptheloser(mu.exp=scenarios[[s]], n1.perarm=sample.sizes[ss,1],n2=sample.sizes[ss,2], requiredfwer=requiredfwer, K=K,mu.null=mu.null,hazard.null=hazard.null, mu.censoring = mu.censoring, followuptime=followuptime, interimanalysisfollowup=interimanalysisfollowup,  recruitmentrate=recruitmentrate)
}
tab.line <- data.frame(cbind(
  sim[[1]]$efwer,
  sim[[2]]$disjunctive.power,
  sim[[3]]$disjunctive.power,
  sim[[4]]$disjunctive.power,
  sim[[5]]$disjunctive.power,
  sim[[6]]$disjunctive.power,
  sim[[7]]$disjunctive.power,
  sim[[8]]$disjunctive.power,
  sim[[9]]$disjunctive.power))
colnames(tab.line) <- paste("scenario ",1:ncol(tab.line))
tab[[ss]] <- tab.line 

tab.sel1[ss,] <- cbind(
  Nt = K*sample.sizes[ss,1]+sample.sizes[ss,2],
  n1.perarm=sample.sizes[ss,1],
  n2=sample.sizes[ss,2],
  scenario=6,
  p.sel.1=prop.table(table(sim[[6]]$selectedarm))[1],
  p.sel.2=prop.table(table(sim[[6]]$selectedarm))[2],
  p.sel.r=prop.table(table(sim[[6]]$selectedarm))[1]/prop.table(table(sim[[6]]$selectedarm))[2]
)

tab.sel2[ss,] <- cbind(
  Nt = K*sample.sizes[ss,1]+sample.sizes[ss,2],
  n1.perarm=sample.sizes[ss,1],
  n2=sample.sizes[ss,2],
  scenario=8,
  p.sel.1=prop.table(table(sim[[8]]$selectedarm))[1],
  p.sel.2=prop.table(table(sim[[8]]$selectedarm))[2],
  p.sel.r=prop.table(table(sim[[8]]$selectedarm))[1]/prop.table(table(sim[[8]]$selectedarm))[2]
)
}
table.power <- cbind(Nt = K*sample.sizes[,1]+sample.sizes[,2], N_1 = sample.sizes[,1], N_2 = sample.sizes[,2], Ratio = round(sample.sizes[,1]/sample.sizes[,2],digits = 1),do.call(rbind,tab) )
head(table.power)

head(table.power[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),] )
head(tab.sel1[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),] )
head(tab.sel2[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),] )

table.power[table.power$Ratio==1 & table.power$Nt!=20 ,-4]
table.power[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),]

tab.sel1[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),] 
tab.sel2[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),] 

summary(tab.sel1[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),7] )
summary(tab.sel2[c(7, 9, 10, 14, 17, 19, 24, 28, 31, 35, 38, 41, 46, 50, 53, 58, 62, 66, 71, 75, 79, 84, 88, 91),7] )

# save
filename <- paste("DHCmisspec_simulation-report_",Sys.Date(), sep="")
save(table.power, tab.sel1,tab.sel2,
     K,mu.censoring,followuptime,interimanalysisfollowup,recruitmentrate,nit,mu.null,
     mu.small,mu.medium,mu.large,q,s,requiredfwer,scenarios,sample.sizes,
     file = paste(filename,"RDATA",sep=".")) 
