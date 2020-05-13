## simulation: drop the loser design with historical control
source(file ="DTLHC_functions.r" ,echo = FALSE)

# simulations parameters
K=3                   # number of exprimental arms
hazard.censoring=0.5  # hazard of censoring (exponential distribution)
followuptime=2        
interimanalysisfollowup=0 # time between inclusion of the last patient in stage j and the jth interim analysis
recruitmentrate=50    
nit=10000 #number of iterations

# specify the treatment effects
hazard.null  = -log(0.4)  / 2
hazard.small = -log(0.50) / 2  # small  effect
hazard.medium= -log(0.58) / 2  # medium effect
hazard.large = -log(0.65) / 2  # large  effect

# specify parameters from historical data 
shape=0.9509 # shape parameter of the weibull model fitted on historical data

## simulation plan
# Nominal fwer
requiredfwer=0.1  #it is a phase II clinical trial

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

# sample size per arm
sample.sizes <- cbind(
  c( 10,15,20,25,30,40,50),
  c( 10,15,20,25,30,40,50)
)

# run simulations
set.seed(2306)
sim = droptheloser(hazard.exp=s3, n1.perarm=30,n2=30, requiredfwer=requiredfwer, K=K,hazard.null=hazard.null, shape=shape, hazard.censoring = hazard.censoring, followuptime=followuptime, interimanalysisfollowup,  recruitmentrate)
tab <- do.call(rbind, sim[c(2,4)])
tab

# run simulations
set.seed(2306)
sim <- tab <-  list()
tab.sel1 <- tab.sel2 <- matrix(0,nrow(sample.sizes),7) 

for (ss in (1:nrow(sample.sizes))){
for (s in 1:length(scenarios)) {
sim[[s]] = droptheloser(hazard.exp=scenarios[[s]], n1.perarm=sample.sizes[ss,1],n2=sample.sizes[ss,2], requiredfwer=requiredfwer, K=K,hazard.null=hazard.null, shape=shape, hazard.censoring = hazard.censoring, followuptime=followuptime, interimanalysisfollowup,  recruitmentrate)
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
table.power
