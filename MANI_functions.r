# This code upload the functions for the fixed design with or without correction for multiplicity
# Please use code MANI_main.R to reproduce results of the simulation study

# A two-stage drop-the-losers design for time-to-event outcome using a historical control arm
# by R. Abbas, J. Wason, S. Michiels, and G. Le Teuff (2020)
# Phase II clinical trials are simulated according to a drop-the-losers and a fixed designs both using historical control arm.


mani=function(seed,K, hazard.null, wei_shape, hazard.exp, hazard.censoring, followuptime, n.total, alpha, recruitmentrate, niterations)
{
  set.seed(seed)
  n.perarm = floor(n.total /K)
  # use simulation to find power:
  rejecth0=Observed=Expected=lr2=matrix(0, niterations, K)
  for(iteration in 1:niterations)
  {
    enrollmenttimes=seq(0, length=n.perarm*K, by=1/recruitmentrate)
    finalanalysistime=enrollmenttimes[length(enrollmenttimes)]+followuptime
    
    #for each arm,  extract event time and type:
    eventtimes.perarm=censoringtimes.perarm=censoringtimes.perarm.final=enrollmenttimes.perarm=time.perarm=event.perarm=matrix(0,n.perarm, K)
    for(k in 1:K)
    { 
      eventtimes.perarm[, k]=rexp(n.perarm, hazard.exp[k])
      censoringtimes.perarm[, k]=rexp(n.perarm, hazard.censoring)
      enrollmenttimes.perarm[, k]=enrollmenttimes[seq(k, length(enrollmenttimes), by=K)]
      #change censoringtimes so maximum is finalanalysistime:
      censoringtimes.perarm.final[, k]=replace(censoringtimes.perarm[, k], which((censoringtimes.perarm[, k]+enrollmenttimes.perarm[, k])>finalanalysistime), (finalanalysistime-enrollmenttimes.perarm[, k])[which((censoringtimes.perarm[, k]+enrollmenttimes.perarm[, k])>finalanalysistime)])
      
      time.perarm[, k]=ifelse(eventtimes.perarm[, k]<censoringtimes.perarm.final[, k], eventtimes.perarm[, k], censoringtimes.perarm.final[, k])
      event.perarm[, k]=ifelse(eventtimes.perarm[, k]<censoringtimes.perarm.final[, k], 1, 0)
    }
    
    # the observed and expected event case numbers
    Observed[iteration, ]= colSums(event.perarm )
    Expected[iteration, ]= colSums(hazard.null*time.perarm^wei_shape) # weibull
    
    lr2[iteration, ]= -(Observed[iteration, ]-Expected[iteration, ])/sqrt((Observed[iteration, ]+Expected[iteration, ])/2)
    rejecth0[iteration, ]=ifelse(lr2[iteration,]>qnorm(1-alpha), 1, 0)
  }
  
return(list( #O=Observed, E=Expected,  
  efwer  = if( sum(hazard.exp == hazard.null) > 1 ) mean(ifelse(rowSums(rejecth0[, which (hazard.exp == hazard.null)]) > 0 ,  1, 0 ) ) else mean( rejecth0[, which (hazard.exp == hazard.null)] ), 
  conjunctive.power  = if (sum(hazard.exp != hazard.null) == 0 ) print(" - ") else {if( sum(hazard.exp != hazard.null) > 1 ) mean( ifelse(rowSums(rejecth0[, which (hazard.exp != hazard.null)]) == length( which (hazard.exp != hazard.null) ),  1, 0 ) ) else mean( rejecth0[, which (hazard.exp != hazard.null)]) }, 
  disjunctive.power = if (sum(hazard.exp != hazard.null) == 0 ) print(" - ") else {if( sum(hazard.exp != hazard.null) > 1 ) mean( ifelse(rowSums(rejecth0[, which (hazard.exp != hazard.null)]) > 0 ,  1, 0 ) ) else mean( rejecth0[, which (hazard.exp != hazard.null)] ) }
  ))   
}


mani_hommel=function(seed,K, hazard.null, wei_shape, hazard.exp, hazard.censoring, followuptime, n.total, alpha, recruitmentrate, niterations)
{
  set.seed(seed)
  n.perarm = floor(n.total /K)
  # use simulation to find power:
  praw=Observed=Expected=lr2=matrix(0, niterations, K)
  for(iteration in 1:niterations)
  {
    enrollmenttimes=seq(0, length=n.perarm*K, by=1/recruitmentrate)
    finalanalysistime=enrollmenttimes[length(enrollmenttimes)]+followuptime
    
    #for each arm,  extract event time and type:
    eventtimes.perarm=censoringtimes.perarm=censoringtimes.perarm.final=enrollmenttimes.perarm=time.perarm=event.perarm=matrix(0,n.perarm, K)
    for(k in 1:K)
    { 
      eventtimes.perarm[, k]=rexp(n.perarm, hazard.exp[k])
      censoringtimes.perarm[, k]=rexp(n.perarm, hazard.censoring)
      enrollmenttimes.perarm[, k]=enrollmenttimes[seq(k, length(enrollmenttimes), by=K)]
      #change censoringtimes so maximum is finalanalysistime:
      censoringtimes.perarm.final[, k]=replace(censoringtimes.perarm[, k], which((censoringtimes.perarm[, k]+enrollmenttimes.perarm[, k])>finalanalysistime), (finalanalysistime-enrollmenttimes.perarm[, k])[which((censoringtimes.perarm[, k]+enrollmenttimes.perarm[, k])>finalanalysistime)])
      
      time.perarm[, k]=ifelse(eventtimes.perarm[, k]<censoringtimes.perarm.final[, k], eventtimes.perarm[, k], censoringtimes.perarm.final[, k])
      event.perarm[, k]=ifelse(eventtimes.perarm[, k]<censoringtimes.perarm.final[, k], 1, 0)
    }
    
    # the observed and expected event case numbers
    Observed[iteration, ]= colSums(event.perarm )
    Expected[iteration, ]= colSums(hazard.null*time.perarm^wei_shape) # weibull
    
    lr2[iteration, ]= -(Observed[iteration, ]-Expected[iteration, ])/sqrt((Observed[iteration, ]+Expected[iteration, ])/2)
    praw[iteration,] = 1-pnorm(lr2[iteration,])
  }
  phom = t(apply(praw,1,p.adjust,method="hommel"))	
  rejecth0 = ifelse(phom < alpha , 1, 0)	
  
  return(list( #O=Observed, E=Expected,  
    efwer  = if( sum(hazard.exp == hazard.null) > 1 ) mean(ifelse(rowSums(rejecth0[, which (hazard.exp == hazard.null)]) > 0 ,  1, 0 ) ) else mean( rejecth0[, which (hazard.exp == hazard.null)] ), 
    conjunctive.power  = if (sum(hazard.exp != hazard.null) == 0 ) print(" - ") else {if( sum(hazard.exp != hazard.null) > 1 ) mean( ifelse(rowSums(rejecth0[, which (hazard.exp != hazard.null)]) == length( which (hazard.exp != hazard.null) ),  1, 0 ) ) else mean( rejecth0[, which (hazard.exp != hazard.null)]) }, 
    disjunctive.power = if (sum(hazard.exp != hazard.null) == 0 ) print(" - ") else {if( sum(hazard.exp != hazard.null) > 1 ) mean( ifelse(rowSums(rejecth0[, which (hazard.exp != hazard.null)]) > 0 ,  1, 0 ) ) else mean( rejecth0[, which (hazard.exp != hazard.null)] ) }
  ))   
}
