# This code upload the functions for the drop-the-losers design for time-to-event outcome using a historical control arm in the misspecified case
# Please use code DTLHC_misspecified_main.R to reproduce results of the simulation study

# A two-stage drop-the-losers design for time-to-event outcome using a historical control arm
# by R. Abbas, J. Wason, S. Michiels, and G. Le Teuff (2020)
# Phase II clinical trials are simulated according to a drop-the-losers and a fixed designs both using historical control arm.

library(mvtnorm)
library(flexsurv)

#function that supports probeventbeforecensoring
integrand.censoring=function(c,mu.c,mu.e)
{
  pgengamma(q =c ,mu=mu.e, Q=-0.4005, sigma = 1.23705 )*dgengamma(x =c ,mu=mu.c, Q=-0.4005, sigma = 1.23705 )
}

#calls integrand.censoring, uses numerical integration to determine the probability of an event occuring before a censoring 
#by time t when the hazards of each are respectively hazard.c and hazard.e
probeventbeforecensoring=function(t,mu.c,mu.e)
{
 if (t==0) {prob=0} 
 else {prob=integrate(integrand.censoring,lower=0,upper=t,mu.c=mu.c,mu.e=mu.e)$value}
 return(prob+pgengamma(q =t ,mu=mu.c, Q=-0.4005, sigma = 1.23705 )*(1-pgengamma(q =t ,mu=mu.c, Q=-0.4005, sigma = 1.23705 )))
}
probeventbeforecensoring=Vectorize(probeventbeforecensoring,vectorize.args="t")

correlation.logrank.new=function(K,n1.perarm=20,n2=20,interimanalysisfollowup=0,recruitmentrate=50,followuptime=2,mu.null=1,mu.true=1,mu.censoring=1)
{ #interimanalysis is done when all stage 1 patients have been recruited and followed up for interimanalysisfollowup time
  enrollmenttimes.stage1=seq(0,length=n1.perarm*K,by=1/recruitmentrate)
  interimanalysistime=enrollmenttimes.stage1[length(enrollmenttimes.stage1)]+interimanalysisfollowup
  prob.event.interim=probeventbeforecensoring(interimanalysistime-enrollmenttimes.stage1,mu.true,mu.censoring)
  exp.event.interim=sum(prob.event.interim[seq(1,n1.perarm*K,by=K)])
  enrollmenttimes.stage2=seq(interimanalysistime,length=n2,by=1/recruitmentrate)
  finalanalysistime=enrollmenttimes.stage2[length(enrollmenttimes.stage2)]+followuptime
  prob.event.final=probeventbeforecensoring(finalanalysistime-c(enrollmenttimes.stage1[seq(1,n1.perarm*K,by=K)],enrollmenttimes.stage2),mu.true,mu.censoring)
  exp.event.final=sum(prob.event.final)
  cov=matrix(c(exp.event.interim,exp.event.interim,exp.event.interim,exp.event.final),2,2)
  return(list(cov=cov,cor=cov2cor(cov)))
}

typeIerrorrate=function(criticalvalue,cov,K,requiredtypeIerrorrate)
{ #test statistic must be better than all others at interim, and then better than critical value at final analysis
  A=matrix(0,K,2*K)
  for(i in 1:(K-1))
  A[1:(K-1),1]=1
  A[1:(K-1),(2:K)]=diag(-1,(K-1))
  A[K,K+1]=1
  mean.transform=as.double(A%*%rep(0,2*K))
  cov.transform=A%*%cov%*%t(A)
  return(K*as.double(pmvnorm(lower=c(rep(0,K-1),criticalvalue),upper=rep(Inf,K),mean=mean.transform,sigma = cov.transform))-requiredtypeIerrorrate)
}

#Arguments as per correlation.logrank, except:
#mu.exp is now a vector with mu for each experimental arm
#requiredfwer is the target family-wise error rate (probability of recommending a treatment when they are all ineffective)
droptheloser=function(seed,
					  K,
                      mu.null, 
                      mu.exp,
                      mu.censoring,
                      hazard.null,followuptime,
                      n1.perarm,
                      n2,
                      requiredfwer,
                      interimanalysisfollowup=0, 
                      recruitmentrate,
                      nit=10000)
{set.seed(seed)
  #get distribution of mean logrank test under global null:
  correlation=correlation.logrank.new(K,n1.perarm,n2,interimanalysisfollowup,recruitmentrate,followuptime,mu.null,mu.null,mu.censoring)$cor[1,2]
  
  #distribution of test statistics under global null:
  
  mean=rep(0,2*K)
  
  cov=matrix(0,2*K,2*K)
  #divide matrix into blocks:
  cov[1:K,1:K]=diag(1,K)
  cov[((K+1):(2*K)),((K+1):(2*K))]=diag(1,K)
  cov[1:K,((K+1):(2*K))]=diag(correlation,K)
  cov[((K+1):(2*K)),1:K]=diag(correlation,K)
  
  #search for critical value that gives required family-wise error rate
  criticalvalue=uniroot(f = typeIerrorrate,lower = -2,upper=5,cov=cov,K=K,requiredtypeIerrorrate=requiredfwer)$root
  
#now use simulation to find power:
  sim=simulation.new(K=K,hazard.null=hazard.null,mu.null=mu.null,mu.exp=mu.exp,mu.censoring=mu.censoring,followuptime=followuptime,n1.perarm=n1.perarm,n2=n2,criticalvalue=criticalvalue,interimanalysisfollowup=interimanalysisfollowup,recruitmentrate=recruitmentrate,niterations=nit)
  rejecth0=sim$r

return(list(#rH0= rejecth0 ,
  efwer  = if( sum(mu.exp == mu.null) > 1 ) mean(ifelse(rowSums(rejecth0[,which (mu.exp == mu.null)]) > 0 , 1,0 ) ) else mean( rejecth0[,which (mu.exp == mu.null)] ),
  disjunctive.power = if (sum(mu.exp != mu.null) == 0 ) print(" - ") else {if( sum(mu.exp != mu.null) > 1 ) mean( ifelse(rowSums(rejecth0[,which (mu.exp != mu.null)]) > 0 , 1,0 ) ) else mean( rejecth0[,which (mu.exp != mu.null)] ) }
  ))   
  
}

simulation.new=function(K,
                        mu.null, 
                        mu.exp,
                        mu.censoring,
                        hazard.null,
                        followuptime,
                        n1.perarm,
                        n2,
                        criticalvalue,
                        interimanalysisfollowup,
                        recruitmentrate,
                        niterations)
{
  selected.arm=stage2.O=stage2.E=stage2.lr1=stage2.lr2=rep(0, niterations)
  rejecth0=stage1.O=stage1.E=stage1.lr1=stage1.lr2=matrix(0, niterations, K)

  for(iteration in 1:niterations)
  {
    enrollmenttimes.stage1=seq(0,length=n1.perarm*K,by=1/recruitmentrate)
    #interim analysis time is after last patient has been followed up for interimanalysisfollowup
    interimanalysistime=enrollmenttimes.stage1[length(enrollmenttimes.stage1)]+interimanalysisfollowup
   
    #for each arm, extract event time and type:
    eventtimes.stage1.perarm=matrix(0,n1.perarm,K)
    censoringtimes.stage1.perarm=matrix(0,n1.perarm,K) 
    censoringtimes.stage1.perarm.interim=matrix(0,n1.perarm,K)
    enrollmenttimes.stage1.perarm=matrix(0,n1.perarm,K)
    time.stage1.perarm=matrix(0,n1.perarm,K)
    event.stage1.perarm=matrix(0,n1.perarm,K)
    for(k in 1:K)
    {
      eventtimes.stage1.perarm[,k]=rgengamma(n =n1.perarm ,mu=mu.exp[k], Q=-0.4005, sigma = 1.23705 )  # here we generate data according to a prespecified distribution
      censoringtimes.stage1.perarm[,k]=rgengamma(n =n1.perarm ,mu=mu.censoring, Q=-0.4005, sigma = 1.23705 ) # here we generate data according to a prespecified distribution
      
      enrollmenttimes.stage1.perarm[,k]=enrollmenttimes.stage1[seq(k,length(enrollmenttimes.stage1),by=K)]
      censoringtimes.stage1.perarm.interim[,k]=replace(censoringtimes.stage1.perarm[,k],which((censoringtimes.stage1.perarm[,k]+enrollmenttimes.stage1.perarm[,k])>interimanalysistime),(interimanalysistime-enrollmenttimes.stage1.perarm[,k])[which((censoringtimes.stage1.perarm[,k]+enrollmenttimes.stage1.perarm[,k])>interimanalysistime)])
    
    #change censoringtimes.stage1 so maximum is interimanalysistime:
    time.stage1.perarm[,k]=ifelse(eventtimes.stage1.perarm[,k]<censoringtimes.stage1.perarm.interim[,k],eventtimes.stage1.perarm[,k],censoringtimes.stage1.perarm.interim[,k])
    event.stage1.perarm[,k]=ifelse(eventtimes.stage1.perarm[,k]<censoringtimes.stage1.perarm.interim[,k],1,0)
    }
    
    # the observed and expected event case numbers
    stage1.O[iteration, ]= colSums(event.stage1.perarm )
    stage1.E[iteration, ]= colSums(hazard.null*time.stage1.perarm) # exponential

    # the "modified" one-sample log-rank test
    stage1.lr2[iteration, ]= -(stage1.O[iteration, ]-stage1.E[iteration, ])/sqrt((stage1.O[iteration, ]+stage1.E[iteration, ])/2)

    #select treatment with maximum test statistic:
    selectedarm=which.max(stage1.lr2[iteration, ])

    enrollmenttimes.stage2=seq(interimanalysistime,length=n2,by=1/recruitmentrate)
    finalanalysistime=enrollmenttimes.stage2[length(enrollmenttimes.stage2)]+followuptime
    eventtimes.stage2=rgengamma(n =n2 ,mu=mu.exp[selectedarm], Q=-0.4005, sigma = 1.23705 )
    censoringtimes.stage2=rgengamma(n =n2 ,mu=mu.censoring, Q=-0.4005, sigma = 1.23705 )
    censoringtimes.stage2=replace(censoringtimes.stage2,which((censoringtimes.stage2+enrollmenttimes.stage2)>finalanalysistime),(finalanalysistime-enrollmenttimes.stage2)[which((censoringtimes.stage2+enrollmenttimes.stage2)>finalanalysistime)])
    #get censoring times of first stage participants on selected arm:
    
    censoringtimes.stage1.final=replace(censoringtimes.stage1.perarm[,selectedarm],which((censoringtimes.stage1.perarm[,selectedarm]+enrollmenttimes.stage1.perarm[,selectedarm])>finalanalysistime),(finalanalysistime-enrollmenttimes.stage1.perarm[,selectedarm])[which((censoringtimes.stage1.perarm[,selectedarm]+enrollmenttimes.stage1.perarm[,selectedarm])>finalanalysistime)])
        
    time.stage1=ifelse(eventtimes.stage1.perarm[,selectedarm]<censoringtimes.stage1.final,eventtimes.stage1.perarm[,selectedarm],censoringtimes.stage1.final)
    event.stage1=ifelse(eventtimes.stage1.perarm[,selectedarm]<censoringtimes.stage1.final,1,0)
    time.stage2=ifelse(eventtimes.stage2<censoringtimes.stage2,eventtimes.stage2,censoringtimes.stage2)
    event.stage2=ifelse(eventtimes.stage2<censoringtimes.stage2,1,0)
    time=c(time.stage1,time.stage2)
    event=c(event.stage1,event.stage2)

    # the observed and expected event case numbers
    stage2.O[iteration]=sum(event )
    stage2.E[iteration]= sum(hazard.null*time) # exponential
 
    # the "modified" one-sample log-rank test
    stage2.lr2[iteration]= -(stage2.O[iteration]-stage2.E[iteration])/sqrt((stage2.O[iteration]+stage2.E[iteration])/2)

  selected.arm[iteration]=selectedarm
  rejecth0[iteration, selectedarm]=ifelse(stage2.lr2[iteration]>criticalvalue, 1, 0)  
  }
  return( list(r=rejecth0, sel=selected.arm))
}

find.parameter= function(S,t,m,q,s)( (1-pgengamma(t, mu=m, Q=q, sigma = s)))-S
