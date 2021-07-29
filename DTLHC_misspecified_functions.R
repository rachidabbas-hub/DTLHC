# This code upload the functions for the drop-the-losers design for time-to-event outcome using a historical control arm in the misspecified case
# Please use code DTLHC_misspecified_main.R to reproduce results of the simulation study

# A two-stage drop-the-losers design for time-to-event outcome using a historical control arm
# by R. Abbas, J. Wason, S. Michiels, and G. Le Teuff (2021)
# Phase II clinical trials are simulated according to a drop-the-losers and a fixed designs both using historical control arm.

library(mvtnorm)
library(flexsurv)

#function that supports probeventbeforecensoring
integrand.censoring=function(c,mu.c,mu.e)
{
  pgengamma(q =c ,mu=mu.e, Q=-0.4005, sigma = 1.23705 )*dgengamma(x =c ,mu=mu.c, Q=-0.4005, sigma = 1.23705 )
}

#calls integrand.censoring, uses numerical integration to determine the probability of an event occuring before a censoring 
probeventbeforecensoring=function(t,mu.c,mu.e)
{
  if (t==0) {prob=0} 
  else {prob=integrate(integrand.censoring,lower=0,upper=t,mu.c=mu.c,mu.e=mu.e)$value}
  return(prob+pgengamma(q =t ,mu=mu.e, Q=-0.4005, sigma = 1.23705 )*(1-pgengamma(q =t ,mu=mu.c, Q=-0.4005, sigma = 1.23705 )))
}
probeventbeforecensoring=Vectorize(probeventbeforecensoring,vectorize.args="t")


#correlation.logrank.new finds the correlation between the interim and final log-rank test taking into account
#follow-up time and recruitment rate
#Arguments:
#K - number of experimental arms
#n1.perarm: number of patients recruited in first stage per arm
#n2: number of patients recruited per arm in the second stage
#interimanalysisfollowup: minimum amount of time that first stage patients are followed up for before interim analysis
#recruitmentrate: number of patients per year recruited
#recruitment rate - number of patients recruited per year
#followuptime: minimum amount of time that patients are followed for
#mu.null: event mu parameter from gamma distribution under the null hypothesis
#mu.true: true event mu parameter from gamma distribution
#mu.censoring: true censoring mu parameter from gamma distribution

#Returns covariance and correlation matrices as a list

correlation.logrank.gamma=function(K,n1.perarm=20,n2=20,interimanalysisfollowup=0,recruitmentrate=50,followuptime=2,mu.null=1,mu.true=1,mu.censoring=1)
{
  #interimanalysis is done when all stage 1 patients have been recruited and followed up for interimanalysisfollowup time
 
  enrollmenttimes.stage1=seq(0,length=n1.perarm*K,by=1/recruitmentrate)
  interimanalysistime=enrollmenttimes.stage1[length(enrollmenttimes.stage1)]+interimanalysisfollowup
  prob.event.interim=probeventbeforecensoring(interimanalysistime-enrollmenttimes.stage1,mu.e=mu.true,mu.c=mu.censoring)
  exp.event.interim=sum(prob.event.interim[seq(1,n1.perarm*K,by=K)])
  enrollmenttimes.stage2=seq(interimanalysistime,length=n2,by=1/recruitmentrate)
  finalanalysistime=enrollmenttimes.stage2[length(enrollmenttimes.stage2)]+followuptime
  prob.event.final=probeventbeforecensoring(finalanalysistime-c(enrollmenttimes.stage1[seq(1,n1.perarm*K,by=K)],enrollmenttimes.stage2),mu.e=mu.true,mu.c=mu.censoring)
  exp.event.final=sum(prob.event.final)
  
  # do the same under H0
  prob.event.interim.H0=probeventbeforecensoring(interimanalysistime-enrollmenttimes.stage1, mu.e=mu.null, mu.c=mu.censoring)
  exp.event.interim.H0=sum(prob.event.interim.H0[seq(1, n1.perarm*K, by=K)])
  prob.event.final.H0=probeventbeforecensoring(finalanalysistime-c(enrollmenttimes.stage1[seq(1, n1.perarm*K, by=K)], enrollmenttimes.stage2), mu.e=mu.null, mu.c=mu.censoring)
  exp.event.final.H0=sum(prob.event.final.H0)

# takes into account the expected number of events as per Wu's paper by averaging under H1 and H0
  Ei=(exp.event.interim + exp.event.interim.H0)/2
  Ef=(exp.event.final + exp.event.final.H0)/2
  cov=matrix(c(Ei, Ei, Ei, Ef), 2, 2)
  
  return(list(cov=cov, cor=cov2cor(cov),sqrtcov=sqrt(cov) ))
}

typeIerrorrate=function(criticalvalue, cov, K, requiredtypeIerrorrate)
{
 #find prob of rejecting hypothesis 1 similar in notation to Wason et al,  SMMR 2017:
 #test statistic must be better than all others at interim,  and then better than critical value at end
  A=matrix(0, K, 2*K)
  A[1,1]=A[2,2]=1
  A[1,2]=A[2,3]=-1
  A[K, K+1]=1
  mean.transform=as.double(A%*%rep(0, 2*K))
  cov.transform=A%*%cov%*%t(A)
return(factorial(K)*as.double(pmvnorm(lower=c(rep(0, K-1), criticalvalue), upper=rep(Inf, K), mean=mean.transform, sigma = cov.transform))-requiredtypeIerrorrate)
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
{
set.seed(seed)

  #get distribution of mean logrank test under global null:
  correlation=correlation.logrank.gamma(K,n1.perarm,n2,interimanalysisfollowup,recruitmentrate,followuptime,mu.null=mu.null,mu.true=mu.null,mu.censoring=mu.censoring)
  
  #distribution of test statistics under global null:
  mean=rep(0,2*K)
  
  cov=matrix(0, 2*K, 2*K)
  cov[1:K, 1:K]=diag(1, K)
  cov[((K+1):(2*K)), ((K+1):(2*K))]=diag(1, K)
  cov[1:K, ((K+1):(2*K))]=diag(correlation$cor[1, 2], K)
  cov[((K+1):(2*K)), 1:K]=diag(correlation$cor[1, 2], K)
  
  #search for critical value that gives required family-wise error rate
  criticalvalue=uniroot(f = typeIerrorrate,lower = -2,upper=5,cov=cov,K=K,requiredtypeIerrorrate=requiredfwer)$root
# Get the theoretic power

  #Figure out the scenario case
  if (sum(mu.exp!=mu.null) == 0) theoric.d.power = 0
else{
##get distribution of mean logrank test under alternative

# Distribution of test statistics under alternative
  cov1=correlation.logrank.gamma(K, n1.perarm, n2, interimanalysisfollowup, recruitmentrate, followuptime, mu.null=mu.null, mu.true=mu.exp[1], mu.censoring=mu.censoring)
  cov2=correlation.logrank.gamma(K, n1.perarm, n2, interimanalysisfollowup, recruitmentrate, followuptime, mu.null=mu.null, mu.true=mu.exp[2], mu.censoring=mu.censoring)
  cov3=correlation.logrank.gamma(K, n1.perarm, n2, interimanalysisfollowup, recruitmentrate, followuptime, mu.null=mu.null, mu.true=mu.exp[3], mu.censoring=mu.censoring)

# it is minus because the OSLRT is reversed: E-O
H.exp =-log(pgengamma(q =2 ,mu=mu.exp,  Q=-0.4005, sigma = 1.23705 ))/2
H.null=-log(pgengamma(q =2 ,mu=mu.null, Q=-0.4005, sigma = 1.23705 ))/2
  mean.alt = rep( log(H.exp/H.null), 2)*c(cov1$sqrtcov[1,2],cov2$sqrtcov[1,2],cov3$sqrtcov[1,2],cov1$sqrtcov[2,2],cov2$sqrtcov[2,2],cov3$sqrtcov[2,2])
  cor.alt = rbind(cov1$cor,cov2$cor,cov3$cor)[seq(2, 2*K, by=2),1]
 
  cov.alt=matrix(0, 2*K, 2*K)
  cov.alt[1:K, 1:K]=diag(1, K)
  cov.alt[((K+1):(2*K)), ((K+1):(2*K))]=diag(1, K)
  cov.alt[1:K, ((K+1):(2*K))]=diag(cor.alt, K)
  cov.alt[((K+1):(2*K)), 1:K]=diag(cor.alt, K)
  
  #A is a matrix used to transform the vector of normal test statistics (partial ordering)
  A1=A2=A3=matrix(0, K, 2*K)
  
  # A1: Treatment 1 is better than all others at interim,  and then better than critical value at the end (partial ordering)
  A1[1,1]=A1[2,1]=A1[K, K+1]=1
  A1[1,2]=A1[2,3]=-1
  
  # A2: Treatment 2 is better than all others at interim,  and then better than critical value at the end (partial ordering)
  A2[1,2]=A2[2,2]=A2[K, K+2]=1
  A2[1,1]=A2[2,3]=-1
  
  # A3: Treatment 3 is better than all others at interim,  and then better than critical value at the end (partial ordering)
  A3[1,3]=A3[2,3]=A3[K, K+3]=1
  A3[1,1]=A3[2,2]=-1
  
### compute all base case probability
  p1 = pmvnorm(lower=c(rep(0, K-1), criticalvalue), upper=rep(Inf, K), mean=as.numeric(A1%*%mean.alt), sigma = A1%*%cov.alt%*%t(A1)) 
  p2 = pmvnorm(lower=c(rep(0, K-1), criticalvalue), upper=rep(Inf, K), mean=as.numeric(A2%*%mean.alt), sigma = A2%*%cov.alt%*%t(A2)) 
  p3 = pmvnorm(lower=c(rep(0, K-1), criticalvalue), upper=rep(Inf, K), mean=as.numeric(A3%*%mean.alt), sigma = A3%*%cov.alt%*%t(A3)) 
  
### sum-up base case probability according to the secnario
  if(sum(mu.exp!=mu.null) == 1) {theoric.d.power = as.numeric(p1) }
  else {
    if(sum(mu.exp!=mu.null) == 2) {theoric.d.power = as.numeric(p1 + p2)} 
    else {
      if(sum(mu.exp!=mu.null) == 3) {theoric.d.power = as.numeric(p1 + p2 + p3)  }
    }
  }
}
  
# now use simulation to find empirical power:
  sim=simulation.new(K=K,hazard.null=hazard.null,mu.null=mu.null,mu.exp=mu.exp,mu.censoring=mu.censoring,followuptime=followuptime,n1.perarm=n1.perarm,n2=n2,criticalvalue=criticalvalue,interimanalysisfollowup=interimanalysisfollowup,recruitmentrate=recruitmentrate,niterations=nit)
  rejecth0=sim$r

return(list(rH0= rejecth0 ,
  efwer  = if( sum(mu.exp == mu.null) > 1 ) mean(ifelse(rowSums(rejecth0[,which (mu.exp == mu.null)]) > 0 , 1,0 ) ) else mean( rejecth0[,which (mu.exp == mu.null)] ),
  disjunctive.power = if (sum(mu.exp != mu.null) == 0 ) print(" - ") else {if( sum(mu.exp != mu.null) > 1 ) mean( ifelse(rowSums(rejecth0[,which (mu.exp != mu.null)]) > 0 , 1,0 ) ) else mean( rejecth0[,which (mu.exp != mu.null)] ) },
  selectedarm=sim$sel,
  theoric.d.power=as.numeric(theoric.d.power)
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
  selected.arm=stage2.O=stage2.E=stage2.lr2=rep(0, niterations)
  rejecth0=stage1.O=stage1.E=stage1.lr2=matrix(0, niterations, K)

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
