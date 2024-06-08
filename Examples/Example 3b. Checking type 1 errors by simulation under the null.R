###################################################################################
#  Purpose: Generate data under the null hypothesis.
#           Then apply the methods on this to verify type 1 error rates etc
#
#  Programmer: 
#  Date:  
###################################################################################

#libraries----------------------------------------------------------------------

library(BITsurv)  
#A lot of the machinery from BITsurv is explicitly written out as code here
#This is because we are constraining ourselve to test the null

library(dplyr)

library(flexsurv)
library(survminer)
library(sinib)






#Some required functions--------------------------------------------------------

#specified interval assignment functions----------------------------------------
lower.int<-function(time){ifelse(time<=min(spec_int),
                                 0,           #forces the lowest interval to have a lower bound of 0
                                 max(spec_int[spec_int<time]))}
lower.int.vec<-function(time){sapply(time, lower.int)}

upper.int<-function(time){ifelse(time>max(spec_int),Inf,   #forces the highest interval to have an upper bound of infinity
                                 min(spec_int[spec_int>=time]))}
upper.int.vec<-function(time){sapply(time, upper.int)}


#specified+censor interval assignment functions---------------------------------
lower.sub.int<-function(time){ifelse(time<=min(sub.ints),
                                     0,   #forces the lowest subinterval to have a lower bound of 0
                                     max(sub.ints[sub.ints<time]))}
lower.sub.int.vec<-function(time){sapply(time, lower.sub.int)}


upper.sub.int<-function(time){ifelse(time>max(sub.ints),
                                     Inf, #forces the highest subinterval to have an upper bound of infinity
                                     min(sub.ints[sub.ints>=time]))}
upper.sub.int.vec<-function(time){sapply(time, upper.sub.int)}











#The model for the simulated data-----------------------------------------------

#Event generating process is
# T~exp(lambda=1/10)

#Censor generating process is
# C1~uniform[0,100]      #general censoring
# C2~uniform[18,22]      #heavy censoring at the end of the curve
# C<-min(C1,C2)          #final censor value    

#100 patients

#times=min(C,T)
#events=1 if T<=C, or 1 otherwise






#Simulation specifications------------------------------------------------------

#number of patients
N.patients<-200    

#expected survival in underlying model = 1/rate.spec
rate.spec<-1/10

#number of simulations
N.sims<-10000

#Interval selection approach
Approach<-'Ten.Fixed.Ints'
#or
#Approach<-'Censor.Ints'





################################################################################
# Run the simulation------------------------------------------------------------
################################################################################
for (i in 1:N.sims){
#i<-1
  
  surv.data<- data.frame(Time.of.Event=rexp(N.patients,rate=rate.spec),
                       C1=runif(N.patients, min=0, max=100),
                       C2=runif(N.patients, min=18, max=22)) %>% 
  mutate(C=pmin(C1,C2)) %>% 
  mutate(event=ifelse(Time.of.Event<C,1,0)) %>% 
  mutate(time=pmin(C,Time.of.Event))%>% 
  select(time, event) 

  
#Decides which specified interval approach to use

censors<-surv.data %>% filter(event==0)

if (Approach=='Ten.Fixed.Ints'){
spec_int<-0.1*max(censors$time)*0:10  
} else{
  spec_int<-censors$time
}
  
#Start of extracted BIT.surv function  
#-------------------------------------------------------------------------------
#The following is just the BIT.surv function for exp----------------------------
#where a small edit has been made to use the true event generating process------
#when calculating p. This edit is highlighted in the code------------------------
#-------------------------------------------------------------------------------
  
  #Grouping and formatting the data-----------------------------------------------
  

  new.data<-rbind(surv.data %>% select(time,event),
                  data.frame(time=spec_int,
                             event=2)) %>% 
    arrange(time) %>% 
    mutate(V.lower=lower.int.vec(time),
           V.upper=upper.int.vec(time))
  
  #define the sub intervals I
  sub.ints<-new.data %>% filter(event!=1)
  sub.ints<-unique(sub.ints$time)
  
  new.data<-new.data %>% 
    filter(event!=2) %>%    #remove the specified intervals as these are no longer of interest
    mutate(I.lower=lower.sub.int.vec(time),
           I.upper=upper.sub.int.vec(time)) %>% 
    mutate(I.lower=ifelse(I.lower==-Inf,0,I.lower))
  
  #So far this is just appropriate grouping of V and I,
  #which can all be checked by reviewing the df
  
  #Obtaining statistics for each interval-----------------------------------------
  
  #We now want the number of observed events in each interval I
  new.data.1<-new.data %>% 
    group_by(I.lower) %>% 
    summarize(I.upper=unique(I.upper),
              V.lower=unique(V.lower),
              V.upper=unique(V.upper),
              Events.obs.I=sum(event),
              Pats.in.I=length(event)) 
  
  #Now we want to calculate the number at risk at the start of interval I
  new.data.1<-new.data.1%>% 
    mutate(N.risk=sum(new.data.1$Pats.in.I)-                      #total number of patients  -
             cumsum(c(0,                                      #patients lost in previous intervals
                      Pats.in.I[1:(nrow(new.data.1)-1)])) ) 
  
  
  
  # Calculating p_{I_j} values----------------------------------------------------
  
  #We also want to calculate the value p for interval I 
  #where p=P(patient experiences an event in I|they make it to I)
  #We define p based on the fitted model that is being used
  
#Edited**************************************************************************
#This has been edited to use the null to calculate p
  
    new.data.2<-new.data.1%>% 
      mutate(p={pexp(q=I.upper, rate=rate.spec)-       #calculate of p for fitted exponential model
          pexp(q=I.lower, rate=rate.spec) }/
            (1-pexp(q=I.lower,rate=rate.spec))) 
# End of edit*******************************************************************  
  
  
  
  #Final summaries----------------------------------------------------------------
  
  #Now we look at summarizing the specified intervals
  
  new.data.3<-new.data.2 %>%
    mutate(expected_E=p*N.risk)%>%    #first add the expect number of events under the fitted model
    group_by(V.lower) %>%
    summarise(V.upper=unique(V.upper),
              Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
              Observed.E_over.V=sum(Events.obs.I),
              V.mid.pval=if(length(V.lower)==1){   #there are cases where sinib breaks in the binomial case
                #to prevents these breaks (p.values outside (0,1)) we use the classical pbinom function
                
                0.5*pbinom(q=as.integer(sum(Events.obs.I)),    #p-value calculation for each V
                           size=as.integer(N.risk),
                           prob=p)+
                  0.5*pbinom(q=as.integer(sum(Events.obs.I)-1),
                             size=as.integer(N.risk),
                             prob=p)
                
              }else{ #for non-binomial situations sinib is used
                
                0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)),    #p-value calculation for each V
                                  size=as.integer(N.risk),
                                  prob=p) +
                  ifelse(sum(Events.obs.I)==0,       #additional step as psinib can break if q<0
                         0,
                         0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)-1),
                                           size=as.integer(N.risk),
                                           prob=p))
              },
              N.risk.at.V.start=max(N.risk),
              E=sum(Events.obs.I))
  
  #In the case that the largest time is an event and the specified V are such that
  #the last event is outside of V then you obtain an interval with an upper bound of
  #resulting in a Prob=1 of observing an event within said interval, making the observed
  #event redundant. As such, such an interval is not of interest and is removed.
  
  #Similarly, if for whatever reason V is specified such that certain events are
  #outside of this. Then they will be assigned to this infinite interval and similarly
  #should be dropped. This is done here
  new.data.3<-new.data.3 %>% filter(V.upper!=Inf)
  
  
  
  #Error catching-------------------------------------------------------------------
  
  #There are a few rare situations (~1%) where sinib does not converge
  #and gives p-values outside of [0,1].
  #We will now implement a step to protect against this
  
  #First check new.data.3 for any p values outside of [0,1]
  
  p.checks<-new.data.3$V.mid.pval
  if(sum(!between(p.checks, left=0, right=1))>=1){ #if one or more p-values fall
    #outside of 0 and 1 then
    #find the v.lower values corresponding to these errors
    error.intervals<-new.data.3$V.lower[!between(p.checks, left=0, right=1)]
    
    #now we want to obtain the corrected values
    #we do this be approximating the probabilities by using the round function
    #this should hopefully provide values that have converged
    ammended.p<-new.data.2 %>%
      filter(V.lower %in% error.intervals) %>%   #select only data related to the error intervals
      mutate(expected_E=p*N.risk)%>%
      group_by(V.lower) %>%
      summarise(V.upper=unique(V.upper),
                Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
                Observed.E_over.V=sum(Events.obs.I),
                V.mid.pval= 0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)),    #p-value calculation for each V
                                              size=as.integer(N.risk),
                                              prob=round(p, digits = 3)) +##################
                ifelse(sum(Events.obs.I)==0,       #additional step as psinib can break if q<0
                       0,
                       0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)-1),
                                         size=as.integer(N.risk),
                                         prob=round(p, digits = 3))),############
                N.risk.at.V.start=max(N.risk),
                E=sum(Events.obs.I) )
    
    
    #Now to replace the error p.val with this updated one
    
    #values outside of range are given by new.data.3$V.mid.pval[which(!between(p.checks, left=0, right=1))]
    #and are updated to ammendment
    new.data.3$V.mid.pval[which(!between(p.checks, left=0, right=1))]<-ammended.p$V.mid.pval
  }
  #End of error catching
  
  
  #Bonferroni and individual test results-----------------------------------------
  #for the specified intervals
  
  information<-new.data.3 %>% 
    mutate(individual.test=ifelse(V.mid.pval<=0.025| V.mid.pval>=0.975,
                                  'Reject', 'Accept')) %>% 
    mutate(bonferroni.test=ifelse(V.mid.pval<=         #I need to derive this 2-tailed result as it is not in my lecture notes
                                    0.025/nrow(new.data.3) |
                                    V.mid.pval>=
                                    1- 0.025/nrow(new.data.3),
                                  'Reject',
                                  'Accept')) %>% 
    select(-E) 

#-------------------------------------------------------------------------------  
#End of extracted BIT.surv function  
#-------------------------------------------------------------------------------  
  
  
#Save results-------------------------------------------------------------------  
  
#Use the information  df to save results of interest for that simulation
  
#We are interested in:
#Bonferroni rejects (that under H0 Prob bonferroni>=1 is approxiately 0.05 or less)  
#TS.PAVSI (should approximately reject 0.05 times) 
#TS.TFT (should approximately reject 0.05 times)
  
sim.results<-data.frame(Bonferroni.rejects=sum(information$bonferroni.test=='Reject'),
                        Individual.rejects=sum(information$individual.test=='Reject'),
                        PAVSI.TS=BIT.TS.PAVSI(information$V.mid.pval,print=FALSE),      
                        TFT.TS=BIT.TS.TFT(information$V.mid.pval, print=FALSE),
                        Intervals=nrow(information)) 
  
#save as final results
if(i==1){
 final.results<-sim.results 
}else{
  final.results<-rbind(final.results,
                       sim.results)
}

  
  
} #end loop of simulations












################################################################################
#Analysis of final simulation results--------------------------------------------
################################################################################
#final.results


#Bonferroni results-------------------------------------------------------------

#We are interested in:
#Bonferroni rejects (that under H0 Prob bonferroni>=1 is approxiately 0.05 or less)  
sum(final.results$Bonferroni.rejects>0)/N.sims
#For Ten.Fixed.Ints' approach with 10,000 sims and 200 patients we get
#0.0435

#For censor intervals with 10,000 sims and 200 patients we get
#0.0213



#PAVSI results------------------------------------------------------------------

#TS.PAVSI (should approximately reject 0.05 times) 
sum(final.results$PAVSI.TS<=0.05)/N.sims
#For Ten.Fixed.Ints' approach with 10,000 sims and 200 patients we get
#0.0700
#0.0732

#For censor intervals with 10,000 sims and 200 patients we get
#0.0041





#TFT results-------------------------------------------------------------

#TS.TFT (should approximately reject 0.05 times)
sum(final.results$TFT.TS<=0.05,rm.na=TRUE)/N.sims
#This often reports NA.
#This is due to a small number of errors in sinib
#Note that numerous steps have been implemented to reduce the prevailence of
#these errors down to around 0.1% 


#Let's analyse the TFT results with this in mind

#How many sinib failures?
#Still some NaN values
sum(final.results$TFT.TS=='NaN')/N.sims
#For Ten.Fixed.Ints' approach with 10,000 sims and 200 patients
#around 0.0024
#For censor intervals with 10,000 sims and 200 patients we get
#around 0.0002

#What is the type 1 error rate when removing sinib failures?
sum(final.results$TFT.TS[final.results$TFT.TS!='NaN']<=0.05)/N.sims
#removing NaN values we get....
#For Ten.Fixed.Ints' approach with 10,000 sims and 200 patients we get
#0.0397

#For censor intervals with 10,000 sims and 200 patients we get
#0.0001




#How are the these TS distributed?
#should also be approximately uniform [0,1]
hist(final.results$TFT.TS, breaks=20, probability = TRUE)
#For Ten.Fixed.Ints' approach with 10,000 sims and 200 patients...
#this is approximately uniform with a very slight upwards skew
#density at x-axis=0 is around 0.8, then from x-axis=0.4 to 1 has a density around 1.1

#For censor intervals with 10,000 sims and 200 patients...
#this is highly skewed upwards with over 95% of the density coming in above
#0.9. That is,  this is not uniform[0,1] and gives a very low type 1 error rate under the null

#As a final note, if you encounter a rare sinib error in your base analysis (not a
#simulation experiment) where my error catching has not resolved this, you can 
#get around this by approximating the probabilities further to give an approximate 
#value for the offending interval (as I have done in the error catching step)
#say run the error catching step with prob=round(p, digits = 2) (instead of 3)