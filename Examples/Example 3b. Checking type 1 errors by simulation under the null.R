###################################################################################
#  Purpose: Generate data under the null hypothesis.
#           Then apply the methods on this to verify type 1 error rates etc
#
#  Coder: Ben Lee
#  Date:  14-05-2024
###################################################################################


#libraries----------------------------------------------------------------------
library(survHE) #Load survHE. Loads flexsurv and survival automatically
library(survminer) #Required for ggsurvplot
library(dplyr)

library(sinib)

base.file<-"C:/R/Simulated KM curves/Final repo"

#source TS functions
source(file.path(base.file, "Functions", "3. Test statistics.R"))

#Some functions-----------------------------------------------------------------

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

N.patients<-1000

#number of simulations
N.sims<-1000

#Interval selection approach
Approach<-'Ten.Fixed.Ints'
#or
#Approach<-'Censor.Ints'

for (i in 1:N.sims){
#i<-1
  
  surv.data<- data.frame(Time.of.Event=rexp(N.patients,rate=1/10),
                       C1=runif(N.patients, min=0, max=100),
                       C2=runif(N.patients, min=18, max=22)) %>% 
  mutate(C=pmin(C1,C2)) %>% 
  mutate(event=ifelse(Time.of.Event<C,1,0)) %>% 
  mutate(time=pmin(C,Time.of.Event))%>% 
  select(time, event) 

  
  
  
#Decide of specified interval approach  

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
      mutate(p={pexp(q=I.upper, rate=0.1)-       #calculate of p for fitted exponential model
          pexp(q=I.lower, rate=0.1) }/
            (1-pexp(q=I.lower,rate=0.1))) 
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
                
                0.5*psinib(q=as.integer(sum(Events.obs.I)),    #p-value calculation for each V     
                           size=as.integer(N.risk), 
                           prob=p) +
                  ifelse(sum(Events.obs.I)==0,       #additional step as psinib can break if q<0
                         0,
                         0.5*psinib(q=as.integer(sum(Events.obs.I)-1),         
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












#Analysis of final simulation results--------------------------------------------

#final.results


#We are interested in:
#Bonferroni rejects (that under H0 Prob bonferroni>=1 is approxiately 0.05 or less)  
sum(final.results$Bonferroni.rejects>0)/N.sims
#For Ten.Fixed.Ints' approach with 10,000 sims we get
#0.0435
#For censor intervals with 1000 sims we get
#0.017

#TS.PAVSI (should approximately reject 0.05 times) 
sum(final.results$PAVSI.TS<=0.05)/N.sims
#For Ten.Fixed.Ints' approach with 10,000 sims we get
#0.059
#For censor intervals
#0.004



#TS.TFT (should approximately reject 0.05 times)
sum(final.results$TFT.TS<=0.05,rm.na=TRUE)/N.sims

#No issues with censor intervals variant
#Rejection at 0.004


#Still some NaN values
sum(final.results$TFT.TS=='NaN')/N.sims
#around 0.0106
sum(final.results$TFT.TS[final.results$TFT.TS!='NaN']<=0.05)/N.sims
#removing NaN values we get around 0.0291

#should also be approximately uniform [0,1]
hist(final.results$TFT.TS, breaks=20, probability = TRUE)
     
     