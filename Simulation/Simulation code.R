###################################################################################
#  Purpose: Generate data under the null hypothesis.
#           Then apply the methods to this to verify type 1 error rates etc
#
#  Programmer: 
#  Date:  09/10/2025
###################################################################################

#Code structure:
# Lines 1-35   : descriptive details of code
# Lines 35-80  : specifying the simulation details - if you want to save
# the results you must add your own output.path
# Lines 80-140 : loading in libraries and required functions
# Lines 140+   : code for running the simulations and saving the results

#Overview
#Please note that this code is more extensive than the other examples.
#With the current specifications and seed, running this code should give 
#the same results as those presented in the simulation section of the paper

#The scenarios that this code explores are summarised by the scenarios dataframe.
#There are 12 scenarios and 4 different numbers of patients, resulting in 48 
#unique cases and thereby 480,000 simulated datasets.

#Runtime
#Running on a standard computer with 10,000 simulations takes around 10 hours. 

#Note: datasets simulated without censors
#For N=50 and Lambda=1/10 about 1 in 20,000 simulated datasets is generated
#without a censor, and therefore breaks the standard process.
#For each such case, a warning is flagged in the console and the dataset is
#rerun. No/minimal impact is expected on the results
#For this seed, this doesn't occur

################################################################################
#### Simulation specifications--------------------------------------------------
################################################################################

#You must specify your own file output directory
#Results are saved in github as excel files if you don't want to run the code
output.path<-"C:/.../BITsurv/Simulation/Results"

#This seed gives the results presented in the paper
#Use a different seed to see impact on the results
set.seed(314159)

#number of simulations
N.sims<-10000            #10,000 takes approximately 10 hours to run


#the different patient numbers considered
pat_vec<-c(50,100,200,500)    
#pat_vec<-100

#Create table of scenarios of interest
#Here I am interested in 12 scenarios
#2 interval selection approaches: 'ten' and 'censors'
#2 p-value approaches: rand and mid
#3 data maturity settings
#(For each of these scenarios the 4 different patients numbers will be applied)
scenarios<-data.frame(Approach=c(rep('Censor.Ints',6), 
                                 rep('Ten.Fixed.Ints',6)),
                      rate.spec=rep(c(1/10, 1/10, 
                                      1/30, 1/30, 
                                      1/70, 1/70), 2), 
                      p_method=rep(c('rand', 'mid'),6))


#say whether to save histograms of p-values
save.figs<-FALSE   #set to false by default for faster runtime
#setting this to true would result in PAVSI and TFT figures for each case 
#resulting in 48*2=96 saved images

#Creates a folder to save histograms of TFT and PAVSI p-values
if(save.figs){
  fig.path<-file.path(output.path, 'Figs')
  if (!dir.exists(fig.path) ) {
    dir.create(fig.path)}}


#libraries----------------------------------------------------------------------

library(BITsurv)  
#A lot of the machinery from BITsurv is explicitly written out as code here
#This is because we are constraining ourselves to test the null

library(dplyr)

library(flexsurv)
library(survminer)
library(PoissonBinomial)

library(openxlsx)


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





















################################################################################
################################################################################
#Running the simulation------------------------------------------------------
################################################################################
################################################################################


#loop through scenarios------------------------------------------------------
 for(k in 1:nrow(scenarios)){
  #k<-1
   
message(paste0('Table: ', scenarios$Approach[k], 
               ', Lambda=1/', 1/scenarios$rate.spec[k], 
               ', p_method=', scenarios$p_method[k]))


#expected survival in underlying model = 1/rate.spec
rate.spec<-scenarios$rate.spec[k]
   
#Interval selection approach
Approach<- scenarios$Approach[k]      #'Ten.Fixed.Ints' or 'Censor.Ints'
  
#p method
p_method<-scenarios$p_method[k]   #'mid for midpoint p-values, or 'rand' for randomised p-values
   






#loop through different numbers of patients-------------------------------------
for(j in 1:length(pat_vec)){     
#j<-1
  
start_time_pats <- Sys.time()

#number of patients
N.patients<-pat_vec[j]
#N.patients<-100




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



################################################################################
# Core simulation------------------------------------------------------------
################################################################################

#Loop through the number of simulations-----------------------------------------
for (i in 1:N.sims){
  #i<-1
  
  #Simulating survival data for the patients in this scenario
  surv.data<- data.frame(Time.of.Event=rexp(N.patients,rate=rate.spec),
                         C1=runif(N.patients, min=0, max=100),
                         C2=runif(N.patients, min=18, max=22)) %>% 
    mutate(C=pmin(C1,C2)) %>% 
    mutate(event=ifelse(Time.of.Event<C,1,0)) %>% 
    mutate(time=pmin(C,Time.of.Event))%>% 
    select(time, event) 
  
#-------------------------------------------------------------------------- 
# For N=50 and Lambda=1/10 about 1 in 20,000 simulated datasets is generated
# without a censor this reruns such datasets and flags them in the console
  
if( sum(surv.data$event==0)==0  ){  #if no censors
  #then rerun
  message("WARNING: A simulated dataset had 0 censors and was resimulated.")
  
  #extremely unlikely to be generated with 0 censors again
  surv.data<- data.frame(Time.of.Event=rexp(N.patients,rate=rate.spec),
                         C1=runif(N.patients, min=0, max=100),
                         C2=runif(N.patients, min=18, max=22)) %>% 
    mutate(C=pmin(C1,C2)) %>% 
    mutate(event=ifelse(Time.of.Event<C,1,0)) %>% 
    mutate(time=pmin(C,Time.of.Event))%>% 
    select(time, event) 
}
#--------------------------------------------------------------------------   
  
  
  #Decides which specified interval approach to use
  
  censors<-surv.data %>% filter(event==0)
  
  if (Approach=='Ten.Fixed.Ints'){
    spec_int<-0.1*max(censors$time)*0:10 
    int_method<-'ten'
  } 
  
  if (Approach=='Censor.Ints'){
    spec_int<-censors$time
    int_method<-'censors'
  }
  
  
  
  
  
  
  
  
  
  
  
  
  #Start of extracted BIT.surv function  
  #-------------------------------------------------------------------------------
  #The following is just the BIT.surv function for exp----------------------------
  #where a small edit has been made to use the true event generating process------
  #when calculating p. This edit is highlighted in the code------------------------
  #-------------------------------------------------------------------------------
  
  
  
  
  #Grouping and formatting the data-----------------------------------------------
  
  
  new.data<-rbind(surv.data %>% select(time,event),    #add the surv data
                  data.frame(time=spec_int,            #add information on the specified intervals
                             event=2)) %>%             #specified intervals are given a special event classification (event=2)
    arrange(time) %>%                       #for easier visualisation
    mutate(V.lower=lower.int.vec(time),     #place each datapoint within one of the specified intervals
           V.upper=upper.int.vec(time))
  
  
  
  
  #define the sub intervals I
  sub.ints<-new.data %>% filter(event!=1)
  sub.ints<-unique(sub.ints$time)       #subintervals are defined by all of the censor times and specified times
  
  
  
  
  
  
  #Now to assign each datapoint to it corresponding subinterval
  new.data<-new.data %>%
    mutate(I.lower=lower.sub.int.vec(time),     #assign datapoint to its subinterval
           I.upper=upper.sub.int.vec(time)) %>% 
    filter( !(I.lower==0 & I.upper==0))     #remove the first superfluos row
  
  #So far this is just a df that includes the appropriate grouping of V and I,
  #which can all be checked by reviewing the df
  
  
  
  
  #Obtaining statistics for each interval-----------------------------------------
  
  #We now want the number of observed events in each interval I
  new.data.1<-new.data %>%
    group_by(I.lower) %>%
    summarize(I.upper=unique(I.upper),
              V.lower=unique(V.lower),
              V.upper=unique(V.upper),
              Events.obs.I=sum(event==1),
              Censors.I=sum(event==0)) %>% 
    mutate(Pats.in.I=Events.obs.I+Censors.I)
  
  #Now we want to calculate the number at risk at the start of interval I

  #Edited**************************************************************************  
  #The proper code does not require the special loop accounting for only 1 censor
  
  if(nrow(new.data.1)==1){       #this catches weird errors when there is only 1 censor simulated
    new.data.1<-new.data.1%>% 
      mutate(N.risk=N.patients)
  }else{
    new.data.1<-new.data.1%>%
      mutate(N.risk=sum(new.data.1$Pats.in.I)-                      #total number of patients  -
               cumsum(c(0,                                      #patients lost in previous intervals
                        Pats.in.I[1:(nrow(new.data.1)-1)])) )
  }
  #End of edit**************************************************************************
  
  
  
  
  
  # Filtering based on considering valid intervals--------------------------------
  
  #First, remove intervals that end at infinity
  #An interval ending at infinity is just a caveat of the upper.int.vec function,
  # telling us that there are events after the last censor time
  #Due to our approach described in the paper, such events are not of interest
  new.data.1<-new.data.1%>%
    filter(!(V.upper==Inf))
  
  #Second, remove intervals for which N_{I_j}<1 as per equation 1
  #This is, remove the final interval if it is known that no events is possible
  #This prevents superfluous midpoint p-values of 0.5
  #Or in the case of randomised p-values, it prevents superfluous noise
  new.data.1<-new.data.1%>%
    filter( !((N.risk-Censors.I)==0)  )
  
  
  
  
  
  
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
  
  
  
  
  #Summarising by calculating p-values----------------------------------------------------------------
  
  #Now we look at summarizing across the specified intervals
  
  
  
  if(p_method=='mid'){
    
    new.data.3<-new.data.2 %>%
      mutate(expected_E=p*(N.risk-Censors.I) )%>%    #first add the expect number of events under the fitted model
      group_by(V.lower) %>%
      summarise(V.upper=unique(V.upper),
                Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
                Observed.E_over.V=sum(Events.obs.I),
                
                V.pval= {0.5*PoissonBinomial::ppbinom(
                  x=as.integer(sum(Events.obs.I)),
                  probs=(rep(p, (N.risk-Censors.I) )))+
                    0.5*PoissonBinomial::ppbinom(
                      x=as.integer(sum(Events.obs.I)-1),
                      probs=(rep(p,(N.risk-Censors.I) )))},
                
                N.risk.at.V.start=max(N.risk),
                E=sum(Events.obs.I))
  } 
  
  
  
  
  if(p_method=='rand'){
    #The only difference in the code compared to p_method=='mid' is how V.pval is defined
    
    #Specifically, we now use randomised p-values 
    #As per Rubin-Delanchy 2019 (Meta-Analysis of Mid-p-Values), randomised Ps can be defined as
    # R = X*P1 + (1-X)*P2      (aside: the direction of the Ps is flipped compared to Delanchy)
    # where X is a uniform 0, 1 r.v
    #rearranging we have: R= X(P1-P2)+P2
    #This is the form used for obtaining randomised p-values when specifying V.pval in our code,
    #and this form is used to avoid storing X values or other code complications 
    
    new.data.3<-new.data.2 %>%
      mutate(expected_E=p*(N.risk-Censors.I) )%>%    #first add the expect number of events under the fitted model
      group_by(V.lower) %>%
      summarise(V.upper=unique(V.upper),
                Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
                Observed.E_over.V=sum(Events.obs.I),
                
                V.pval= { runif(1)*
                    (PoissonBinomial::ppbinom(x=as.integer(sum(Events.obs.I)),
                                              probs=(rep(p, (N.risk-Censors.I) )))  -
                       PoissonBinomial::ppbinom(x=as.integer(sum(Events.obs.I)-1),
                                                probs=(rep(p, (N.risk-Censors.I) )))) +
                    PoissonBinomial::ppbinom(x=as.integer(sum(Events.obs.I)-1),
                                             probs=(rep(p, (N.risk-Censors.I) )))  },
                
                N.risk.at.V.start=max(N.risk),
                E=sum(Events.obs.I))
    
  }
  
  
 
  
  
  #Bonferroni and individual test results-----------------------------------------
  #for the specified intervals
  
  information<-new.data.3 %>%
    mutate(individual.test=ifelse(V.pval<=0.025| V.pval>=0.975,
                                  'Reject', 'Accept')) %>%
    mutate(bonferroni.test=ifelse(V.pval<=         
                                    0.025/nrow(new.data.3) |
                                    V.pval>=
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
                          PAVSI.TS=BIT.TS.PAVSI(information$V.pval,p_method=p_method,print=FALSE),      
                          TFT.TS=BIT.TS.TFT(information$V.pval,  print=FALSE),
                          Intervals=nrow(information)) 
  
  #save as final results
  if(i==1){
    final.results<-sim.results 
  }else{
    final.results<-rbind(final.results,
                         sim.results)
  }
  
  
  
} #end loop of simulations------------------------------------------------------




#Save p-values histograms
if(save.figs){

#Save TFT histograms
png(file.path(fig.path,
              paste0('TFT_',
                     scenarios$Approach[k],
                     '_',
                     'Recip.Lambda.',
                      1/scenarios$rate.spec[k],
                     '_',
                     p_method,
                     '_pats',
                     N.patients,
                     '.png')), 
              width = 800, height = 600)  
  
hist(final.results$TFT.TS, breaks=20)

dev.off()



#Save PAVSI histograms
png(file.path(fig.path,
              paste0('PAVSI_',
                     scenarios$Approach[k],
                     '_',
                     'Recip.Lambda.',
                     1/scenarios$rate.spec[k],
                     '_',
                     p_method,
                     '_pats',
                     N.patients,
                     '.png')), 
    width = 800, height = 600)  

hist(final.results$PAVSI.TS, breaks=20)

dev.off()
  

} #end of if loop to run figures

















################################################################################
#Saving table of simulation scenarios results for different N--------------------------------
################################################################################

runtime <- Sys.time() - start_time_pats

#formatting
runtime<-paste(round(as.numeric(runtime),1),
      attr(runtime, "units"))
      


sim.sum.results<-data.frame(N=N.patients,
                            Bonferroni.rate=sum(final.results$Bonferroni.rejects>0)/N.sims,
                            TFT.rate=sum(final.results$TFT.TS<=0.05,rm.na=TRUE)/N.sims,
                            PAVSI.rate=sum(final.results$PAVSI.TS<=0.05)/N.sims,   
                            Lambda=rate.spec,
                            Approach=Approach,
                            p_method=p_method,
                            N.sims=N.sims,
                            Runtime=runtime) 



#save as final results
if(j==1){
  output.results<-sim.sum.results 
}else{
  output.results<-rbind(output.results,
                        sim.sum.results)
}





message(paste0('  -', N.sims, ' simulations run for ', 
               N.patients, ' patients scenario. Runtime: ',
               runtime))


}
#end of patient loops







#save results for patient loop-------------------------------------------
file.name<-paste0( 'p.', p_method,
                   '_', Approach, 
                    '_mean.surv.', 1/rate.spec,
                   '.xlsx')

#write dataframe to excel file
write.xlsx(
  list(Sheet1 = output.results),
  file = file.path(output.path, file.name)
)





 }
#end of loop through scenarios
