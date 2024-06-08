###################################################################################
#  Purpose: Perform the binomial interval tests for the seven standard parametric
#  survival distributions
#
#  Programmer: 
#  Date:  
###################################################################################

#libraries----------------------------------------------------------------------
library(BITsurv)
library(dplyr)

#set file paths, load data and format data--------------------------------------

#melanoma data from https://github.com/SCFreeman/Melanoma_NMA 
#this is saved within the package
melanoma.data

#select just 1 arm from this data
surv.data<-melanoma.data %>% 
  filter(arm=="COMBI-d - TRT 5") %>%         #select arm and treatment
  mutate(USUBJID=as.character(patid)) %>%    #data formatting
  select(USUBJID,time, event)


#Quick initial plot to view the data--------------------------------------------
Fit.curve.plot(surv.data=surv.data,
               Distribution='exp',
               break.time=5)

#Perform the binomial intevral test (BIT)--------------------------------------------

#Specify intervals of interest--------------------------------------------------

#Decide which approach to use:
#here we go with the 10 evenly-spaced interval approach
#replace lines 44 and 45 with uncommented lines 40 and 41 to use the censor interval approach 

#The censor interval approach
# censors<-surv.data %>% filter(event==0)
# spec_int<-censors$time

#The 10 evenly-spaced intervals approach
censors<-surv.data %>% filter(event==0)
spec_int<-0.1*max(censors$time)*0:10


#Pick a survival model that you would like to check------------------------------

#Options are: exp, weibull, gompertz, llogis, lnorm, gamma, gengamma   
#We are interested in all of these models, so we will loop through them
#and save the results
Dist.loop=c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma", "gengamma")


#Run analyses using original censor intervals approach--------------------------

#set up empty dataframe to save results
TS.results<-data.frame(fitted.dist=rep(NA,length(Dist.loop)),
                       PAVSI=rep(NA,length(Dist.loop)),
                       TFT=rep(NA,length(Dist.loop)))


#Produce the BIT tables for each distribution-----------------------------------

for (i in 1:length(Dist.loop)){ #loop through distributions

#run the BIT
BIT.table<-BIT.surv(surv.data, Distribution=Dist.loop[i], spec_int)

#Save the TS results
TS.results$fitted.dist[i]<- Dist.loop[i]
TS.results$PAVSI[i]<- BIT.TS.PAVSI(BIT.table$V.mid.pval, print=FALSE) 
TS.results$TFT[i]<- BIT.TS.TFT(BIT.table$V.mid.pval, print=FALSE)

#Save the BIT results
if (i==1){ #save first result
 BIT.tab.final<-BIT.table  
  } else{    #bind results together
BIT.tab.final<-rbind(BIT.tab.final,
                       BIT.table)
}
  
} #end of loop through distributions


#Summarize the results----------------------------------------------------------
Results.Ints<-BIT.tab.final %>% 
  group_by(fitted.dist) %>% 
  summarize(Intervals=length(fitted.dist),
            Ints.Bonferroni=sum(bonferroni.test=='Reject'),
            Ints.Individual=sum(individual.test=='Reject')) %>% 
  left_join(TS.results,
            by='fitted.dist')



#Collect the AIC and BIC values for each distribution---------------------------

for (i in 1:length(Dist.loop)){ #loop through distributions
  
  #fit the model to the data
  par.est<-flexsurv::flexsurvreg(survival::Surv(time, event) ~ 1, 
                       data = surv.data, dist=Dist.loop[i]) 
  
  #save the AIC values
  if (i==1){ #save first result
    AIC.table<-data.frame(fitted.dist=Dist.loop[i], AIC=AIC(par.est), BIC=BIC(par.est))
  } else{    #bind results together
    AIC.table<-rbind(AIC.table,
                     data.frame(fitted.dist=Dist.loop[i], AIC=AIC(par.est), BIC=BIC(par.est)))
  }
  
}#end of loop through distributions


#Bind the tables together-------------------------------------------------------

Final.table<-Results.Ints %>% 
               left_join(AIC.table, by='fitted.dist')









#Additional individual plot checks as required---------------------------------------------

#lnorm perform well on AIC and BIC but has a Bonferroni rejection for censor intervals

spec.dist<-'gengamma'

BIT.spec.table<-BIT.surv(surv.data, 
                   Distribution=spec.dist, 
                   spec_int=0.1*max(censors$time)*0:10)

# BIT.spec.table<-BITsurv(surv.data, 
#                         Distribution=spec.dist, 
#                         spec_int=censors$time)

BIT.plot(surv.data, 
         BIT.table=BIT.spec.table, 
         break.time=5)

#-------------------------------------------------------------------------------
spec.dist<-'lnorm'

BIT.spec.table<-BIT.surv(surv.data, 
                        Distribution=spec.dist, 
                        spec_int=0.1*max(censors$time)*0:10)

BIT.plot(surv.data, 
         BIT.table=BIT.spec.table, 
         break.time=5)
