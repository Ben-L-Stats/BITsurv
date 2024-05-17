###################################################################################
#  Purpose: Generate data under the null hypothesis.
#           Then plot this to visualize a single simulation
#
#  Coder: Ben Lee
#  Date:  14-05-2024
###################################################################################


#libraries----------------------------------------------------------------------
library(survHE) #Load survHE. Loads flexsurv and survival automatically
library(survminer) #Required for ggsurvplot
library(dplyr)

base.file<-"C:/R/Simulated KM curves/Final repo"

#Produce simulated data---------------------------------------------------------

#Event generating process is
# T~exp(lambda=1/10)

#Censor generating process is
# C1~uniform[0,100]      #general censoring
# C2~uniform[18,22]      #heavy censoring at the end of the curve
# C<-min(C1,C2)          #final censor value    

#100 patients

#Example simulation
N.patients<-100

simulation<- data.frame(T=rexp(N.patients,rate=1/10),
                       C1=runif(N.patients, min=0, max=100),
                       C2=runif(N.patients, min=18, max=22)) %>% 
  mutate(C=pmin(C1,C2)) %>% 
  mutate(event=ifelse(T<C,1,0)) %>% 
  mutate(time=pmin(C,T))

#Format as survival data
surv.data<-simulation %>% 
  select(time, event) 



#Produce plot-------------------------------------------------------------------

Distribution="exp"

#Fit the model to the data
par.est<-flexsurvreg(Surv(time, event) ~ 1, data = surv.data, 
                     dist=Distribution) 

#define the survival data y, for the fitted curve
#x is specified in as 
#the resolution and boundaries of the fitted curve
x<-round(max(surv.data$time))*(1:1000)/1000   #1000 times up to the last observed time
#number at risk itervals are defaulted to every 10 time units

#Exp curve using MLE for lambda from simulated data
y.MLE<-1-pexp(q=x,rate=exp(par.est$coefficients))

#Exp curve using ture lambda value
y.true<-1-pexp(q=x,rate=1/10)


#Plotting-------------------------------------------------------------------------------
#get data in appropriate form for initial plot
KM.est <- survfit(Surv(time,event) ~ 1, 
                  data = surv.data)
#initial plot
plot<-ggsurvplot(KM.est, data = surv.data,
                 title= paste0('Fitted Survival curve - ', Distribution),
                 risk.table = TRUE,
                 palette = c("purple"),
                 break.time.by = 5)

#Plot with exp curves over the top-----------------------------------

#KM curve with fitted curve on top
plot$plot+
   geom_line(data=data.frame(x=x,y=y.MLE), 
                           aes(x,y.MLE),
                           size=1,
                           color='black')+
  geom_line(data=data.frame(x=x,y=y.true), 
                           aes(x,y.true),
                           size=1,
                           color='red')

#The red curve shows the true event generating process
#The black curves shows the estimate of the event generating process using MLE and the exp model

#When testing the null by simulation, we will obtain all of our estimates using only
#the simulated data and the red curve (the black curve is not of interest)

