################################################################################
#  Purpose: Perform the binomial interval test (BIT) for a given model
#
#  Programmer: 
#  Date:  
################################################################################

#libraries----------------------------------------------------------------------
library(dplyr)
library(BITsurv)

#create the survival data-------------------------------------------------------

simple.data<-data.frame(T=rexp(100, rate=1/10),        #the underlying event process
                        cens=runif(100, 0, 20)) %>%    #the underlying censor process
            mutate(event=ifelse(T<=cens,0, 1),
                   time=ifelse(T<=cens,T,cens)) %>%
            select(time,event)

#Specify intervals of interest--------------------------------------------------

#if not specifying intervals, then use the original approach. That is
# censors<-simple.data %>% filter(event==0)
# spec_int<-censors$time

#Recommended specifications are
censors<-simple.data %>% filter(event==0)
spec_int<-0.1*max(censors$time)*0:10


#Pick a survival model that you would like to check------------------------------

#Options are: exp, weibull, gompertz, llogis, lnorm, gamma, gengamma  
#All options are 2 parameter models, with the exception of exp (1 parameter) and
# gengamma (3 parameters)

Distribution="exp"

#An initial plot of the data and the fitted survival model----------------------

#?Fit.curve.plot
Fit.curve.plot(surv.data=simple.data,
               Distribution='exp',
               break.time=5)

#Produce BIT table--------------------------------------------------------------

#?BIT.surv
BIT.table<-BIT.surv(surv.data=simple.data,
                    Distribution=Distribution,
                    spec_int=spec_int)

#Run BIT plot function---------------------------------------------------------

#?BIT.plot
BIT.plot(surv.data=simple.data,
         BIT.table=BIT.table,
         break.time=5)

#Run test statistics------------------------------------------------------------
#Obtain an approximate overall p-values

#Uses the number of individual rejections to determine the overall p-value
#?BIT.TS.PAVSI
BIT.TS.PAVSI(BIT.table$V.mid.pval)       #protection against very small intervals (PAVSI)

#Uses all of the information to determine the overall p-value
#This TS is not recommended if using intervals determined by censoring as it becomes very conservative
#?BIT.TS.TFT
BIT.TS.TFT(BIT.table$V.mid.pval )       #transformed fisher test (TFT)
