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


#Pick a survival model that you would like to check------------------------------

#Options are: exp, weibull, gompertz, llogis, lnorm, gamma, gengamma  
#All options are 2 parameter models, with the exception of exp (1 parameter) and
# gengamma (3 parameters)
Distribution="exp"

#An initial plot of the data and the fitted survival model----------------------

#?Fit.curve.plot
Fit.curve.plot(surv.data=simple.data,
               Distribution=Distribution,
               break.time=5)

#Produce BIT table--------------------------------------------------------------

#Decide on the interval selection approach
#Either 'ten' (generally recommended) or 'censors'
#?BIT.surv
BIT.table<-BIT.surv(surv.data=simple.data,
                    Distribution=Distribution,
                    int_method='ten')
#result
BIT.table


#Run BIT plot function---------------------------------------------------------

#?BIT.plot
BIT.plot(surv.data=simple.data,
         BIT.table=BIT.table,
         break.time=5)


#Run test statistics------------------------------------------------------------
#Obtain an approximate overall p-values
#See simulation results in the paper to understand the limitations 

#Uses the number of individual rejections to determine the overall p-value
#?BIT.TS.PAVSI
BIT.TS.PAVSI(BIT.table$V.pval)       #protection against very small intervals (PAVSI)

#Uses all of the information to determine the overall p-value
#?BIT.TS.TFT
BIT.TS.TFT(BIT.table$V.pval )       #transformed fisher test (TFT)
