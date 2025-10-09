###################################################################################
#  Purpose: Perform the binomial interval test for exponential model for the 
#  BREAK-3 TRT 4 arm
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
  filter(arm=="BREAK-3 - TRT 4") %>%         #select arm and treatment
  mutate(USUBJID=as.character(patid)) %>%   #data formatting
  select(USUBJID,time, event)

#Perform the binomial interval test (BIT)--------------------------------------------

#Pick a survival model that you would like to check------------------------------

#Options are: exp, weibull, gompertz, llogis, lnorm, gamma, gengamma   
#We are only interested in exp for this example
Distribution="exp"

#Source and run BIT function--------------------------------------------------

#For example 1 we are interested in the censor approach
BIT.table<-BIT.surv(surv.data=surv.data,
                    Distribution=Distribution,
                    int_method='censors')


#Source and run plot function---------------------------------------------------

BIT.plot(surv.data, BIT.table)

#The test statistics------------------------------------------------------------

BIT.TS.TFT(BIT.table$V.pval )       #transformed fisher test (TFT)

BIT.TS.PAVSI(BIT.table$V.pval)       #protection against very small intervals (PAVSI)
