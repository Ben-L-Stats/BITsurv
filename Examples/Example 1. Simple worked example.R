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

#Specify intervals of interest--------------------------------------------------

#We are using the original approach where the censors define the intervals
censors<-surv.data %>% filter(event==0)
spec_int<-censors$time

#Pick a survival model that you would like to check------------------------------

#Options are: exp, weibull, gompertz, llogis, lnorm, gamma, gengamma   
#We are only interested in exp for this example
Distribution="exp"

#Source and run BIT function--------------------------------------------------

BIT.table<-BIT.surv(surv.data, Distribution, spec_int)


#Source and run plot function---------------------------------------------------

BIT.plot(surv.data, BIT.table)

#The test statistics------------------------------------------------------------

BIT.TS.PAVSI(BIT.table$V.mid.pval)       #protection against very small intervals (PAVSI)

BIT.TS.TFT(BIT.table$V.mid.pval )       #transformed fisher test (TFT)
