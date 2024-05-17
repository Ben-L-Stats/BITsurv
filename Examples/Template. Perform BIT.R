###################################################################################
#  Purpose: Perform the binomial interval test for a given model
#
#  Coder: Ben Lee
#  Date:  14-05-2024
###################################################################################


#libraries----------------------------------------------------------------------
library(survHE) #Load survHE. Loads flexsurv and survival automatically
library(survminer) #Required for ggsurvplot
library(dplyr)

#set file paths, load data and format data--------------------------------------

base.file<-"C:/R/Simulated KM curves/Final repo"

#melanoma data from https://github.com/SCFreeman/Melanoma_NMA
melanoma.data<-read.csv(file.path(base.file,'Data','melanoma_ipd.csv'))

melanoma.data <-melanoma.data %>% 
  mutate(arm=paste0(study, ' - TRT ',txCode ))

#unique(melanoma.data$arm)
uniq.arms<-"BREAK-3 - TRT 4"
j<-1

#select just 1 arm
surv.data<-melanoma.data %>% 
  filter(arm==uniq.arms[j]) %>%             #select arm and treatment
  mutate(USUBJID=as.character(patid)) %>%   #data formatting
  select(USUBJID,time, event)


#Perform the binomial intevral test (BIT)--------------------------------------------

#Specify intervals of interest--------------------------------------------------

#if not specifying intervals, then use the original approach. That is
# censors<-surv.data %>% filter(event==0)
# spec_int<-censors$time

#Recommended specifications are
censors<-surv.data %>% filter(event==0)
spec_int<-0.1*max(censors$time)*0:10


#Pick a survival model that you would like to check------------------------------

#exp        (1 parameter)
#weibull    (2 parameters)
#gompertz   .
#llogis     .
#lnorm      .
#gamma      .
#gengamma   (3 parameters)

Distribution="gengamma"


#Source and run BIT function--------------------------------------------------
#requires the sinib package

source(file.path(base.file, "Functions", "1. BITsurv.R"))
BIT.table<-BITsurv(surv.data, Distribution, spec_int)



#Source and run plot function---------------------------------------------------

source(file.path(base.file, "Functions", "2. BITplot.R"))
BIT.plot(surv.data, Distribution, BIT.table, X)



#Source and run test statistics-------------------------------------------------
source(file.path(base.file, "Functions", "3. Test statistics.R"))

#Obtain an approximate overall p-values
#That is, we assume the mid-p vals for each interval are distributed uniform[0,1]

#Returns a midpoint overall p-value
#Only uses the number of individual rejections to determine the overall p-value
BIT.TS.PAVSI(BIT.table$V.mid.pval)       #protection against very small intervals (PAVSI)

#Returns a continuous overall p-value
#Uses all of the information to determine the overall p-value
BIT.TS.TFT(BIT.table$V.mid.pval )       #transformed fisher test (TFT)



