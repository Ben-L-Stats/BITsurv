###################################################################################
#  Purpose: Perform the binomial interval test for exponential model for the 
#  BREAK-3 TRT 4 arm
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


#select just 1 arm
surv.data<-melanoma.data %>% 
  filter(arm=="BREAK-3 - TRT 4") %>%         #select arm and treatment
  mutate(USUBJID=as.character(patid)) %>%   #data formatting
  select(USUBJID,time, event)


#Perform the binomial intevral test (BIT)--------------------------------------------

#Specify intervals of interest--------------------------------------------------

#We are using the original approach where the censors define the intervals
censors<-surv.data %>% filter(event==0)
spec_int<-censors$time


#Pick a survival model that you would like to check------------------------------

#Options are: exp, weibull, gompertz, llogis, lnorm, gamma, gengamma   
#We are only interested in exp for this example
Distribution="exp"


#Source and run BIT function--------------------------------------------------

source(file.path(base.file, "Functions", "1. BITsurv.R"))
BIT.table<-BITsurv(surv.data, Distribution, spec_int)


#Source and run plot function---------------------------------------------------

source(file.path(base.file, "Functions", "2. BITplot.R"))
BIT.plot(surv.data, Distribution, BIT.table, X)
