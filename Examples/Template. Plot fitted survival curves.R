###################################################################################
#  Purpose: Fit and plot survival curves
#
#  Coder: Ben Lee
#  Date:  13-05-2024
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


#Fit model to the data---------------------------------------------------------

#Pick a survival model

#exp        (1 parameter)
#weibull    (2 parameters)
#gompertz   .
#llogis     .
#lnorm      .
#gamma      .
#gengamma   (3 parameters)


Distribution="gengamma"


#Source and run plot function
source(file.path(base.file, "Functions", "0. Fit curves and plot.R"))
#defaults to number at risk interval of 10 time units
Fit.curve.plot(surv.data, Distribution)

