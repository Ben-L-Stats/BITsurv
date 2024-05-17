library(devtools)
library(dplyr)

simple.data<-data.frame(T=rexp(100, rate=1/10),
                        cens=runif(100, 0, 20)) %>%
  mutate(event=ifelse(T<=cens,
                      0, 1),
         time=ifelse(T<=cens,
                     T,cens))

load_all()

#--------------------------------------------------------------------------------
?Fit.curve.plot
Fit.curve.plot(surv.data=simple.data,
               Distribution='exp',
               break.time = 5)

#--------------------------------------------------------------------------------
?BIT.surv

#Here we use the 10 evenly spaced interval approach:
censors<-simple.data %>% filter(event==0)
spec_int<-0.1*max(censors$time)*0:10

#Run BITsurv
BIT.table<-BIT.surv(surv.data=simple.data,
        Distribution='exp',
        spec_int=spec_int)

#--------------------------------------------------------------------------------
?BIT.plot

BIT.plot(surv.data=simple.data,
         Distribution='exp',
         BIT.table=BIT.table,
         break.time=5)


#--------------------------------------------------------------------------------
?BIT.TS.PAVSI
BIT.TS.PAVSI(p.vals=BIT.table$V.mid.pval)
BIT.TS.PAVSI(p.vals=0.1*c(1:9))

?BIT.TS.TFT
BIT.TS.TFT(p.vals=BIT.table$V.mid.pval)
BIT.TS.TFT(p.vals=0.1*c(1:9))

