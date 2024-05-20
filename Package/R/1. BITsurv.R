#Fit exponential model to single arm data and produce plot of this

#Input is:
# surv.data
# Distribution
# time.breaks

# surv.data with numeric time and event values

#and requires libraries
# library(survHE) #Load survHE. Loads flexsurv and survival automatically
# library(survminer) #Required for ggsurvplot
# library(dplyr)
# library(sinib)





#' Fits a parametric survival curve to data and returns a binomial interval test (BIT) summary
#'
#' This function produces a dataframe.
#'
#' @param surv.data A dataframe with 2 columns labelled time and event.
#' @param Distribution Pick your parametric survival distribution. Available options are:
#' 'exp', 'weibull', 'gompertz', 'llogis','lnorm', 'gamma', and 'gengamma'. See
#' flexsurv for further details on these distributions.
#' @param spec_int The specified intervals: Provide a vector of values that represent your chosen intervals,
#' where the sequential values in this vector form an interval.
#'
#' @return A binomial interval test (BIT) summary in the form of a dataframe
#' @export
#'
#' @examples
#' #The dplyr package is used to make the simple.data
#' simple.data<-data.frame(T=rexp(100, rate=1/10),        #the underlying event process
#'                         cens=runif(100, 0, 20)) %>%    #the underlying censor process
#'      mutate(event=ifelse(T<=cens,0, 1),
#'             time=ifelse(T<=cens,T,cens)) %>%
#'      select(time,event)
#'
#' #Here we use the 10 evenly spaced interval approach:
#'   censors<-simple.data %>% filter(event==0)
#'   spec_int<-0.1*max(censors$time)*0:10
#'
#' #Run BITsurv
#' BIT.surv(surv.data=simple.data,
#'         Distribution='exp',
#'         spec_int=spec_int)
#'
BIT.surv<-function(surv.data, Distribution, spec_int){
#start of function------------------------------------------------------------
#-----------------------------------------------------------------------------

#Fit the  model to the data-----------------------------------------------------

par.est<-flexsurv::flexsurvreg(survival::Surv(time, event) ~ 1,
                               data = surv.data,
                               dist=Distribution)


#Grouping and formatting the data-----------------------------------------------

  lower.int<-function(time){ifelse(time<=min(spec_int),
                                   0,           #forces the lowest interval to have a lower bound of 0
                                   max(spec_int[spec_int<time]))}
  lower.int.vec<-function(time){sapply(time, lower.int)}

  upper.int<-function(time){ifelse(time>max(spec_int),Inf,   #forces the highest interval to have an upper bound of infinity
                                   min(spec_int[spec_int>=time]))}
  upper.int.vec<-function(time){sapply(time, upper.int)}


#start to create the dataframe describing the intervals
new.data<-rbind(surv.data %>% select(time,event),    #add the surv data
                data.frame(time=spec_int,            #add information on the specified intervals
                           event=2)) %>%             #specified intervals are given a special event classification (event=2)
  arrange(time) %>%                       #for easier visualisation
  mutate(V.lower=lower.int.vec(time),     #place each datapoint within one of the specified intervals
         V.upper=upper.int.vec(time))

#define the sub intervals I
sub.ints<-new.data %>% filter(event!=1)
sub.ints<-unique(sub.ints$time)       #subintervals are defined by all of the censor times and specified times

lower.sub.int<-function(time){ifelse(time<=min(sub.ints),
                                     0,   #forces the lowest subinterval to have a lower bound of 0
                                     max(sub.ints[sub.ints<time]))}
lower.sub.int.vec<-function(time){sapply(time, lower.sub.int)}


upper.sub.int<-function(time){ifelse(time>max(sub.ints),
                                     Inf, #forces the highest subinterval to have an upper bound of infinity
                                     min(sub.ints[sub.ints>=time]))}
upper.sub.int.vec<-function(time){sapply(time, upper.sub.int)}

#Now to assign each datapoint to it corresponding subinterval
new.data<-new.data %>%
  filter(event!=2) %>%    #remove the specified intervals as these are no longer of interest
  mutate(I.lower=lower.sub.int.vec(time),     #assign datapoint to its subinterval
         I.upper=upper.sub.int.vec(time))

#So far this is just a df that includes the appropriate grouping of V and I,
#which can all be checked by reviewing the df







#Obtaining statistics for each interval-----------------------------------------

#We now want the number of observed events in each interval I
new.data.1<-new.data %>%
  group_by(I.lower) %>%
  summarize(I.upper=unique(I.upper),
            V.lower=unique(V.lower),
            V.upper=unique(V.upper),
            Events.obs.I=sum(event),
            Pats.in.I=length(event))

#Now we want to calculate the number at risk at the start of interval I
new.data.1<-new.data.1%>%
  mutate(N.risk=sum(new.data.1$Pats.in.I)-                      #total number of patients  -
           cumsum(c(0,                                      #patients lost in previous intervals
                    Pats.in.I[1:(nrow(new.data.1)-1)])) )



# Calculating p_{I_j} values----------------------------------------------------

#We also want to calculate the value p for interval I
#where p=P(patient experiences an event in I|they make it to I)
#We define p based on the fitted model that is being used

if (Distribution=='exp'){
new.data.2<-new.data.1%>%
  mutate(p={pexp(q=I.upper, rate=exp(par.est$coefficients))-       #calculate of p for fitted exponential model
            pexp(q=I.lower, rate=exp(par.est$coefficients)) }/
           (1-pexp(q=I.lower,rate=exp(par.est$coefficients))))
                         }

if (Distribution=='weibull'){
  new.data.2<-new.data.1%>%
    mutate(p={pweibull(q=I.upper,
                       shape=exp(par.est$coefficients['shape']),
                       scale=exp(par.est$coefficients['scale']))-       #calculate of p
              pweibull(q=I.lower,
                       shape=exp(par.est$coefficients['shape']),
                       scale=exp(par.est$coefficients['scale'])) }/
          (1-pweibull(q=I.lower,
                      shape=exp(par.est$coefficients['shape']),
                      scale=exp(par.est$coefficients['scale'])))   )
                            }

if (Distribution=='gompertz'){
  new.data.2<-new.data.1%>%
    mutate(p={flexsurv::pgompertz(q=I.upper,
                       shape=par.est$coefficients['shape'],
                       rate=exp(par.est$coefficients['rate']))-       #calculate of p
        flexsurv::pgompertz(q=I.lower,
                       shape=par.est$coefficients['shape'],
                       rate=exp(par.est$coefficients['rate'])) }/
          (1-flexsurv::pgompertz(q=I.lower,
                       shape=par.est$coefficients['shape'],
                       rate=exp(par.est$coefficients['rate'])))   )
                             }


if (Distribution=='llogis'){
  new.data.2<-new.data.1%>%
    mutate(p={flexsurv::pllogis(q=I.upper,
                       shape=exp(par.est$coefficients['shape']),
                       scale=exp(par.est$coefficients['scale']))-       #calculate of p
        flexsurv::pllogis(q=I.lower,
                 shape=exp(par.est$coefficients['shape']),
                 scale=exp(par.est$coefficients['scale'])) }/
          (1-flexsurv::pllogis(q=I.lower,
                      shape=exp(par.est$coefficients['shape']),
                      scale=exp(par.est$coefficients['scale'])))   )
                           }


if (Distribution=='lnorm'){
  new.data.2<-new.data.1%>%
    mutate(p={plnorm(q=I.upper,
                     meanlog=par.est$coefficients['meanlog'],
                     sdlog=exp(par.est$coefficients['sdlog']))-       #calculate of p
        plnorm(q=I.lower,
               meanlog=par.est$coefficients['meanlog'],
               sdlog=exp(par.est$coefficients['sdlog'])) }/
          (1-plnorm(q=I.lower,
                    meanlog=par.est$coefficients['meanlog'],
                    sdlog=exp(par.est$coefficients['sdlog'])))   )
}

if (Distribution=='gamma'){
  new.data.2<-new.data.1%>%
    mutate(p={pgamma(q=I.upper,
                      shape=exp(par.est$coefficients['shape']),
                      rate=exp(par.est$coefficients['rate']))-       #calculate of p
        pgamma(q=I.lower,
                shape=exp(par.est$coefficients['shape']),
                rate=exp(par.est$coefficients['rate'])) }/
          (1-pgamma(q=I.lower,
                     shape=exp(par.est$coefficients['shape']),
                     rate=exp(par.est$coefficients['rate'])))   )
}


if (Distribution=='gengamma'){
  new.data.2<-new.data.1%>%
    mutate(p={flexsurv::pgengamma(q=I.upper,
                     mu=par.est$coefficients['mu'],
                     sigma=exp(par.est$coefficients['sigma']),
                     Q=par.est$coefficients['Q'])-       #calculate of p
        flexsurv::pgengamma(q=I.lower,
               mu=par.est$coefficients['mu'],
               sigma=exp(par.est$coefficients['sigma']),
               Q=par.est$coefficients['Q']) }/
          (1-flexsurv::pgengamma(q=I.lower,
                    mu=par.est$coefficients['mu'],
                    sigma=exp(par.est$coefficients['sigma']),
                    Q=par.est$coefficients['Q']))   )
}








#Final summaries----------------------------------------------------------------

#Now we look at summarizing the specified intervals

new.data.3<-new.data.2 %>%
  mutate(expected_E=p*N.risk)%>%    #first add the expect number of events under the fitted model
  group_by(V.lower) %>%
  summarise(V.upper=unique(V.upper),
            Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
            Observed.E_over.V=sum(Events.obs.I),
            V.mid.pval=if(length(V.lower)==1){   #there are cases where sinib breaks in the binomial case
              #to prevents these breaks (p.values outside (0,1)) we use the classical pbinom function

              0.5*pbinom(q=as.integer(sum(Events.obs.I)),    #p-value calculation for each V
                         size=as.integer(N.risk),
                         prob=p)+
                0.5*pbinom(q=as.integer(sum(Events.obs.I)-1),
                           size=as.integer(N.risk),
                           prob=p)

            }else{ #for non-binomial situations sinib is used

              0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)),    #p-value calculation for each V
                         size=as.integer(N.risk),
                         prob=p) +
                ifelse(sum(Events.obs.I)==0,       #additional step as psinib can break if q<0
                       0,
                       0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)-1),
                                  size=as.integer(N.risk),
                                  prob=p))
            },
            N.risk.at.V.start=max(N.risk),
            E=sum(Events.obs.I))

#In the case that the largest time is an event and the specified V are such that
#the last event is outside of V then you obtain an interval with an upper bound of
#resulting in a Prob=1 of observing an event within said interval, making the observed
#event redundant. As such, such an interval is not of interest and is removed.

#Similarly, if for whatever reason V is specified such that certain events are
#outside of this. Then they will be assigned to this infinite interval and similarly
#should be dropped. This is done here
new.data.3<-new.data.3 %>% filter(V.upper!=Inf)



#Error catching-------------------------------------------------------------------

#There are a few rare situations (~1%) where sinib does not converge
#and gives p-values outside of [0,1].
#We will now implement a step to protect against this

#First check new.data.3 for any p values outside of [0,1]

p.checks<-new.data.3$V.mid.pval
if(sum(!between(p.checks, left=0, right=1))>=1){ #if one or more p-values fall
  #outside of 0 and 1 then
  #find the v.lower values corresponding to these errors
  error.intervals<-new.data.3$V.lower[!between(p.checks, left=0, right=1)]

  #now we want to obtain the corrected values
  #we do this be approximating the probabilities by using the round function
  #this should hopefully provide values that have converged
  ammended.p<-new.data.2 %>%
    filter(V.lower %in% error.intervals) %>%   #select only data related to the error intervals
    mutate(expected_E=p*N.risk)%>%
    group_by(V.lower) %>%
    summarise(V.upper=unique(V.upper),
              Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
              Observed.E_over.V=sum(Events.obs.I),
              V.mid.pval= 0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)),    #p-value calculation for each V
                                     size=as.integer(N.risk),
                                     prob=round(p, digits = 3)) +##################
              ifelse(sum(Events.obs.I)==0,       #additional step as psinib can break if q<0
                     0,
                     0.5*sinib::psinib(q=as.integer(sum(Events.obs.I)-1),
                                size=as.integer(N.risk),
                                prob=round(p, digits = 3))),############
              N.risk.at.V.start=max(N.risk),
              E=sum(Events.obs.I) )


  #Now to replace the error p.val with this updated one

  #values outside of range are given by new.data.3$V.mid.pval[which(!between(p.checks, left=0, right=1))]
  #and are updated to ammendment
  new.data.3$V.mid.pval[which(!between(p.checks, left=0, right=1))]<-ammended.p$V.mid.pval
}
#End of error catching





#Bonferroni and individual test results-----------------------------------------
#for the specified intervals

information<-new.data.3 %>%
  mutate(individual.test=ifelse(V.mid.pval<=0.025| V.mid.pval>=0.975,
                                'Reject', 'Accept')) %>%
  mutate(bonferroni.test=ifelse(V.mid.pval<=         #I need to derive this 2-tailed result as it is not in my lecture notes
                                  0.025/nrow(new.data.3) |
                                V.mid.pval>=
                                  1- 0.025/nrow(new.data.3),
                                'Reject',
                                'Accept')) %>%
  select(-E) %>%
  mutate(fitted.dist=Distribution)


return(information)

#end of function------------------------------------------------------------
#-----------------------------------------------------------------------------
}

