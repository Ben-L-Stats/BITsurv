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
# library(PoissonBinomial)





#' Fits a parametric survival curve to data and returns a binomial interval test (BIT) summary
#'
#' This function produces a dataframe.
#'
#' @param surv.data A dataframe with 2 columns labelled time and event.
#' @param Distribution Pick your parametric survival distribution. Available options are:
#' 'exp', 'weibull', 'gompertz', 'llogis','lnorm', 'gamma', and 'gengamma'. See
#' flexsurv for further details on these distributions.
#' @param int_method The interval method. Options are 'ten' for 10 evenly-spaced
#' intervals up to the last censor and 'censors'. To specify custom intervals 
#' please edit the BIT.surv function directly.
#' @param p_method The method of obtaining p-values: Default is 'mid', which is recommended for
#' practical applications. Options are either 'mid' for midpoint p-values,
#' or 'rand' for randomised p-values. 
#'
#' @return A binomial interval test (BIT) summary in the form of a dataframe
#' @export
#'
#' @examples
#' #The dplyr package is required for BIT.surv to work
#' library(dplyr)
#' 
#' simple.data<-data.frame(T=rexp(100, rate=1/10),        #the underlying event process
#'                         cens=runif(100, 0, 20)) %>%    #the underlying censor process
#'      mutate(event=ifelse(T<=cens,0, 1),
#'             time=ifelse(T<=cens,T,cens)) %>%
#'      select(time,event)
#'
#' #Here we use the 10 evenly-spaced interval approach for the exponential model
#' BIT.surv(surv.data=simple.data,
#'         Distribution='exp',
#'         int_method='ten')
#'






BIT.surv<-function(surv.data, 
                   Distribution, 
                   int_method,
                   p_method='mid'){
#start of function------------------------------------------------------------
#-----------------------------------------------------------------------------

  # Catch errors in specifications
  if (!(Distribution %in%  c('exp', 'weibull', 'gompertz', 'llogis','lnorm', 'gamma', 'gengamma'))) {
    stop("The distribution you have specified is not an available option.")}
  
  if (!(p_method %in%  c('mid', 'rand'))) {
    stop("The p_method you have specified is not an available option.")}
  
  if (!(int_method %in%  c('ten', 'censors'))) {
    stop("The int_method you have specified is not an available option.")}
  

#Define the specified intervals-----------------------------------------------------

if (int_method =='censors'){
  censors<-surv.data %>% filter(event==0)
  spec_int<-censors$time
}

if (int_method =='ten'){
  censors<-surv.data %>% filter(event==0)
  spec_int<-0.1*max(censors$time)*0:10
}
  
  
    
  
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
  mutate(I.lower=lower.sub.int.vec(time),     #assign datapoint to its subinterval
         I.upper=upper.sub.int.vec(time)) %>% 
  filter( !(I.lower==0 & I.upper==0))     #remove the first superfluos row

#So far this is just a df that includes the appropriate grouping of V and I,
#which can all be checked by reviewing the df







#Obtaining statistics for each interval-----------------------------------------

#We now want the number of observed events in each interval I
new.data.1<-new.data %>%
  group_by(I.lower) %>%
  summarize(I.upper=unique(I.upper),
            V.lower=unique(V.lower),
            V.upper=unique(V.upper),
            Events.obs.I=sum(event==1),
            Censors.I=sum(event==0)) %>% 
  mutate(Pats.in.I=Events.obs.I+Censors.I)

#Now we want to calculate the number at risk at the start of interval I
new.data.1<-new.data.1%>%
  mutate(N.risk=sum(new.data.1$Pats.in.I)-                      #total number of patients  -
           cumsum(c(0,                                      #patients lost in previous intervals
                    Pats.in.I[1:(nrow(new.data.1)-1)])) )


# Filtering based on considering valid intervals--------------------------------

#First, remove intervals that end at infinity
#An interval ending at infinity is just a caveat of the upper.int.vec function,
# telling us that there are events after the last censor time
#Due to our approach described in the paper, such events are not of interest
new.data.1<-new.data.1%>%
  filter(!(V.upper==Inf))

#Second, remove intervals for which N_{I_j}<1 as per equation 1
#This is, remove the final interval if it is known that no events is possible
#This prevents superfluous midpoint p-values of 0.5
#Or in the case of randomised p-values, it prevents superfluous noise
new.data.1<-new.data.1%>%
  filter( !((N.risk-Censors.I)==0)  )



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








#Summarising by calculating p-values----------------------------------------------------------------

#Now we look at summarizing the specified intervals


if(p_method=='mid'){
  
  new.data.3<-new.data.2 %>%
    mutate(expected_E=p*(N.risk-Censors.I) )%>%    #first add the expect number of events under the fitted model
    group_by(V.lower) %>%
    summarise(V.upper=unique(V.upper),
              Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
              Observed.E_over.V=sum(Events.obs.I),
              
              V.pval= {0.5*PoissonBinomial::ppbinom(
                x=as.integer(sum(Events.obs.I)),
                probs=(rep(p, (N.risk-Censors.I) )))+
                  0.5*PoissonBinomial::ppbinom(
                    x=as.integer(sum(Events.obs.I)-1),
                    probs=(rep(p,(N.risk-Censors.I) )))},
              
              N.risk.at.V.start=max(N.risk),
              E=sum(Events.obs.I))
} 




if(p_method=='rand'){
  #The only difference in the code compared to p_method=='mid' is how V.pval is defined
  
  #Specifically, we now use randomised p-values 
  #As per Rubin-Delanchy 2019 (Meta-Analysis of Mid-p-Values), randomised Ps can be defined as
  # R = X*P1 + (1-X)*P2      (aside: the direction of the Ps is flipped compared to Delanchy)
  # where X is a uniform 0, 1 r.v
  #rearranging we have: R= X(P1-P2)+P2
  #This is the form used for obtaining randomised p-values when specifying V.pval in our code,
  #and this form is used to avoid storing X values or other code complications 
  
  new.data.3<-new.data.2 %>%
    mutate(expected_E=p*(N.risk-Censors.I) )%>%    #first add the expect number of events under the fitted model
    group_by(V.lower) %>%
    summarise(V.upper=unique(V.upper),
              Expect.E_over.V=sum(expected_E),             #summary statistics for each interval V
              Observed.E_over.V=sum(Events.obs.I),
              
              V.pval= { runif(1)*
                  (PoissonBinomial::ppbinom(x=as.integer(sum(Events.obs.I)),
                                            probs=(rep(p, (N.risk-Censors.I) )))  -
                     PoissonBinomial::ppbinom(x=as.integer(sum(Events.obs.I)-1),
                                              probs=(rep(p, (N.risk-Censors.I) )))) +
                  PoissonBinomial::ppbinom(x=as.integer(sum(Events.obs.I)-1),
                                           probs=(rep(p, (N.risk-Censors.I) )))  },
              
              N.risk.at.V.start=max(N.risk),
              E=sum(Events.obs.I))
  
}





#Bonferroni and individual test results-----------------------------------------
#for the specified intervals

information<-new.data.3 %>%
  mutate(individual.test=ifelse(V.pval<=0.025| V.pval>=0.975,
                                'Reject', 'Accept')) %>%
  mutate(bonferroni.test=ifelse(V.pval<=         
                                  0.025/nrow(new.data.3) |
                                V.pval>=
                                  1- 0.025/nrow(new.data.3),
                                'Reject',
                                'Accept')) %>%
  select(-E) %>%
  mutate(fitted.dist=Distribution)


return(information)

#end of function------------------------------------------------------------
#-----------------------------------------------------------------------------
}

