#Fit exponential model to single arm data and produce plot of this

#Input is:
# surv.data
# Distribution
# break.time

# surv.data with numeric time and event values

#and requires libraries
# library(survHE) #Load survHE. Loads flexsurv and survival automatically
# library(survminer) #Required for ggsurvplot
# library(dplyr)


#' Fits a parametric survival curve to data and plots the results
#'
#' This function produces a plot
#'
#' @param surv.data A dataframe with 2 columns labelled time and event.
#' @param Distribution Pick your parametric survival distribution. Available options are:
#' 'exp', 'weibull', 'gompertz', 'llogis','lnorm', 'gamma', and 'gengamma'. See
#' flexsurv for further details on these distributions.
#' @param break.time Gives the labelled breaks in time on the survival plot and number and risk plot.
#'
#' @return A plot: KM curve (in purple), an overlaid  fitted parametric curve (in black), and
#' number at risk.
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
#' #Produce the plot
#' Fit.curve.plot(surv.data=simple.data,
#'                Distribution='exp',
#'                break.time=5)
#'
Fit.curve.plot<-function(surv.data, Distribution, break.time=10){
  #Begin function----------------------------------------------------------------
  #------------------------------------------------------------------------------

#Fit the model to the data
par.est<-flexsurv::flexsurvreg(survival::Surv(time, event) ~ 1, data = surv.data,
                     dist=Distribution)

#define the survival data y, for the fitted curve
#x is specified in as
#the resolution and boundaries of the fitted curve
x<-round(max(surv.data$time))*(1:1000)/1000   #1000 times up to the last observed time
#number at risk itervals are defaulted to every 10 time units


#note that the transoformation of the par.est values here are either exp(..) or
#the identity. The transform to use is specified in
#par.est$dlist$transforms

#Obtain survival values for the fitted curves-----------------------------------
if (Distribution=='exp')    {y<-1-pexp(q=x,
                                     rate=exp(par.est$coefficients))}
if (Distribution=='weibull'){y<-1-pweibull(q=x,
                                         shape=exp(par.est$coefficients['shape']),
                                         scale=exp(par.est$coefficients['scale']))}
if (Distribution=='gompertz'){y<-1-flexsurv::pgompertz(q=x,
                                         shape=par.est$coefficients['shape'],
                                         rate=exp(par.est$coefficients['rate']))}
if (Distribution=='llogis'){y<-1-flexsurv::pllogis(q=x,
                                         shape=exp(par.est$coefficients['shape']),
                                         scale=exp(par.est$coefficients['scale']))}
if (Distribution=='lnorm'){y<-1-plnorm(q=x,
                                      meanlog=par.est$coefficients['meanlog'],
                                      sdlog=exp(par.est$coefficients['sdlog']))}
if (Distribution=='gamma'){y<-1-pgamma(q=x,
                                       shape=exp(par.est$coefficients['shape']),
                                       rate=exp(par.est$coefficients['rate']))}
if (Distribution=='gengamma'){y<-1-flexsurv::pgengamma(q=x,
                                      mu=par.est$coefficients['mu'],
                                      sigma=exp(par.est$coefficients['sigma']),
                                      Q=par.est$coefficients['Q'])}

#Plotting-------------------------------------------------------------------------------


#get data in appropriate form for initial plot
KM.est <- survival::survfit(survival::Surv(time,event) ~ 1,
                            data = surv.data)
#initial plot
plot<-survminer::ggsurvplot(KM.est, data = surv.data,
                 title= paste0('Fitted Survival curve - ', Distribution),
                 risk.table = TRUE,
                 palette = c("purple"),
                 break.time.by = break.time)

#stitch it all together-----------------------------------

#KM curve with fitted curve on top
g1<-ggplot2::ggplotGrob(plot$plot+
         ggplot2::geom_line(data=data.frame(x=x,y=y),
                            ggplot2::aes(x,y)))

#At risk table
g2<-ggplot2::ggplotGrob(plot$table)

min_ncol <- min(ncol(g1), ncol(g2))

g <- gridExtra::gtable_rbind(g1[, 1:min_ncol],
                             g2[, 1:min_ncol],
                             size="last")

g$widths <- grid::unit.pmax(g1$widths, g2$widths)
grid::grid.newpage()



print(grid::grid.draw(g))


#End of function----------------------------------------------------------------
#------------------------------------------------------------------------------
}
