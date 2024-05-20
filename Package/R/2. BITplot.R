#Fit exponential model to single arm data and produce plot of this

#Input is:
# surv.data
# x values for the plot

# surv.data with numeric time and event values

#and libraries
# library(survHE) #Load survHE. Loads flexsurv and survival automatically
# library(survminer) #Required for ggsurvplot
# library(dplyr)




#' Takes a BITsurv dataframe and returns a BIT plot
#'
#' This function produces a plot
#'
#' @param surv.data A dataframe with 2 columns labelled time and event.
#' @param Distribution Pick your parametric survival distribution. Available options are:
#' 'exp', 'weibull', 'gompertz', 'llogis','lnorm', 'gamma', and 'gengamma'. See
#' flexsurv for further details on these distributions.
#' @param BIT.table A dataframe which is the output of the BITsurv function.
#' This dataframe also states the parametric survival distribution of interest, where available options are:
#' 'exp', 'weibull', 'gompertz', 'llogis','lnorm', 'gamma', and 'gengamma'.
#' @param break.time Gives the labelled breaks in time on the survival plot and number and risk plot.
#'
#' @return A binomial interval test (BIT) plot
#' @export
#'
#' @examples
#'
#' #The dplyr package is used to make the simple.data
#' simple.data<-data.frame(T=rexp(100, rate=1/10),        #the underlying event process
#'                         cens=runif(100, 0, 20)) %>%    #the underlying censor process
#'      mutate(event=ifelse(T<=cens,0, 1),
#'             time=ifelse(T<=cens,T,cens)) %>%
#'      select(time,event)
#'
#'#Obtain a BIT.table dataframe using the BIT.surv(...) function
#'
#' #Produce the plot
#' BIT.plot(surv.data=simple.data,
#'          BIT.table=BIT.table,
#'          break.time=5)
#'
BIT.plot<-function(surv.data, BIT.table, break.time=10){
  #Begin function----------------------------------------------------------------
  #------------------------------------------------------------------------------

  #Specify the distribution of interest
  Distribution<-unique(BIT.table$fitted.dist)

  #Fit the model to the data
  par.est<-flexsurv::flexsurvreg(survival::Surv(time, event) ~ 1,
                                 data = surv.data,
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
  if (Distribution=='gompertz'){y<-1-pgompertz(q=x,
                                               shape=par.est$coefficients['shape'],
                                               rate=exp(par.est$coefficients['rate']))}
  if (Distribution=='llogis'){y<-1-pllogis(q=x,
                                           shape=exp(par.est$coefficients['shape']),
                                           scale=exp(par.est$coefficients['scale']))}
  if (Distribution=='lnorm'){y<-1-plnorm(q=x,
                                         meanlog=par.est$coefficients['meanlog'],
                                         sdlog=exp(par.est$coefficients['sdlog']))}
  if (Distribution=='gamma'){y<-1-pgamma(q=x,
                                         shape=exp(par.est$coefficients['shape']),
                                         rate=exp(par.est$coefficients['rate']))}
  if (Distribution=='gengamma'){y<-1-pgengamma(q=x,
                                               mu=par.est$coefficients['mu'],
                                               sigma=exp(par.est$coefficients['sigma']),
                                               Q=par.est$coefficients['Q'])}

  #Survival plotting-------------------------------------------------------------------------------


  #get data in appropriate form for initial plot
  KM.est <- survival::survfit(survival::Surv(time,event) ~ 1,
                              data = surv.data)
  #initial plot
  plot<-survminer::ggsurvplot(KM.est,
                              data = surv.data,
                              title= paste0('Fitted Survival curve - ', Distribution),
                              risk.table = TRUE,
                              palette = c("purple"),
                              break.time.by = break.time)

#Interval plotting----------------------------------------------------------------------

  plot.data<-BIT.table %>%
    filter(individual.test=='Reject') %>%
    select(V.lower,V.upper,bonferroni.test,individual.test) %>%
    mutate(rejection.type=ifelse(bonferroni.test=='Reject',
                                 'Bonferroni',
                                 'Individual'))




base.plot<-survminer::ggsurvplot(KM.est,
                                 data = surv.data,
                                 risk.table = TRUE,
                                 linetype= 0,
                                 conf.int = FALSE,
                                 censor = FALSE,
                                 legend.title = '',
                                 legend.labs = '',
                                 break.time.by = break.time,
                                 break.y.by = 1,
                                 ylab=" ")



BIT.plot<-base.plot$plot+
    ggplot2::geom_rect(ggplot2::aes(NULL, NULL,
                           xmin = V.lower, xmax = V.upper,
                           fill = rejection.type),
              ymin = 0, ymax = 1, data = plot.data)  +
   ggplot2::scale_fill_manual(values = ggplot2::alpha(c('Bonferroni'="red",
                                                        'Individual'="darkgrey"), .5))+
   ggplot2::geom_vline(xintercept=unique(c(BIT.table$V.upper,BIT.table$V.lower)),
               colour='black', linetype='dashed')


  #stitch it all together-----------------------------------

  #KM curve with fitted curve on top
  g1<-ggplot2::ggplotGrob(plot$plot+
                   ggplot2::geom_line(data=data.frame(x=x,y=y),
                                      ggplot2::aes(x,y))+
                   ggplot2::geom_vline(xintercept=unique(censors$time),
                              colour='grey', linetype='dashed'))

  #BIT plot
  g2 <- ggplot2::ggplotGrob(BIT.plot)

  #At risk table
  g3<-ggplot2::ggplotGrob(plot$table)

  min_ncol <- min(ncol(g1), ncol(g2), ncol(g3))

  g <- gridExtra::gtable_rbind(g1[, 1:min_ncol],
                               g2[, 1:min_ncol],
                               g3[, 1:min_ncol],
                               size="last")

  g$widths <- grid::unit.pmax(g1$widths, g2$widths, g3$widths)
  grid::grid.newpage()



  print(grid::grid.draw(g))


  #End of function----------------------------------------------------------------
  #------------------------------------------------------------------------------
}
