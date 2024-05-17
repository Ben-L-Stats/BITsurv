#tWO Functions to obtain approximate overall test statistics:
# PAVSI and TFT




#' Protection Against Very Small Intervals (PAVSI)
#'
#' This function produces an overall p-value using the PAVSI approach.
#'
#' @param p.vals A vector of p-values
#' @param print A simple command of whether to return the results of the test under H0.
#' Default is to return the result (print=TRUE).
#'
#' @return An overall p-value
#' @export
#'
#' @examples
#'  #In practice, this will just be
#'  BIT.TS.PAVSI(p.vals=BIT.table$V.mid.pval)
#'  #where BIT.table is the dataframe returned by the BIT.surv function.
#'
#'  #An example that can be run immediately
#'  BIT.TS.PAVSI(p.vals=0.1*c(1:9))
#'
BIT.TS.PAVSI<-function(p.vals, print=TRUE){

rejects<-sum(p.vals<=0.025)+   #num of intervals that fail the individual test
          sum(p.vals>=0.975)

size=length(p.vals)            #number of intervals

#calculate the midpoint p-value
overall.p<-0.5*
    (pbinom(q=rejects-1, size=size, prob=0.05,lower.tail=FALSE)+   #Prob(T>=rejects)=Prob(T>rejects-1)
     pbinom(q=rejects, size=size, prob=0.05, lower.tail=FALSE))    #Prob(T>rejects)
#Note that 0.05 is probability of failing an individual test under the null

if(print==TRUE){
if(overall.p<0.05)
{
  print('Reject H0: Overall PAVSI p-value<0.05')
} else{
  print('Accept H0: Overall PAVSI p-value>=0.05')
}
}

overall.p
}











#----------------------------------------------------------------------------------


#New function






#' Transformed Fisher Test (TFT)
#'
#' This function produces an overall p-value using the TFT approach.
#'
#' @param p.vals A vector of p-values
#' @param print A simple command of whether to return the results of the test under H0.
#' Default is to return the result (print=TRUE).
#'
#' @return An overall p-value
#' @export
#'
#' @examples
#'  #In practice, this will just be
#'  BIT.TS.TFT(p.vals=BIT.table$V.mid.pval)
#'  #where BIT.table is the dataframe returned by the BIT.surv function.
#'
#'  #An example that can be run immediately
#'  BIT.TS.TFT(p.vals=0.1*c(1:9))
#'
BIT.TS.TFT<-function(p.vals, print=TRUE){

  #number of intervals
  size=length(p.vals)

  #apply the transformation
  df<-data.frame(p.vals) %>%
    mutate(u.prime=ifelse(p.vals<=0.5,
                          2*p.vals,
                          2*(1-p.vals)))

  #obtain the test statistic
  Trans.Fisher.TS<- -2*sum(log(df$u.prime))

  #obtain the overall p-value
  #no mid-p here as this is a continuous test statistic
  overall.p<-pchisq(Trans.Fisher.TS, 2*size, lower.tail = FALSE)  #=Prob(T>t)

  if(print==TRUE){
  if(overall.p<0.05)
  {
    print('Reject H0: Overall TFT p-value<0.05')
  } else{
    print('Accept H0: Overall TFT p-value>=0.05')
  }
  }

  overall.p
}

