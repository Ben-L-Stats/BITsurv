# R Package: BITsurv
Binomial interval testing for fitted survival models.

## Description
This package uses the methodology described in the following arxiv paper ...

The examples provided in the paper are available in the examples folder.
 
The package facilitates interval checking for the 7 standard parametric survival models. (Exponential, Gamma, Generalised gammma, Gompertz, Log-logistic, Log-normal, and Weibull.)
 
## Installation and use
When installing BITsurv for the first time run the following 4 lines:
```
install.packages('dplyr')         #this installs dplyr
install.packages('devtools')      #provides the install_github function
library(devtools)
install_github('Ben-L-Stats/BITsurv/tree/main/Package')  #this installs BITsurv
```

To use BITsurv, run the following:
```
library(dplyr)
library(BITsurv)
```
Note that the dplyr package must be loaded in order for BITsurv functions to work. 

## Examples
We recommend opening the examples folder and using these to get started with the package. 

## Documentation
The documentation with examples can be accessed in R 
```
#The two main functions are:
?BIT.surv
?BIT.plot

#Additional functions include:
?BIT.TS.PAVSI
?BIT.TS.TFT
?Fit.curve.plot
