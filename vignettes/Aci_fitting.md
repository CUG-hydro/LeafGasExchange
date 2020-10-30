Aci\_fitting
================
Julien LAMOUR
10/16/2020

# Fitting an Aci curve

The objectif of this tutorial is to show how to fit an Aci curve using
the LeafGasExchange package. In this tutorial, an Aci curve is first
simulated with known photosynthetic parameters and noise. This curve is
then fitted to retrieve the parameters.

## Simulating an Aci curve

For this exemple we first simulate a photosynthesis curve, but it would
work the same if the data were not simulated but measured. The data
simulation is done using the function f.A. This function needs a list of
photosynthetic parameters which are produced using the function
f.make.param() and a list of input variables (CO2 at the surface of the
leaf, leaf temperature, incident light, RH). Too have more information
on the function f.make.param, you can use the command ?f.make.param in R
console.

``` r
param=f.make.param(VcmaxRef = 50,JmaxRef=50*1.7,TpRef=50/10)
CO2=seq(50,1500,50)
Tleaf=30+273.16
Tair=27+273.16
PAR=1800
RH=80
simul=f.A(PFD = PAR,cs = CO2,Tleaf = Tleaf,Tair = Tair,RH = RH,param = param)

# Here we include a normal error with a standard deviation proportionnal to the gross photosynthesis (as it is often observed)

noise=unlist(lapply(X = simul$Ag,FUN = function(x){rnorm(n=1,mean=0,sd=0.035*(x))}))
simul$A=simul$A+noise
measures=data.frame(Tleaf=Tleaf,Ci=simul$ci,PARi=PAR,Photo=simul$A)
```

We display this simulated curve using the function f.plot

``` r
f.plot(measures = measures,type = 'Aci',list_legend = param[c('VcmaxRef','JmaxRef','TpRef','RdRef')],param = param)
```

![](Aci_fitting_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Fitting an Aci curve

To fit an Aci curve, it is necessary to detail the parameter that we
want to estimate. All the parameters present in f.make.param can
potentially be fitted even if it would not always make sense. We do a
first fitting with only the parameters VcmaxRef, JmaxRef and RdRef.
Those parameters have to be given in the list Start, with initial
values. The method will look for different initial values around those
values so it is not necessary to give very good ones, just not too
stupid ones. The photosynthetic parameters have to be given in the list
param. This is used to determine what should be the parameters for the
temperature dependence, for the leaf absorbance, theta, etc. By default,
the equations and parameters used in the TBM FATES to simulate the
photosynthesis are used. In this example, we also give a very high value
to TpRef so the TPU limitation is not considered when fitting the curve.

``` r
fitting1=f.fitting(measures = measures,Start = list(JmaxRef = 30, VcmaxRef = 50, RdRef = 1),param=f.make.param(),modify.init=TRUE,do.plot=TRUE,type='Aci')
```

    ## $par
    ##  JmaxRef VcmaxRef    RdRef 
    ## 78.51244 48.36841  1.16680 
    ## 
    ## $value
    ## [1] 17.35969
    ## 
    ## $counts
    ## function gradient 
    ##      128       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL
    ## 
    ## [1] "sd 0.760694608423429"
    ## Length  Class   Mode 
    ##      1   mle2     S4

![](Aci_fitting_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

In a second example we now also fit TpRef

``` r
fitting2=f.fitting(measures = measures,Start = list(JmaxRef = 30, VcmaxRef = 50, RdRef = 1,TpRef=9),param=f.make.param(),modify.init=TRUE,do.plot=TRUE,type='Aci')
```

    ## $par
    ##   JmaxRef  VcmaxRef     RdRef     TpRef 
    ## 86.628962 50.379243  1.513889  4.973832 
    ## 
    ## $value
    ## [1] 9.426514
    ## 
    ## $counts
    ## function gradient 
    ##      387       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL
    ## 
    ## [1] "sd 0.560550733535672"
    ## Length  Class   Mode 
    ##      1   mle2     S4

![](Aci_fitting_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

The fitting corresponds to a list of 2 objects. The first object
corresponds to the fitting using a minimum square function whereas the
second object corresponds to a maximum likelihood which is made using
the package mle2. This latter method is usefull because it allowes us to
calculate the confidence interval of the parameters. The mean parameters
used to simulate the curves (VcmaxRef = 50,JmaxRef=85,TpRef=5,
RdRef=1.43 and sigma\_b= 0.035) should be in the confidence interval.

``` r
confint(fitting2[[2]])
```

    ##                2.5 %      97.5 %
    ## sigma_b   0.03089982  0.05151112
    ## JmaxRef  82.39958570          NA
    ## VcmaxRef 47.48716414 51.37171767
    ## TpRef     4.79290861  5.09148706
    ## RdRef     1.28516419  1.50849007

It is possible to compare the AIC of the two models using the base
function AIC. The lower AIC corresponds to the best model, showing that
in this case and as expected, the TPU limitation is usefull to improve
the fit of the model.

``` r
AIC(fitting1[[2]])
```

    ## [1] 62.20528

``` r
AIC(fitting2[[2]])
```

    ## [1] 48.57756

It is also possible to calculate the interval of confidence and
prediction of the Aci curve from the outputs of the fitting.

``` r
var_cov=fitting2[[2]]@vcov
random_param=rmvnorm(1000,mean=fitting2[[2]]@coef[c('sigma_b','JmaxRef','VcmaxRef','TpRef','RdRef')],sigma = var_cov)

random_simul=matrix(data = NA,nrow = 1000,ncol = nrow(measures))
for(i in 1:1000){
  random_simul[i,]=f.Aci(PFD = measures$PARi,ci = measures$Ci,Tleaf =measures$Tleaf,param = f.make.param(VcmaxRef=random_param[i,'VcmaxRef'],JmaxRef=random_param[i,'JmaxRef'],TpRef=random_param[i,'TpRef'],RdRef=random_param[i,'RdRef']))$A
}

fit_pred=f.Aci(PFD = measures$PARi,ci = measures$Ci,Tleaf =measures$Tleaf,param = f.make.param(VcmaxRef=fitting2[[2]]@coef[c('VcmaxRef')],JmaxRef=fitting2[[2]]@coef[c('JmaxRef')],TpRef=fitting2[[2]]@coef[c('TpRef')],RdRef=fitting2[[2]]@coef[c('RdRef')]))

sd_mean=apply(X = random_simul,MARGIN = 2,FUN =sd)
sd_res=fitting2[[2]]@coef['sigma_b']*fit_pred$Ag
sd_tot=sqrt(sd_mean^2+sd_res^2)
simul_confint=rbind(fit_pred$A-1.96*sd_mean,fit_pred$A+1.96*sd_mean)
simul_pred=rbind(fit_pred$A-1.96*sd_tot,fit_pred$A+1.96*sd_tot)

f.plot(measures = measures,type = 'Aci',list_legend = param[c('VcmaxRef','JmaxRef','TpRef','RdRef')],param = param)
polygon(c(measures$Ci ,rev(measures$Ci)),c(simul_pred[1,], rev(simul_pred[2,])),
        col=adjustcolor("lightgrey",alpha.f=0.5),border=NA)
polygon(c(measures$Ci ,rev(measures$Ci)),c(simul_confint[1,], rev(simul_confint[2,])),
        col=adjustcolor("#99CC99",alpha.f=0.5),border=NA)
```

![](Aci_fitting_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->