



#' @title Function to describe the light levels inside the canopy
#' @param meteo_hourly Hourly weather data frame with at least the column time (time in numeric, for example 0 to 23),Tair (air temperature in degree C), RH (humidity in pc from 0 to 100), cs the CO2 concentration and PFD the total PFD in micro mol m-2 s-1.
#' @param lat Latitude of the canopy to model.
#' @param DOY Day of Year .
#' @param nlayers Number of layers inside the canopy (max = 50).
#' @param dLAI LAI of each one of the n layers of vegetation in the canopy.
#' @param LAI Cumulated LAI in the midle of each layer.
#' @param Rho Leaf reflectance in the visible wavelengths.
#' @param Tau Leaf transmittance in the visible wavelengths.
#' @param Rho_soil Soil reflectance in the visible wavelengths.
#' @param Rho_NIR Leaf reflectance in the NIR wavelengths.
#' @param Tau_NIR Leaf transmittance in the NIR wavelengths.
#' @param Rho_soil_NIR Soil reflectance in the NIR wavelengths.
#' @param chil Index of departure of the leaf angles from a spherical distribution. -0.4 < chil < 0.6.
#' @param clumpfac Clumping factor, index of non random spatial distribution of leaves. = 1 for randomly spaced leaves, <1 for clumed leaves (Chen et al. 2012).

#' @param model Model for the radiation interception model, default is Norman (only Norman implemented so far).
#'
#' @return
#' @export
#'
#' @examples
#' ##Simulation of the vegetation
#' LAItot = 6
#' nlayers=20
#' dLAI=rep(6/nlayers,nlayers)
#' LAI=cumsum(dLAI)-dLAI/2 # LAI in the midle of each layer
#'##Simulation of the weather
#' meteo_hourly=data.frame(time=0:23,RH=80,Tair=25,cs=400,PFD=dnorm(x = seq(0,23,1),mean = 12,sd = 2.5)/0.16*2000,Tleaf=25)
#' ##Simulation of position and moment of the simulation
#' lat=9.2801048
#' DOY = 60
#' ##Representation of the light interception inside the canopy
#' canopy=f.canopy.interception(meteo_hourly=meteo_hourly,lat = lat,DOY = DOY,nlayers = nlayers,dLAI = dLAI,LAI=LAI)
f.canopy.interception=function(meteo_hourly,lat,DOY,nlayers,dLAI,LAI,Rho=0.11,Tau=0.06,Rho_soil=0.1,Rho_NIR=0.46,Tau_NIR=0.33,Rho_soil_NIR=0.33,chil=0.32,clumpfac=0.85,model='Norman'){
  if(is.null(meteo_hourly$NIR)){meteo_hourly$NIR=meteo_hourly$PFD/4.57}
  time=meteo_hourly$time
  ##Calculation of cosz according to Miguel et al. 2009
  phi=lat*pi/180
  delta=-23.5*cos(360*(DOY+10)/365*pi/180)*pi/180
  h=15*(time-12)*pi/180
  cosz=pmax(0.01,sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(h))
 
  ## Calculation of the proportion of direct and diffuse light using the empirical equation from CLM5 (33.7)
  a0=0.17639
  a1=0.00380
  a2=-9.0039*10^-6
  a3=8.1351*10^-9
  PFD_W=meteo_hourly$PFD/4.57 # Conversion from micro mol m-2 s-1 to W m-2 to be consistent with CLM5 eqn
  prop_dir=pmax(0.01,pmin(0.99,a0+a1*PFD_W+a2*(PFD_W)^2+a3*(PFD_W)^3))
  PFD_dir=prop_dir*meteo_hourly$PFD
  PFD_dif=meteo_hourly$PFD-PFD_dir
  
  b0=0.29548
  b1=0.00504
  b2=-1.4957*10-5
  b3=1.4881*10-8
  
  prop_NIR_dir=pmax(0.01,pmin(0.99,b0+b1*meteo_hourly$NIR+b2*(meteo_hourly$NIR)^2+b3*(meteo_hourly$NIR)^3))
  NIR_dir=prop_NIR_dir*meteo_hourly$NIR
  NIR_dif=meteo_hourly$NIR-NIR_dir
   
  #########################
  plot(x=time,y=meteo_hourly$PFD,type="l",xlab="Time of the day",ylab=expression(Light~intensity~(mu~mol~m^-2~s^-1)),col='black')
  lines(x=time,y=PFD_dif,col="blue")
  lines(x=time,y=PFD_dir,col='red')
  legend('topleft',c('Total','Diffuse light','Direct light'),col=c('black','blue','red'),lty=c(1,1,1))
  ### Creation of matrices with 50 vertical layers and 24 hours
  Canopy_time_dir=Canopy_time_dif=Canopy_time_NIR_dir=Canopy_time_NIR_dif=Canopy_time_tot=Photosynthesis_rate_dir=Photosynthesis_rate_dif=f_sun=f_shade=Temp_leaf_dir=Temp_leaf_dif=gs_dir=gs_dif=matrix(data = NA,nrow = nlayers,ncol = length(time),dimnames = list(Layer=1:nlayers,time=time))
  for(i in 1:length(time)){
    Light_Profile=f.Norman.Radiation(PARdir = PFD_dir[i],PARdif = PFD_dif[i], dLAI = dLAI,nlayers = nlayers,cosz = cosz[i],chil=chil,clumpfac = clumpfac,Rho = Rho,Tau = Tau,Rho_soil_dif =Rho_soil,Rho_soil_dir=Rho_soil)
    Light_Profile_NIR=f.Norman.Radiation(PARdir = NIR_dir[i],PARdif = NIR_dif[i], dLAI = dLAI,nlayers = nlayers,cosz = cosz[i],chil=chil,clumpfac = clumpfac,Rho =Rho_NIR,Tau=Tau_NIR,Rho_soil_dif =Rho_soil_NIR,Rho_soil_dir=Rho_soil_NIR)
    Canopy_time_dir[,i]=(Light_Profile$PARsun)
    Canopy_time_dif[,i]=(Light_Profile$PARsha)
    Canopy_time_NIR_dir[,i]=(Light_Profile_NIR$PARsun)
    Canopy_time_NIR_dif[,i]=(Light_Profile_NIR$PARsha)
    f_sun[,i]=(Light_Profile$fracsun)
    f_shade[,i]=(Light_Profile$fracsha)
  }
  
  ## VIS
  Light=Canopy_time_dir*f_sun+Canopy_time_dif*(1-f_sun)
  figure_light_dir=melt(Canopy_time_dir)
  print(ggplot(data=figure_light_dir,aes(x=time,y=Layer,fill=value))
       +scale_y_reverse()
       +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
       +labs(fill=expression(Direct~light~(mu~mol~m^-2~s^-1))))
  
  figure_light_dif=melt(Canopy_time_dif)
  print(ggplot(data=figure_light_dif,aes(x=time,y=Layer,fill=value))
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(Diffuse~light~(mu~mol~m^-2~s^-1))))
  
  figure_f_sun=melt(f_sun)
  print(ggplot(data=figure_f_sun,aes(x=time,y=Layer,fill=value))+geom_raster()
   +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
    +ggtitle('Fraction of leaves in direct light')
  +labs(fill='f_sun %'))
  
  figure_light_tot=melt(Light)
  print(ggplot(data=figure_light_tot,aes(x=time,y=Layer,fill=value))
        +ggtitle('Mean PFD of an average leaf')
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(PFD~(mu~mol~m^-2~s^-1))))
  
  ## NIR
  NIR=Canopy_time_NIR_dir*f_sun+Canopy_time_NIR_dif*(1-f_sun)
  figure_NIR_dir=melt(Canopy_time_NIR_dir)
  print(ggplot(data=figure_NIR_dir,aes(x=time,y=Layer,fill=value))
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(Direct~NIR~(W~m^-2))))
  
  figure_NIR_dif=melt(Canopy_time_NIR_dif)
  print(ggplot(data=figure_NIR_dif,aes(x=time,y=Layer,fill=value))
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(Diffuse~NIR~(W~m^-2))))
  
  figure_NIR_tot=melt(NIR)
  print(ggplot(data=figure_NIR_tot,aes(x=time,y=Layer,fill=value))
        +ggtitle('Mean NIR of an average leaf')
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(NIR~(W~m^-2))))
  
  return(list(Canopy_time_dir=Canopy_time_dir,Canopy_time_tot=Canopy_time_tot,Canopy_time_dif=Canopy_time_dif,Canopy_time_NIR_dir=Canopy_time_NIR_dir,Canopy_time_NIR_dif=Canopy_time_NIR_dif,f_sun=f_sun,f_shade=f_shade,Light_Profile=Light_Profile,LAItot=sum(dLAI),LAI=LAI))
}



#' @title Gradients of photosynthetic parameters
#' @description Several versions of gradients can be found in the litterature, see for example Lloyd et al. 2010 (Fig. 10 and equation A2), but also the equation A14 from Krinner et al. 2005 and the equation 33 from Clark et al. 2011
#' The simpler model describing the gradients is Vcmax(LAI)=Vcmax0 x exp(-kn x LAI) with Vcmax0 Vcmax at the top of the canopy
#' kn can be also calculated as a function of Vcmax0: kn=exp(alpha x Vcmax0+beta).
#' If kn is NULL, then the function will use the default alpha and beta to calculate kn. If, on the contrary, kn is given, this specific one will be used to calculate the gradients.
#' Krinner et al use a slightly different version of this equation with the parameter lambda: Vcmax(LAI)=Vcmax0 x (1-lambda x (1-exp(-kn*LAI))). The previous equation is a particular case of this one for lambda = 1.
#' @param alpha Slope of the relationship between Vcmax0 and log(kn), see Lloyd et al. 2010.
#' @param beta Intercept of the relationship between Vcmax0 and log(kn), see Lloyd et al. 2010.
#' @param Vcmax0 Vcmax at 25 degree C at the top of the canopy.
#' @param LAI Vector of Leaf Area Index (or depth within the canopy see Clark et al. 2011).
#' @param kn Exponential decrease.
#' @param lambda Asymptot of the decrease (see Krinner et al. 2005).
#' @references Krinner, G., Viovy, N., de Noblet-Ducoudr?, N., Og?e, J., Polcher, J., Friedlingstein, P., . Prentice, I. C. (2005). A dynamic global vegetation model for studies of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19(1). doi:10.1029/2003gb002199.
#' Clark, D. B., Mercado, L. M., Sitch, S., Jones, C. D., Gedney, N., Best, M. J., . Cox, P. M. (2011). The Joint UK Land Environment Simulator (JULES), model description - Part 2: Carbon fluxes and vegetation dynamics. Geoscientific Model Development, 4(3), 701-722. doi:10.5194/gmd-4-701-2011.
#' Lloyd, J., Pati?o, S., Paiva, R. Q., Nardoto, G. B., Quesada, C. A., Santos, A. J. B., . Mercado, L. M. (2010). Optimisation of photosynthetic carbon gain and within-canopy gradients of associated foliar traits for Amazon forest trees. Biogeosciences, 7(6), 1833-1859. doi:10.5194/bg-7-1833-2010.
#' @return Vector of Vcmax (or any other parameter) at the different LAI specified in the call of the function
#' @export
#'
#' @examples
#' LAI=seq(0,6.2,6.2/49)
#' Vcmax=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70)
#' Vcmax2=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70,lambda=0.7)
#' plot(Vcmax)
#' lines(Vcmax2)
f.VcmaxRef.LAI=function(alpha=0.00963,beta=-2.43,Vcmax0=50,LAI=0:8,kn=NULL,lambda=1){
  if(is.null(kn)){kn=exp(alpha*Vcmax0+beta)
  print(paste('kn is',kn))}
  return(Vcmax0*(1-lambda*(1-exp(-kn*LAI))))
}



#' @title Canopy scale GPP calculation
#' @description Generic function to calculate the GPP within a forest (Here GPP = sum of Anet at the canopy level, so it takes into account the leaf mitochondrial respiration)
#' @param meteo_hourly See f.canopy.interception doc. In addition to the requirement for f.canopy.interception, the leaf temperature has to be informed within the column Tleaf.
#' @param Vcmax_Profile Vector of the values of Vcmax at the reference temperature at each layer of the canopy.
#' @param Jmax_Profile Vector of the values of Jmax at the reference temperature at each layer of the canopy.
#' @param Rd_Profile Vector of the values of Rd at the reference temperature at each layer of the canopy.
#' @param Tp_Profile Vector of the values of Tp at the reference temperature at each layer of the canopy.
#' @param g0_Profile Vector of the values of g0 at the reference temperature at each layer of the canopy.
#' @param g1_Profile Vector of the values of g1 at the reference temperature at each layer of the canopy.
#' @param gsmin Minimum stomatal conductance for water to consider. This value will be used as the minimum conductance value to avoid 0 and negative values obtained from the coupled assimilation and conductance models.
#' @param canopy Description of the canopy interception (see canopy_interception function).
#' @param Patm Atmospheric pressure (used to calculate the transpiration).
#' @param ... Other parameters of the photosynthetic model, without gradients, for example curvature factor, quantum yield.. see the help of f.make.param().
#'
#' @return
#' @export
#'
#' @examples
#' #See vignettes on github
f.GPP<-function(meteo_hourly,Vcmax_Profile,Jmax_Profile,Rd_Profile,Tp_Profile,g0_Profile,g1_Profile,gsmin,canopy,Patm=100,...){
  if(length(Vcmax_Profile)!=nrow(canopy$Canopy_time_dir)){print(paste('Are you sure you want to use',length(Vcmax_Profile),'different Vcmax but ',nrow(canopy$Canopy_time_dir),'vertical canopy layers ?'))}
  VpdL_dir=VpdL_dif=Photosynthesis_rate_dir=Photosynthesis_rate_dif=gs_dir=gs_dif=canopy$Canopy_time_dir
  nlayer=nrow(canopy$Canopy_time_dir)
  g1_min=-1 ## This trick is used to fix gsw to gswmin. 
  for(Layer in 1:nlayer){
    res_dir=f.A(PFD = canopy$Canopy_time_dir[Layer,],
                cs = meteo_hourly[,"cs"],
                Tair = meteo_hourly[,"Tair"]+273.15,
                Tleaf= meteo_hourly[,"Tleaf"]+273.15,
                RH = meteo_hourly[,"RH"],
                param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],
                                     abso=1,...
                ))
    ls.gs=which(res_dir$gs<gsmin)
    res_dir$gs[ls.gs]=gsmin
    
    res_dir$A[ls.gs]=f.A(PFD = canopy$Canopy_time_dir[Layer,],cs = meteo_hourly[,"cs"],Tleaf = meteo_hourly[,"Tleaf"]+273.15,Tair = meteo_hourly[,"Tair"]+273.15,RH = meteo_hourly[,"RH"],param = f.make.param(
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=gsmin,
                                                                                                                                                             g1=g1_min,
                                                                                                                                                             abso=1,...
    ))$A[ls.gs]
    Photosynthesis_rate_dir[Layer,]=res_dir$A
    VpdL_dir[Layer,]=res_dir$ds/1000
    gs_dir[Layer,]=res_dir$gs
    res_dif=f.A(PFD = canopy$Canopy_time_dif[Layer,],
                cs = meteo_hourly[,"cs"],
                Tair = meteo_hourly[,"Tair"]+273.15,
                Tleaf= meteo_hourly[,"Tleaf"]+273.15,
                RH = meteo_hourly[,"RH"],
                param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],
                                     abso=1,...
                ))
    ls.gs=which(res_dif$gs<gsmin)
    res_dif$gs[ls.gs]=gsmin
    res_dif$A[ls.gs]=f.A(PFD = canopy$Canopy_time_dif[Layer,],cs =meteo_hourly[,"cs"],Tleaf = meteo_hourly[,"Tleaf"]+273.15,Tair = meteo_hourly[,"Tair"]+273.15,RH = meteo_hourly[,"RH"],param = f.make.param(
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=gsmin,
                                                                                                                                                             g1=g1_min,
                                                                                                                                                             abso=1,...
    ))$A[ls.gs]
    Photosynthesis_rate_dif[Layer,]=res_dif$A
    gs_dif[Layer,]=res_dif$gs
    VpdL_dif[Layer,]=res_dif$ds/1000
  }
  Photosynthesis_rate=(Photosynthesis_rate_dir*canopy$f_sun+Photosynthesis_rate_dif*(1-canopy$f_sun))
  figure_photosynthesis=melt(Photosynthesis_rate)
  a=(ggplot(data=figure_photosynthesis,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(A~(mu~mol~m^-2~s^-1))))
  
  Conductance_rate=(gs_dir*canopy$f_sun+gs_dif*(1-canopy$f_sun))
  Trans=(gs_dir*VpdL_dir/Patm*canopy$f_sun+gs_dif*VpdL_dif/Patm*(1-canopy$f_sun))
  figure_conductance=melt(Conductance_rate)
  b=(ggplot(data=figure_conductance,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(g[sw]~(mol~m^-2~s^-1))))
  print(a)
  print(b)
  totalGPP= sum(Photosynthesis_rate*dLAI,na.rm=TRUE)*365*3600*44/10^6
  totalET= sum(Trans*dLAI,na.rm=TRUE)*365*3600*18*10^-3
  print(paste("GPP = ",totalGPP,"g CO2 m-2 Ground Y-1"))
  print(paste("ET = ",totalET,"L H20 m-2 Ground Y-1"))
  return(list(A=Photosynthesis_rate,gs=Conductance_rate,A_dir=Photosynthesis_rate_dir,gs_dir=gs_dir,A_dif=Photosynthesis_rate_dif,gs_dif=gs_dif,Trans=Trans,GPP=totalGPP,ET=totalET,fig_A=a,fig_gs=b))
}

#' @title Canopy scale GPP calculation, with leaf energy budget
#' @description Generic function to calculate the GPP within a forest (Here GPP = sum of Anet at the canopy level, so it takes into account the leaf mitochondrial respiration)
#' @param meteo_hourly See f.canopy.interception
#' @param Vcmax_Profile Vector of the values of Vcmax at the reference temperature at each layer of the canopy
#' @param Jmax_Profile Vector of the values of Jmax at the reference temperature at each layer of the canopy
#' @param Rd_Profile Vector of the values of Rd at the reference temperature at each layer of the canopy
#' @param Tp_Profile Vector of the values of Tp at the reference temperature at each layer of the canopy
#' @param g0_Profile Vector of the values of g0 at the reference temperature at each layer of the canopy
#' @param g1_Profile Vector of the values of g1 at the reference temperature at each layer of the canopy
#' @param gsmin Minimum stomatal conductance for water to consider. This value will be used as the minimum conductance value to avoid 0 and negative values obtained from the coupled assimilation and conductance models
#' @param canopy Description of the canopy interception (see canopy_interception function)
#' @param Patm Atmospheric pressure (used to calculate the transpiration)
#' @param ... Other parameters of the photosynthetic model, without gradients, for example curvature factor, quantum yield.. see the help of f.make.param()
#'
#' @return
#' @export
#'
#' @examples
#' #See github vignettes
f.GPPT<-function(meteo_hourly,Vcmax_Profile,Jmax_Profile,Rd_Profile,Tp_Profile,g0_Profile,g1_Profile,gsmin,canopy,Patm=100,...){
  if(length(Vcmax_Profile)!=nrow(canopy$Canopy_time_dir)){print(paste('Are you sure you want to use',length(Vcmax_Profile),'different Vcmax but ',nrow(canopy$Canopy_time_dir),'vertical canopy layers ?'))}
  VpdL_dir=VpdL_dif=Photosynthesis_rate_dir=Photosynthesis_rate_dif=gs_dir=gs_dif=rd_dir=rd_dif=Tleaf_dir=Tleaf_dif=RHs_dir=RHs_dif=cs_dir=cs_dif=canopy$Canopy_time_dir
  nlayer=nrow(canopy$Canopy_time_dir)
  param=f.make.param()
 g1_min=-1  #This trick is used to fix gsw to gswmin
  for(Layer in 1:nlayer){
    print(paste('Layer',Layer,'of', nrow(canopy$Canopy_time_dir),'layers'))
    res_dir=f.AT(PFD = canopy$Canopy_time_dir[Layer,],
                 NIR= canopy$Canopy_time_NIR_dir[Layer,],
                 ca = meteo_hourly[,"cs"],
                 Tair = meteo_hourly[,"Tair"]+273.15,
                 wind= meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),
                 RHa = meteo_hourly[,"RH"],
                 abso_s=1,
                 param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                      RdRef = Rd_Profile[Layer],
                                      JmaxRef=Jmax_Profile[Layer],
                                      TpRef=Tp_Profile[Layer],
                                      g0=g0_Profile[Layer],
                                      g1=g1_Profile[Layer],abso=1,...
                 ))
    ls.gs=which(res_dir$gs<gsmin)
    res_dir$gs[ls.gs]=gsmin
    res_dir$A[ls.gs]=f.AT(PFD = canopy$Canopy_time_dir[Layer,],NIR = canopy$Canopy_time_NIR_dir[Layer,],ca = meteo_hourly[,"cs"],Tair = meteo_hourly[,"Tair"]+273.15,RHa = meteo_hourly[,"RH"],wind=meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),abso_s=1,param = f.make.param(
                                                                                                                                                      VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                      RdRef = Rd_Profile[Layer],
                                                                                                                                                      JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                      TpRef=Tp_Profile[Layer],
                                                                                                                                                      g0=gsmin,
                                                                                                                                                      g1=g1_min,abso=1,...
    ))$A[ls.gs]
    #((-g0_Profile[Layer])*400*sqrt(f.ds(Tleaf = meteo_hourly[,"tl"]+273.15,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"RH"])/1000)/(1.6*g1_Profile[Layer]))[ls.gs]
    Photosynthesis_rate_dir[Layer,]=res_dir$A
    VpdL_dir[Layer,]=res_dir$ds/1000
    gs_dir[Layer,]=res_dir$gs
    rd_dir[Layer,]=res_dir$Rd
    Tleaf_dir[Layer,]=res_dir$Tleaf
    RHs_dir[Layer,]=res_dir$RHs
    cs_dir[Layer,]=res_dir$cs
    res_dif=f.AT(PFD = canopy$Canopy_time_dif[Layer,],
                 NIR= canopy$Canopy_time_NIR_dif[Layer,],
                 ca = meteo_hourly[,"cs"],
                 Tair = meteo_hourly[,"Tair"]+273.15,
                 wind=meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),
                 RHa = meteo_hourly[,"RH"],
                 abso_s=1,
                 param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                      RdRef = Rd_Profile[Layer],
                                      JmaxRef=Jmax_Profile[Layer],
                                      TpRef=Tp_Profile[Layer],
                                      g0=g0_Profile[Layer],
                                      g1=g1_Profile[Layer],abso=1,...
                 ))
    ls.gs=which(res_dif$gs<gsmin)
    res_dif$gs[ls.gs]=gsmin
    res_dif$A[ls.gs]=f.AT(PFD = canopy$Canopy_time_dif[Layer,],NIR = canopy$Canopy_time_NIR_dif[Layer,],ca = meteo_hourly[,"cs"],Tair = meteo_hourly[,"Tair"]+273.15,RHa = meteo_hourly[,"RH"],wind=meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),abso_s=1,param = f.make.param(
                                                                                                                                                      VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                      RdRef = Rd_Profile[Layer],
                                                                                                                                                      JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                      TpRef=Tp_Profile[Layer],
                                                                                                                                                      g0=gsmin,
                                                                                                                                                      g1=g1_min,abso=1,...
    ))$A[ls.gs]
    Photosynthesis_rate_dif[Layer,]=res_dif$A
    gs_dif[Layer,]=res_dif$gs
    rd_dif[Layer,]=res_dif$Rd
    VpdL_dif[Layer,]=res_dif$ds/1000
    Tleaf_dif[Layer,]=res_dif$Tleaf
    RHs_dif[Layer,]=res_dif$RHs
    cs_dif[Layer,]=res_dif$cs
  }
  Photosynthesis_rate=(Photosynthesis_rate_dir*canopy$f_sun+Photosynthesis_rate_dif*(1-canopy$f_sun))
  figure_photosynthesis=melt(Photosynthesis_rate)
  a=(ggplot(data=figure_photosynthesis,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab("Vertical level (0= top, 50 = ground)")
     +labs(fill=expression(A~(mu~mol~m^-2~s^-1))))
  
  Conductance_rate=(gs_dir*canopy$f_sun+gs_dif*(1-canopy$f_sun))
  Trans=(gs_dir*VpdL_dir/Patm*canopy$f_sun+gs_dif*VpdL_dif/Patm*(1-canopy$f_sun))
  Tleaf=(Tleaf_dir*canopy$f_sun+Tleaf_dif*(1-canopy$f_sun))
  figure_conductance=melt(Conductance_rate)
  b=(ggplot(data=figure_conductance,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(g[sw]~(mol~m^-2~s^-1))))
  figure_Tleaf=melt(Tleaf)
  c=(ggplot(data=figure_Tleaf,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(Tleaf~(K))))
  print(a)
  print(b)
  print(c)
  totalGPP= sum(Photosynthesis_rate*dLAI,na.rm=TRUE)*365*3600*44/10^6
  totalET= sum(Trans*dLAI,na.rm=TRUE)*365*3600*18*10^-3
  print(paste("GPP = ",totalGPP,"g CO2 m-2 Ground Y-1"))
  print(paste("ET = ",totalET,"L H20 m-2 Ground Y-1"))
  return(list(A=Photosynthesis_rate,gs=Conductance_rate,A_dir=Photosynthesis_rate_dir,gs_dir=gs_dir,A_dif=Photosynthesis_rate_dif,gs_dif=gs_dif,Tleaf_dir=Tleaf_dir,Tleaf_dif=Tleaf_dif,Tleaf=Tleaf,Rd_dir=rd_dir,Rd_dif=rd_dif,Trans=Trans,GPP=totalGPP,ET=totalET,VpdL_dif=VpdL_dif,VpdL_dir=VpdL_dir,RHs_dif=RHs_dif,RHs_dir=RHs_dir,cs_dif=cs_dif,cs_dir=cs_dir,fig_A=a,fig_gs=b,fig_Tleaf=c))
}


