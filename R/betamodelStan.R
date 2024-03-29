#'Zero-and-one inflated Beta model'
#'@export
betamod<-function(count,mgmt,hab,species=c("YCHUB","BSHINER")){
  #Error bounds
  if(length(species)>1|missing(species))stop("'species' must contain only one value",call.=FALSE)
  if(missing(species))stop("must specify species",call.=FALSE)
  if(missing(count))stop("must specify count data",call.=FALSE)
  if(missing(mgmt))stop("must specify mgmt (management) activity data",call.=FALSE)
  if(missing(hab))stop("must specify hab (habitat) activity data",call=FALSE)
  options(warn=-1)

  if(species=="YCHUB"){
  #Add pname field to count, mgmt, and hab .csv files
  #pname is a numeric wetland pond field
  #yr is a numeric year field

  count = count %>% filter(include == 1)
  mgmt = mgmt %>% filter(include == 1)
  hab = hab %>% filter(include == 1)

  count$pname<-as.numeric(as.factor(count$pond_name))
  count$yr<-as.numeric(as.factor(count$year))
  mgmt$pname<-as.numeric(as.factor(mgmt$pond_name))
  mgmt$yr<-as.numeric(as.factor(mgmt$year))
  hab$pname<-as.numeric(as.factor(hab$pond_name))
  hab$yr<-as.numeric(as.factor(hab$year))
  count$YC_small[is.na(count$YC_small)] <- 0
  count$BS_small[is.na(count$BS_small)] <- 0
  count$YT_small[is.na(count$YT_small)] <- 0

  #Reorganize count data by Site, Wetland Pond, Year, and Species

    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,site=count$site,y=count$YC_small)
    newdat<-newdat %>% group_by(pname,year,site) %>% summarise(y=mean(y))
    countdata<-newdat[order(newdat$year,newdat$pname,newdat$site),]


  #Define 4-dimensional array dimensions
  #nsite = the twenty sites within a wetland pond
  #npond= 11 wetland ponds
  #nday = three survey days
  #nyear = the number of years that each wetland pond was surveyed using this protocol
  nsite=max(countdata$site)
  npond=max(countdata$pname)
  nyear=max(countdata$year)

  #Define and create 4 dimensional array
  #Read observed count data into 4-d array
  #y=array(as.numeric(countdata$y)+.00001,c(nsite,npond,nyear))
  #y[y>1]=0.99
  #Reorganize management activity data by Site, Wetland Pond, and Year
  newmgmt<-data.frame(pname=mgmt$pname,year=mgmt$yr,ychubrem=mgmt$YCHUB_removed,
                      bshinerrem=mgmt$BSHINER_removed,ytoprem=mgmt$YTOP_removed,ychubstock=mgmt$YCHUB_stocked,
                      bshinerstock=mgmt$BSHINER_stocked,ytopstock=mgmt$YTOP_stocked,include=mgmt$include)

  #Reorder management activity data
  mgmtdata<-newmgmt[order(newmgmt$year,newmgmt$pname),]





  #Reorganize habitat data by Site, Wetland Pond, and Year
  newhab<-data.frame(pname=hab$pname,year=hab$yr,site=hab$site,pH=hab$pH,wtemp=hab$wtemp,
                     doxygen=hab$doxygen,wcond=hab$wcond,veg=hab$veg,
                     wdepth=hab$wdepth,include=hab$include)

  #Reorder habitat data
  habdata<-newhab[order(newhab$year,newhab$pname,newhab$site),]
  habdata<-data.table(pname=habdata$pname,year=habdata$year,site=habdata$site,
                      pH=as.numeric(habdata$pH),
                      wtemp=as.numeric(habdata$wtemp),doxygen=as.numeric(habdata$doxygen),
                      wcond=as.numeric(habdata$wcond),
                      veg=as.numeric(habdata$veg),wdepth=as.numeric(habdata$wdepth))

  habdata<-setnafill(habdata, type = "locf")

  newhab<-habdata %>% group_by(pname,year,site) %>% summarise(veg=mean(veg),depth=mean(wdepth),temp=mean(wtemp),oxy=mean(doxygen),ph=mean(pH))


  print("Initiate Bayesian formulation of the Beta regression using Stan.",quote=FALSE)

  model = stan_model("zoibeta.stan")

  nsite=length(countdata$y)
  veg=scale(newhab$veg)
  veg=veg[1:nsite]
  depth=scale(newhab$depth)
  depth=depth[1:nsite]
  temp=scale(newhab$temp)
  temp=temp[1:nsite]
  ph=scale(newhab$ph)
  ph=ph[1:nsite]
  oxy=scale(newhab$oxy)
  oxy=oxy[1:nsite]

  stan_data = list(n = length(countdata$y),
                   y = countdata$y,
                   veg = veg,
                   depth = depth,
                   temp = temp,
                   oxy = oxy,
                   ph = ph)

  fit = sampling(
    model,
    data = stan_data,
    thin = 1,iter=5000,
    verbose = FALSE
  )

  #plot(fit, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
  #plot(fit, plotfun = "trace", pars = c("coef_m","coef_p"), inc_warmup = TRUE)
  #plot(fit, plotfun = "rhat") + ggtitle("Example of adding title to plot")

    #Capture and Write results to working directory (R Data Files)
    #saveRDS(fit,"YCHUB_rds_file")
    print("//////////////////////////////////////////////",quote=FALSE)
    print("//////////////////////////////////////////////",quote=FALSE)
    print("Save these results to your working directory directly or by cut-and-paste into a txt file.",quote=FALSE)

    print(fit, pars=c("mu_coef"), probs=c(.1,.5,.9))

    print("mu_coef[1] = Intercept,
          mu_coef[2] = Vegetation,
          mu_coef[3] = Depth,
          mu_coef[4] = Water Temperature",quote=FALSE)

    #fit_summary=summary(fit)
    #write.csv(print(fit_summary$fit),"PSmall_YCHUB_HabitatParameters.csv",row.names=F)

    print(stan_plot(fit, pars=c("mu_coef"), include = TRUE))
    ggsave("PSmall_YCHUB.png",width = 5, height = 4, dpi = 300, units = "in")

  }else if(species=="BSHINER"){

    count = count %>% filter(include == 1 & onrefuge == 1)
    mgmt = mgmt %>% filter(include == 1 & onrefuge == 1)
    hab = hab %>% filter(include == 1 & onrefuge  == 1)

    count$pname<-as.numeric(as.factor(count$pond_name))
    count$yr<-as.numeric(as.factor(count$year))
    mgmt$pname<-as.numeric(as.factor(mgmt$pond_name))
    mgmt$yr<-as.numeric(as.factor(mgmt$year))
    hab$pname<-as.numeric(as.factor(hab$pond_name))
    hab$yr<-as.numeric(as.factor(hab$year))
    count$YC_small[is.na(count$YC_small)] <- 0
    count$BS_small[is.na(count$BS_small)] <- 0
    count$YT_small[is.na(count$YT_small)] <- 0



    #Reorganize count data by Site, Wetland Pond, Year, and Species

    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,site=count$site,y=count$YC_small)
    newdat<-newdat %>% group_by(pname,year,site) %>% summarise(y=mean(y))
    countdata<-newdat[order(newdat$year,newdat$pname,newdat$site),]


    #Define 4-dimensional array dimensions
    #nsite = the twenty sites within a wetland pond
    #npond= 11 wetland ponds
    #nday = three survey days
    #nyear = the number of years that each wetland pond was surveyed using this protocol
    nsite=max(countdata$site)
    npond=max(countdata$pname)
    nyear=max(countdata$year)

    #Define and create 4 dimensional array
    #Read observed count data into 4-d array
    #y=array(as.numeric(countdata$y)+.00001,c(nsite,npond,nyear))
    #y[y>1]=0.99
    #Reorganize management activity data by Site, Wetland Pond, and Year
    newmgmt<-data.frame(pname=mgmt$pname,year=mgmt$yr,ychubrem=mgmt$YCHUB_removed,
                        bshinerrem=mgmt$BSHINER_removed,ytoprem=mgmt$YTOP_removed,ychubstock=mgmt$YCHUB_stocked,
                        bshinerstock=mgmt$BSHINER_stocked,ytopstock=mgmt$YTOP_stocked,include=mgmt$include)

    #Reorder management activity data
    mgmtdata<-newmgmt[order(newmgmt$year,newmgmt$pname),]





    #Reorganize habitat data by Site, Wetland Pond, and Year
    newhab<-data.frame(pname=hab$pname,year=hab$yr,site=hab$site,pH=hab$pH,wtemp=hab$wtemp,
                       doxygen=hab$doxygen,wcond=hab$wcond,veg=hab$veg,
                       wdepth=hab$wdepth,include=hab$include)

    #Reorder habitat data
    habdata<-newhab[order(newhab$year,newhab$pname,newhab$site),]
    habdata<-data.table(pname=habdata$pname,year=habdata$year,site=habdata$site,
                        pH=as.numeric(habdata$pH),
                        wtemp=as.numeric(habdata$wtemp),doxygen=as.numeric(habdata$doxygen),
                        wcond=as.numeric(habdata$wcond),
                        veg=as.numeric(habdata$veg),wdepth=as.numeric(habdata$wdepth))

    habdata<-setnafill(habdata, type = "locf")

    newhab<-habdata %>% group_by(pname,year,site) %>% summarise(veg=mean(veg),depth=mean(wdepth),temp=mean(wtemp),oxy=mean(doxygen),ph=mean(pH))


    print("Initiate Bayesian formulation of the Beta regression using Stan.",quote=FALSE)

    model = stan_model("zoibeta.stan")

    nsite=length(countdata$y)
    veg=scale(newhab$veg)
    veg=veg[1:nsite]
    depth=scale(newhab$depth)
    depth=depth[1:nsite]
    temp=scale(newhab$temp)
    temp=temp[1:nsite]
    ph=scale(newhab$ph)
    ph=ph[1:nsite]
    oxy=scale(newhab$oxy)
    oxy=oxy[1:nsite]

    stan_data = list(n = length(countdata$y),
                     y = countdata$y,
                     veg = veg,
                     depth = depth,
                     temp = temp,
                     oxy = oxy,
                     ph = ph)

    fit = sampling(
      model,
      data = stan_data,
      thin = 1,iter=5000,
      verbose = FALSE
    )


    #plot(fit, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
    #plot(fit, plotfun = "trace", pars = c("coef_m","coef_p"), inc_warmup = TRUE)
    #plot(fit, plotfun = "rhat") + ggtitle("Example of adding title to plot")

      #Capture and Write results to working directory (R Data Files)
      #saveRDS(fit,"BSHINER_rds_file")
      print("//////////////////////////////////////////////",quote=FALSE)
      print("//////////////////////////////////////////////",quote=FALSE)
      print("Save these results to your working directory directly or by cut-and-paste into a txt file.",quote=FALSE)

      print(fit, pars=c("mu_coef"), probs=c(.1,.5,.9))
      #fit_summary=summary(fit,pars=c("mu_coef"))
      #write.csv(print(fit_summary$summary),"PSmall_BSHINER_HabitatParameters.csv",row.names=F)
      print("mu_coef[1] = Intercept,
          mu_coef[2] = Vegetation,
          mu_coef[3] = Depth,
          mu_coef[4] = Water Temperature",quote=FALSE)

      print(stan_plot(fit, pars=c("mu_coef"), include = TRUE))
      ggsave("PSmall_BSHINER.png",width = 5, height = 4, dpi = 300, units = "in")

    }

  print("Consult with the Regional Statistician (Dr. David R. Stewart) and the Regional Data Management Team once complete (if needed) or if you have any concerns.")
}
