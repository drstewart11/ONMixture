#'Hierarchical Bayesian models'
#'@export
countmix<-function(count,mgmt,hab,species=c("YCHUB","BSHINER")){
  #Error bounds
  if(length(species)>1|missing(species))stop("'species' must contain
                                             only one value",call.=FALSE)
  if(missing(species))stop("must specify species",call.=FALSE)
  if(missing(count))stop("must specify count data",call.=FALSE)
  if(missing(mgmt))stop("must specify mgmt (management) activity
                        data",call.=FALSE)
  if(missing(hab))stop("must specify hab (habitat) activity data",call=FALSE)
  options(warn=-1)

  #Add pname field to count, mgmt, and hab .csv files
  #pname is a numeric wetland pond field
  #yr is a numeric year field
  count$pname<-as.numeric(as.factor(count$pond_name))
  count$yr<-as.numeric(as.factor(count$year))
  mgmt$pname<-as.numeric(as.factor(mgmt$pond_name))
  mgmt$yr<-as.numeric(as.factor(mgmt$year))
  hab$pname<-as.numeric(as.factor(hab$pond_name))
  hab$yr<-as.numeric(as.factor(hab$year))

  #Reorganize count data by Site, Wetland Pond, Year, and Species
  if(species=="YCHUB"){
    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,
                       site=count$site,y=count$YCHUB)
    countdata<-newdat[order(newdat$year,newdat$pname,newdat$day,newdat$site),]
  }else if(species=="BSHINER"){
    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,
                       site=count$site,y=count$BSHINER)
    countdata<-newdat[order(newdat$year,newdat$pname,newdat$day,newdat$site),]
  }

  #Define 4-dimensional array dimensions
  #nsite = the twenty sites within a wetland pond
  #npond= 11 wetland ponds
  #nday = three survey days
  #nyear = the number of survey years
  nsite=max(countdata$site)
  npond=max(countdata$pname)
  nday=max(countdata$day)
  nyear=max(countdata$year)

  #Define and create 4 dimensional array
  #Read observed count data into 4-d array
  y=array(as.numeric(countdata$y),c(nsite,nday,npond,nyear))

  #Reorganize management activity data by Site, Wetland Pond, and Year
  newmgmt<-data.frame(pname=mgmt$pname,year=mgmt$yr,ychubrem=mgmt$YCHUB_removed,
                      bshinerrem=mgmt$BSHINER_removed,ytoprem=mgmt$YTOP_removed,
                      ychubstock=mgmt$YCHUB_stocked,
                      bshinerstock=mgmt$BSHINER_stocked,
                      ytopstock=mgmt$YTOP_stocked)

  #Reorder management activity data
  mgmtdata<-newmgmt[order(newmgmt$year,newmgmt$pname),]
  
  #Expand data.frame
  mgmtdat<-mgmtdata %>% group_by(pname,year,ychubrem,bshinerrem,ytoprem,
                                 ychubstock,bshinerstock,ytopstock) %>% expand(sites = 1:10)

  #Number of Yaqui chub removed
  ychubrem=scale(mgmtdat$ychubrem)
  ychubrem=array(ychubrem,c(nsite,npond,nyear))

  #Number of beautiful shiner removed
  bshinerem=scale(mgmtdat$bshinerrem)
  bshinerem=array(bshinerem,c(nsite,npond,nyear))

  #Number of Yaqui topminnow removed
  ytoprem=scale(mgmtdat$ytoprem)
  ytoprem=array(ytoprem,c(nsite,npond,nyear))

  #Number of Yaqui chub stocked
  ychubstock=scale(mgmt$YCHUB_stocked)
  ychubstock=array(ychubstock,c(nsite,npond,nyear))

  #Number of beautiful shiner stocked
  bshinestock=scale(mgmtdat$bshinerstock)
  bshinestock=array(bshinestock,c(nsite,npond,nyear))

  #Number of Yaqui topminnow stocked
  ytopstock=scale(mgmtdat$ytopstock)
  ytopstock=array(ytopstock,c(nsite,npond,nyear))


  #Reorganize habitat data by Site, Wetland Pond, and Year
  newhab<-data.frame(pname=hab$pname,year=hab$yr,site=hab$site,
                     pH=hab$pH,wtemp=hab$wtemp,
                     doxygen=hab$doxygen,wcond=hab$wcond,
                     ntu=hab$ntu,algal=hab$algal,veg=hab$veg,
                     wdepth=hab$wdepth)

  #Reorder habitat data
  habdata<-newhab[order(newhab$year,newhab$pname,newhab$site),]
  habdata<-data.table(pH=as.numeric(habdata$pH),
                      wtemp=as.numeric(habdata$wtemp),
                      doxygen=as.numeric(habdata$doxygen),
                      wcond=as.numeric(habdata$wcond),
                      ntu=as.numeric(habdata$ntu),
                      algal=as.numeric(habdata$algal),
                      veg=as.numeric(habdata$veg),
                      wdepth=as.numeric(habdata$wdepth))

  habdata<-setnafill(habdata, type = "locf")

  #Scale habitat variables with mean 0 and input into a 3 dimensional array
  #pH
  pH=scale(habdata$pH)
  pH=array(pH,c(nsite,npond,nyear))

  #Water temperature
  wtemp=scale(habdata$wtemp)
  wtemp=array(wtemp,c(nsite,npond,nyear))

  #Dissolved oxygen
  doxygen=scale(habdata$doxygen)
  doxygen=array(doxygen,c(nsite,npond,nyear))

  #Water conductivity
  wcond=scale(habdata$wcond)
  wcond=array(wcond,c(nsite,npond,nyear))

  #Turbidity
  ntu=scale(habdata$ntu)
  ntu=array(ntu,c(nsite,npond,nyear))

  #Algae concentration
  #algal=scale(habdata$algal)
  #algal=array(algal,c(nsite,npond,nyear))

  #Percent submergent aquatic vegetation
  veg=scale(habdata$veg)
  veg=array(veg,c(nsite,npond,nyear))

  #Water depth
  wdepth=scale(habdata$wdepth)
  wdepth=array(wdepth,c(nsite,npond,nyear))

  print("Initiate Bayesian Population Model. This may take several minutes to hours.",quote=FALSE)

  if(species=="YCHUB"){
  modelFilename="Bayesian.Population.Model.txt"
  cat("
      model{

      phi~dunif(0,100)
      #alpha_veg~dnorm(0,0.5)
      #alpha_depth~dnorm(0,0.5)
      #alpha_temp~dnorm(0,0.5)

      for(k in 1:npond){
      beta[k]~dnorm(0,0.01)
      r[k]~dunif(0,5)
      K[k]~dunif(50,5000)
      eta[k]~dgamma(phi,phi)
      #alpha[k]~dnorm(0,0.5)
      }


      for(i in 1:nsite){
      for(k in 1:npond){
      N[i,k,1]~dpois(lambda[i,k])
      lambda[i,k]<-muL[i,k]*eta[k]
      log(muL[i,k])<-beta[k]

      for(t in 2:nyear){
      mu[i,k,t-1]<-N[i,k,t-1]*exp(r[k]*(1-(N[i,k,t-1]/K[k])))
      #mu[i,k,t-1]<-r[k]*N[i,k,t-1]*((K[k]-N[i,k,t-1])/K[k])
      N[i,k,t]~dpois(mu[i,k,t-1])
      }
      }
      }

      for(i in 1:nsite){
      for(j in 1:nrep){
      for(k in 1:npond){
      for(t in 1:nyear){
      y[i,j,k,t]~dbin(q[i,j,k,t],N[i,k,t])

      #Detection probabilities
      q[i,j,k,t]~dbeta(2.5,6) #Slightly informative prior
      #logit(q[i,j,k,t])<-alpha[k] + alpha_veg*veg[i,k,t] + alpha_depth*wdepth[i,k,t] + alpha_temp*wtemp[i,k,t]
      }}}}

      for(k in 1:npond){
      for(t in 1:nyear){
      N.total[k,t]<-sum(N[,k,t])
      det.mean[k,t]<-mean(q[,,k,t])
      }
      }

      }",fill=TRUE,file=modelFilename)
}else if(species=="BSHINER"){
     modelFilename="Bayesian.Population.Model.txt"
  cat("
      model{

      phi~dunif(0,100)
      #alpha_veg~dnorm(0,0.5)
      alpha_depth~dnorm(0,0.5)
      #alpha_temp~dnorm(0,0.5)

      for(k in 1:npond){
      beta[k]~dnorm(0,0.01)
      r[k]~dunif(0,5)
      K[k]~dunif(50,5000)
      eta[k]~dgamma(phi,phi)
      alpha[k]~dnorm(0,0.5)
      }


      for(i in 1:nsite){
      for(k in 1:npond){
      N[i,k,1]~dpois(lambda[i,k])
      lambda[i,k]<-muL[i,k]*eta[k]
      log(muL[i,k])<-beta[k]

      for(t in 2:nyear){
      mu[i,k,t-1]<-N[i,k,t-1]*exp(r[k]*(1-(N[i,k,t-1]/K[k])))
      #mu[i,k,t-1]<-r[k]*N[i,k,t-1]*((K[k]-N[i,k,t-1])/K[k])
      N[i,k,t]~dpois(mu[i,k,t-1])
      }
      }
      }

      for(i in 1:nsite){
      for(j in 1:nrep){
      for(k in 1:npond){
      #alpha[i,j,k]~dnorm(0,0.5)
      for(t in 1:nyear){
      y[i,j,k,t]~dbin(q[i,j,k,t],N[i,k,t])

      #Detection probabilities
      #q[i,j,k,t]~dbeta(10,9) #Slightly informative prior
      logit(q[i,j,k,t])<-alpha[k] + alpha_depth*wdepth[i,k,t] 
      }}}}

      for(k in 1:npond){
      for(t in 1:nyear){
      N.total[k,t]<-sum(N[,k,t])
      det.mean[k,t]<-mean(q[,,k,t])
      }
      }

      }",fill=TRUE,file=modelFilename)
    }
    
    
    
  #Initial values
  Nst<-apply(y,c(1,3,4),sum,na.rm=T)+1
  jags.inits=function()list(N=Nst)

  #Bundle data
  jags.data=list(y=y,nsite=nsite,nrep=nday,npond=npond,nyear=nyear,veg=veg,wdepth=wdepth,wtemp=wtemp)

  #Parameters monitored
  jags.params=c("phi","det.mean","K","r",
                "N.total")

  #MCMC settings
  nc=4; nt=1; nb=15000; ni=500000

  out<-jags(jags.data,jags.inits,parameters.to.save=jags.params,
            model.file=modelFilename,n.chains=nc,
            n.iter=ni,n.burnin=nb,n.thin=nt,parallel=TRUE,verbose=TRUE)

  #Create Year labels
  yrlab<-seq(min(count$year),max(count$year),by=1)
  yrlab<-sort(rep(yrlab,npond))

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(count$pond_name)))),nyear)

  #Summarize posteriors for abundance
  N.total<-round(unlist(out$mean$N.total))
  N.total<-as.vector(N.total)
  N.lower<-round(unlist(out$q2.5$N.total))
  N.lower<-as.vector(N.lower)
  N.upper<-round(unlist(out$q97.5$N.total))
  N.upper<-as.vector(N.upper)

  #Summarize posteriors for detection
  q.mean<-round(unlist(out$mean$det.mean),2)
  q.mean<-as.vector(q.mean)
  q.lower<-round(unlist(out$q2.5$det.mean),2)
  q.lower<-as.vector(q.lower)
  q.upper<-round(unlist(out$q97.5$det.mean),2)
  q.upper<-as.vector(q.upper)

  #Summarize posteriors for lambda
  lambda<-round(unlist(out$mean$r),2)
  lambda<-as.vector(lambda)
  lambda.lower<-round(unlist(out$q2.5$r),2)
  lambda.lower<-as.vector(lambda.lower)
  lambda.upper<-round(unlist(out$q97.5$r),2)
  lambda.upper<-as.vector(lambda.upper)

  #Summarize posteriors for stable equilibrium
  k.mean<-round(unlist(out$mean$K))
  k.mean<-as.vector(k.mean)
  k.lower<-round(unlist(out$q2.5$K))
  k.lower<-as.vector(k.lower)
  k.upper<-round(unlist(out$q97.5$K))
  k.upper<-as.vector(k.upper)

  #Use data.frame to package N and q results and save to working directory
  NQresults<-data.frame(WetlandPond=pond.name,Year=yrlab,NLower95=N.lower,
                        Pop_estimate=N.total,NUpper95=N.upper,QLower95=q.lower,
                        Detection_estimate=q.mean,QUpper95=q.upper)
  NQresults<-NQresults[with(NQresults,order(NQresults$WetlandPond,
                                            NQresults$Year)),
  ]

  #Recreate Wetland pond labels
  pond.name2<-as.character(unique(unlist(sort(count$pond_name))))

  #Use data.frame to package lambda and carrying capacity results and
  #save to working directory
  GRresults<-data.frame(WetlandPond=pond.name2,RLower95=lambda.lower,
                        Inst_Pop_G_Rate=lambda,RUpper95=lambda.upper,
                        KLower95=k.lower,Stable_Equilibrium=k.mean,
                        KUpper95=k.upper)


  if(species=="YCHUB"){

    #Capture and Write results to working directory (R Data Files)
    write.csv(NQresults,"YaquiChubPondAbundanceDetection.csv",row.names=F)
    write.csv(GRresults,"YaquiChubPondPopGrowthRate.csv",row.names=F)

    plot<-ggplot(NQresults,aes(Year,Pop_estimate,colour=factor(WetlandPond)))+
      geom_point(size=4)+
      geom_errorbar(aes(ymin=NLower95,ymax=NUpper95),width=.2,size=1.5)+
      facet_wrap(~WetlandPond,ncol=7)+
      guides(colour="none")+
      theme_bw()+
      xlab("Year")+
      ylab("Abundance")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=16),
            strip.text.x=element_text(size=12))
    print(plot)
    ggsave("YaquiChubWetlandPondAbundanceFigure.tiff",plot=plot,
           width=15,height=10,dpi=300)
  }else if(species=="BSHINER"){
    #Capture and Write results to working directory (R Data Files)
    write.csv(NQresults,"BeautifulShinerPondAbundanceDetection.csv",row.names=F)
    write.csv(GRresults,"BeautifulShinerPondPopGrowthRate.csv",row.names=F)

    plot<-ggplot(NQresults,aes(Year,Pop_estimate,colour=factor(WetlandPond)))+
      geom_point(size=4)+
      geom_errorbar(aes(ymin=NLower95,ymax=NUpper95),width=.2,size=1.5)+
      facet_wrap(~WetlandPond,ncol=7)+
      guides(colour="none")+
      theme_bw()+
      xlab("Year")+
      ylab("Abundance")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=16),
            strip.text.x=element_text(size=12))
    print(plot)
    ggsave("BeautifulShinerWetlandPondAbundanceFigure.tiff",plot=plot,
           width=15,height=10,dpi=300)
  }

  print("Bayesian Population Model Complete",quote=FALSE)
  print("Executing JAGS model to assess relationships with select habitat variables. This may take several minutes to hours.",quote=FALSE)

  #Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Habitat.Model.txt"
  cat("
      model{
      phi~dunif(0,100)

      #Abundance parameters
      beta.depth~dnorm(0,0.01)
      beta.oxygen~dnorm(0,0.01)
      beta.cond~dnorm(0,0.01)
      beta.ntu~dnorm(0,0.01)
      #beta.algal~dnorm(0,0.01)

      #Detection parameters
      alpha~dnorm(0,0.01)
      alpha.depth~dnorm(0,0.01)
      alpha.temp~dnorm(0,0.01)
      alpha.veg~dnorm(0,0.01)

      for(k in 1:npond){
      beta[k]~dnorm(mu_pond,tau_pond) #Random pond effect
      }

      mu_pond~dnorm(0,0.01)
      tau_pond~dgamma(0.01,0.01)
      sd_pond<-1/sqrt(tau_pond)

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)
      for(k in 1:npond){
      for(t in 1:nyear){
      N[i,k,t]~dpois(lambda[i,k,t])
      lambda[i,k,t]<-muL[i,k,t]*eta[i]
      log(muL[i,k,t])<-beta[k] + beta.depth*wdepth[i,k,t] +
      beta.oxygen*doxygen[i,k,t] + beta.cond*wcond[i,k,t] +
      beta.ntu*ntu[i,k,t] 
      }}}

      for(i in 1:nsite){
      for(j in 1:nrep){
      for(k in 1:npond){
      for(t in 1:nyear){
      y[i,j,k,t]~dbin(q[i,j,k,t],N[i,k,t])

      #Detection probabilities
      logit(q[i,j,k,t])<-alpha + alpha.depth*wdepth[i,k,t] +
      alpha.temp*wtemp[i,k,t] + alpha.veg*veg[i,k,t]
      }}}}

      }",fill=TRUE,file=modelFilename)

  #Initial values
  Nst<-apply(y,c(1,3,4),sum,na.rm=T)+1
  inits=function()list(N=Nst)

  #Bundle data
  ndata=list(y=y,nsite=nsite,nrep=nday,npond=npond,nyear=nyear,
             wdepth=wdepth,veg=veg,wtemp=wtemp,wcond=wcond,doxygen=doxygen,
             ntu=ntu)

  #Parameters monitored
  params=c("phi","beta","beta.depth","beta.oxygen","beta.cond","beta.ntu",
           "alpha.depth","alpha.temp","alpha.veg","sd_pond")

  #MCMC settings
  nc=4; nt=1; nb=15000; ni=500000

  #Call JAGS
  out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,n.chains=nc,
             n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,n.cores=nc,DIC=TRUE)

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdat$pond_name)))),nyear)

  #Summarize posteriors for detection (alpha parameters)
  Det.depth<-round(unlist(out2$mean$alpha.depth),2)
  Det.depth<-as.vector(Det.depth)
  DDepth.lower<-round(unlist(out2$q2.5$alpha.depth),2)
  DDepth.lower<-as.vector(DDepth.lower)
  DDepth.upper<-round(unlist(out2$q97.5$alpha.depth),2)
  DDepth.upper<-as.vector(DDepth.upper)

  Det.temp<-round(unlist(out2$mean$alpha.temp),2)
  Det.temp<-as.vector(Det.temp)
  DTemp.lower<-round(unlist(out2$q2.5$alpha.temp),2)
  DTemp.lower<-as.vector(DTemp.lower)
  DTemp.upper<-round(unlist(out2$q97.5$alpha.temp),2)
  DTemp.upper<-as.vector(DTemp.upper)

  Det.veg<-round(unlist(out2$mean$alpha.veg),2)
  Det.veg<-as.vector(Det.veg)
  DVeg.lower<-round(unlist(out2$q2.5$alpha.veg),2)
  DVeg.lower<-as.vector(DVeg.lower)
  DVeg.upper<-round(unlist(out2$q97.5$alpha.veg),2)
  DVeg.upper<-as.vector(DVeg.upper)

  Det.mean<-c(Det.depth,Det.veg,Det.temp)
  Det.lower<-c(DDepth.lower,DVeg.lower,DTemp.lower)
  Det.upper<-c(DDepth.upper,DVeg.upper,DTemp.upper)
  Variable<-c("Net Depth","Vegetation","Water Temperature")
  Parameter<-rep(c("Detection"),3)

  res2.D<-data.frame(Parameter=Parameter,Variable=Variable,Lower=Det.lower,
                     Mean=Det.mean,Upper=Det.upper)

  #Summarize posteriors for abundance (beta parameters)
  N.depth<-round(unlist(out2$mean$beta.depth),2)
  N.depth<-as.vector(N.depth)
  NDepth.lower<-round(unlist(out2$q2.5$beta.depth),2)
  NDepth.lower<-as.vector(NDepth.lower)
  NDepth.upper<-round(unlist(out2$q97.5$beta.depth),2)
  NDepth.upper<-as.vector(NDepth.upper)

  N.oxygen<-round(unlist(out2$mean$beta.oxygen),2)
  N.oxygen<-as.vector(N.oxygen)
  Noxygen.lower<-round(unlist(out2$q2.5$beta.oxygen),2)
  Noxygen.lower<-as.vector(Noxygen.lower)
  Noxygen.upper<-round(unlist(out2$q97.5$beta.oxygen),2)
  Noxygen.upper<-as.vector(Noxygen.upper)

  N.cond<-round(unlist(out2$mean$beta.cond),2)
  N.cond<-as.vector(N.cond)
  Ncond.lower<-round(unlist(out2$q2.5$beta.cond),2)
  Ncond.lower<-as.vector(Ncond.lower)
  Ncond.upper<-round(unlist(out2$q97.5$beta.cond),2)
  Ncond.upper<-as.vector(Ncond.upper)

  N.ntu<-round(unlist(out2$mean$beta.ntu),2)
  N.ntu<-as.vector(N.ntu)
  Nntu.lower<-round(unlist(out2$q2.5$beta.ntu),2)
  Nntu.lower<-as.vector(Nntu.lower)
  Nntu.upper<-round(unlist(out2$q97.5$beta.ntu),2)
  Nntu.upper<-as.vector(Nntu.upper)

  #N.algal<-round(unlist(out2$mean$beta.algal),2)
  #N.algal<-as.vector(N.algal)
  #Nalgal.lower<-round(unlist(out2$q2.5$beta.algal),2)
  #Nalgal.lower<-as.vector(Nalgal.lower)
  #Nalgal.upper<-round(unlist(out2$q97.5$beta.algal),2)
  #Nalgal.upper<-as.vector(Nalgal.upper)

  N.mean<-c(N.depth,N.oxygen,N.cond,N.ntu)
  N.lower<-c(NDepth.lower,Noxygen.lower,Ncond.lower,Nntu.lower)
  N.upper<-c(NDepth.upper,Noxygen.upper,Ncond.upper,Nntu.upper)
  Variable<-c("Water Depth","Dissolved Oxygen","Water Conductivity",
              "NTU")
  Parameter<-rep(c("Abundance"),4)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  phi.mu<-round(unlist(out2$mean$phi),2)
  phi.mu<-as.vector(phi.mu)
  phi.lower<-round(unlist(out2$q2.5$phi),2)
  phi.lower<-as.vector(phi.lower)
  phi.upper<-round(unlist(out2$q97.5$phi),2)
  phi.upper<-as.vector(phi.upper)

  sdpond.mu<-round(unlist(out2$mean$sd_pond),2)
  sdpond.mu<-as.vector(sdpond.mu)
  sdpond.lower<-round(unlist(out2$q2.5$sd_pond),2)
  sdpond.lower<-as.vector(sdpond.lower)
  sdpond.upper<-round(unlist(out2$q97.5$sd_pond),2)
  sdpond.upper<-as.vector(sdpond.upper)

  N.lower=c(phi.lower,sdpond.lower)
  N.upper=c(phi.upper,sdpond.upper)
  N.mean<-c(phi.mu,sdpond.mu)

  Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
  Parameter<-rep(c("Random Effect"),2)
  res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N,res2.D,res2.R)


  if(species=="YCHUB"){

    #Capture and Write results to working directory (R Data Files)
    write.csv(res2,"YaquiChubSHabitatParameters.csv",row.names=F)

    plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
      geom_point(size=4)+
      geom_vline(aes(xintercept=0.0),color="black",size=2)+
      geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.3,size=1.5)+
      facet_wrap(~Parameter,ncol=2,scales="free_x")+
      guides(colour="none")+
      theme_bw()+
      xlab("Parameter Estimate")+
      ylab("")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=16),
            strip.text.x=element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(plot)
    ggsave("YaquiChubWetlandPondParameterFigure.tiff",plot=plot,
           width=9,height=7,dpi=300)
  }else if(species=="BSHINER"){
    #Capture and Write results to working directory (R Data Files)
    write.csv(res2,"BeautifulShinerHabitatParameters.csv",row.names=F)

    plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
      geom_point(size=4)+
      geom_vline(aes(xintercept=0.0),color="black",size=2)+
      geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.3,size=1.5)+
      facet_wrap(~Parameter,ncol=2,scales="free_x")+
      guides(colour="none")+
      theme_bw()+
      xlab("Parameter Estimate")+
      ylab("")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=16),
            strip.text.x=element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(plot)
    ggsave("BeautifulShinerWetlandPondParameterFigure.tiff",plot=plot,width=7,
           height=9,dpi=300)
  }
  print("Results for the Habitat Model are saved and stored in your working directory.",quote=FALSE)
  print("Executing JAGS model to assess relationships with select Management variables. This may take several minutes to hours.",quote=FALSE)

#Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Management.Model.txt"
  cat("
      model{
      phi~dunif(0,100)

      #Abundance parameters
      beta~dnorm(0,0.01)
      beta.ytoprem~dnorm(0,0.01)
      beta.ytopstock~dnorm(0,0.01)
      beta.ychubrem~dnorm(0,0.01)
      beta.ychubstock~dnorm(0,0.01)
      beta.bshinerem~dnorm(0,0.01)
      beta.bshinestock~dnorm(0,0.01)

      #Detection parameters
      for(k in 1:npond){
      alpha[k]~dnorm(mu.alpha,tau.alpha)
      }
      mu.alpha~dnorm(0,0.01)
      tau.alpha~dgamma(0.1,0.1)
      sd_pond<-1/tau.alpha

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)
      for(k in 1:npond){
      for(t in 1:nyear){
      N[i,k,t]~dpois(lambda[i,k,t])
      lambda[i,k,t]<-muL[i,k,t]*eta[i]
      log(muL[i,k,t])<-beta + beta.ytoprem*ytoprem[i,k,t] +
      beta.ytopstock*ytopstock[i,k,t] + beta.ychubrem*ychubrem[i,k,t] + beta.ychubstock*ychubstock[i,k,t] +
      beta.bshinerem*bshinerem[i,k,t] + beta.bshinestock*bshinestock[i,k,t] 
      }}}

      for(i in 1:nsite){
      for(j in 1:nrep){
      for(k in 1:npond){
      for(t in 1:nyear){
      y[i,j,k,t]~dbin(q[i,j,k,t],N[i,k,t])

      #Detection probabilities
      logit(q[i,j,k,t])<-alpha[k]   
      }}}}

      }",fill=TRUE,file=modelFilename)

  #Initial values
  Nst<-apply(y,c(1,3,4),sum,na.rm=T)+1
  inits=function()list(N=Nst)

  #Bundle data
  ndata=list(y=y,nsite=nsite,nrep=nday,npond=npond,nyear=nyear,
             ytoprem=ytoprem,ytopstock=ytopstock,ychubrem=ychubrem,ychubstock=ychubstock,
            bshinerem=bshinerem,bshinestock=bshinestock)

  #Parameters monitored
  params=c("phi","beta","beta.ytoprem","beta.ytopstock","beta.ychubrem","beta.ychubstock",
           "beta.bshinerem","beta.bshinestock","alpha","sd_pond")

  #MCMC settings
  nc=4; nt=1; nb=15000; ni=500000

  #Call JAGS
  out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,n.chains=nc,
             n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,n.cores=nc,DIC=TRUE)

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdat$pond_name)))),nyear)

  #Summarize posteriors for abundance (beta parameters)
  N.ytopstock<-round(unlist(out2$mean$beta.ytopstock),2)
  N.ytopstock<-as.vector(N.ytopstock)
  Nytopstock.lower<-round(unlist(out2$q2.5$beta.ytopstock),2)
  Nytopstock.lower<-as.vector(Nytopstock.lower)
  Nytopstock.upper<-round(unlist(out2$q97.5$beta.ytopstock),2)
  Nytopstock.upper<-as.vector(Nytopstock.upper)
  
  N.ychubstock<-round(unlist(out2$mean$beta.ychubstock),2)
  N.ychubstock<-as.vector(N.ychubstock)
  Nychubstock.lower<-round(unlist(out2$q2.5$beta.ychubstock),2)
  Nychubstock.lower<-as.vector(Nychubstock.lower)
  Nychubstock.upper<-round(unlist(out2$q97.5$beta.ychubstock),2)
  Nychubstock.upper<-as.vector(Nychubstock.upper)
  
  N.bshinestock<-round(unlist(out2$mean$beta.bshinestock),2)
  N.bshinestock<-as.vector(N.bshinestock)
  Nbshinestock.lower<-round(unlist(out2$q2.5$beta.bshinestock),2)
  Nbshinestock.lower<-as.vector(Nbshinestock.lower)
  Nbshinestock.upper<-round(unlist(out2$q97.5$beta.bshinestock),2)
  Nbshinestock.upper<-as.vector(Nbshinestock.upper)
  
  N.ytoprem<-round(unlist(out2$mean$beta.ytoprem),2)
  N.ytoprem<-as.vector(N.ytoprem)
  Nytoprem.lower<-round(unlist(out2$q2.5$beta.ytoprem),2)
  Nytoprem.lower<-as.vector(Nytoprem.lower)
  Nytoprem.upper<-round(unlist(out2$q97.5$beta.ytoprem),2)
  Nytoprem.upper<-as.vector(Nytoprem.upper)
  
  N.ychubrem<-round(unlist(out2$mean$beta.ychubrem),2)
  N.ychubrem<-as.vector(N.ychubrem)
  Nychubrem.lower<-round(unlist(out2$q2.5$beta.ychubrem),2)
  Nychubrem.lower<-as.vector(Nychubrem.lower)
  Nychubrem.upper<-round(unlist(out2$q97.5$beta.ychubrem),2)
  Nychubrem.upper<-as.vector(Nychubrem.upper)
  
  N.bshinerem<-round(unlist(out2$mean$beta.bshinerem),2)
  N.bshinerem<-as.vector(N.bshinerem)
  Nbshinerem.lower<-round(unlist(out2$q2.5$beta.bshinerem),2)
  Nbshinerem.lower<-as.vector(Nbshinerem.lower)
  Nbshinerem.upper<-round(unlist(out2$q97.5$beta.bshinerem),2)
  Nbshinerem.upper<-as.vector(Nbshinerem.upper)
  
  N.mean<-c(N.ytopstock,N.ytoprem,N.ychubstock,N.ychubrem,N.bshinestock,N.bshinerem)
  N.lower<-c(Nytopstock.lower,Nytoprem.lower,Nychubstock.lower,Nychubrem.lower,Nbshinestock.lower,Nbshinerem.lower)
  N.upper<-c(Nytopstock.upper,Nytoprem.upper,Nychubstock.upper,Nychubrem.upper,Nbshinestock.upper,Nbshinerem.upper)
  Variable<-c("Yaqui Topminnow Stocked","Yaqui Topminnow Removed","Yaqui Chub Stocked",
             "Yaqui Chub Removed","Beautiful Shiner Stocked","Beautiful Shiner Removed")
  Parameter<-rep(c("Abundance"),6)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  phi.mu<-round(unlist(out2$mean$phi),2)
  phi.mu<-as.vector(phi.mu)
  phi.lower<-round(unlist(out2$q2.5$phi),2)
  phi.lower<-as.vector(phi.lower)
  phi.upper<-round(unlist(out2$q97.5$phi),2)
  phi.upper<-as.vector(phi.upper)

  sdpond.mu<-round(unlist(out2$mean$sd_pond),2)
  sdpond.mu<-as.vector(sdpond.mu)
  sdpond.lower<-round(unlist(out2$q2.5$sd_pond),2)
  sdpond.lower<-as.vector(sdpond.lower)
  sdpond.upper<-round(unlist(out2$q97.5$sd_pond),2)
  sdpond.upper<-as.vector(sdpond.upper)

  N.lower=c(phi.lower,sdpond.lower)
  N.upper=c(phi.upper,sdpond.upper)
  N.mean<-c(phi.mu,sdpond.mu)

  Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
  Parameter<-rep(c("Random Effect"),2)
  res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N,res2.R)


  if(species=="YCHUB"){

    #Capture and Write results to working directory (R Data Files)
    write.csv(res2,"YaquiChubManagementParameters.csv",row.names=F)

    plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
      geom_point(size=4)+
      geom_vline(aes(xintercept=0.0),color="black",size=2)+
      geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.3,size=1.5)+
      facet_wrap(~Parameter,ncol=2,scales="free_x")+
      guides(colour="none")+
      theme_bw()+
      xlab("Parameter Estimate")+
      ylab("")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=16),
            strip.text.x=element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(plot)
    ggsave("YaquiChubWetlandManagementParameterFigure.tiff",plot=plot,
           width=9,height=7,dpi=300)
  }else if(species=="BSHINER"){
    #Capture and Write results to working directory (R Data Files)
    write.csv(res2,"BeautifulShinerManagementParameters.csv",row.names=F)

    plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
      geom_point(size=4)+
      geom_vline(aes(xintercept=0.0),color="black",size=2)+
      geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.3,size=1.5)+
      facet_wrap(~Parameter,ncol=2,scales="free_x")+
      guides(colour="none")+
      theme_bw()+
      xlab("Parameter Estimate")+
      ylab("")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=16),
            strip.text.x=element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(plot)
    ggsave("BeautifulShinerWetlandManagementParameterFigure.tiff",plot=plot,width=9,
           height=7,dpi=300)
  }
  print("Results of the Management Model are saved and stored in your working directory.",quote=FALSE)




  print("If needed please consult with the Regional Statistician (Dr. David R. Stewart) or the Regional Data Management Team once complete and if you have any concerns.")
}
