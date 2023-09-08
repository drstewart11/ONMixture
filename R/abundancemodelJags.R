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
  #nsite = sites within a wetland pond
  #npond= # wetland ponds
  #nday = two survey days
  #nyear = the number of survey years
  nsite=max(countdata$site)
  npond=max(countdata$pname)
  nday=max(countdata$day)
  nyear=max(countdata$year)

  #Define and create 4 dimensional array
  #Read observed count data into 4-d array
  y=array(as.numeric(countdata$y),c(nsite,nday,npond,nyear))

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
      q[i,j,k,t]~dbeta(5,9) #Slightly informative prior
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
      #alpha_depth~dnorm(0,0.5)
      #alpha_temp~dnorm(0,0.5)
      alpha~dnorm(0,0.5)
      for(k in 1:npond){
      beta[k]~dnorm(0,0.01)
      r[k]~dunif(0,5)
      K[k]~dunif(50,5000)
      eta[k]~dgamma(phi,phi)
      }
      for(i in 1:nsite){
      for(k in 1:npond){
      N[i,k,1]~dpois(lambda[i,k])
      lambda[i,k]<-muL[i,k]*eta[k]
      log(muL[i,k])<-beta[k]
      for(t in 2:nyear){
      #mu[i,k,t-1]<-N[i,k,t-1]*exp(r[k]*(1-(N[i,k,t-1]/K[k])))
      mu[i,k,t-1]<-r[k]*N[i,k,t-1]*((K[k]-N[i,k,t-1])/K[k])
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
      q[i,j,k,t]~dbeta(2,4) #Slightly informative prior
      #logit(q[i,j,k,t])<-alpha + alpha_veg*veg[i,k,t]
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
  jags.data=list(y=y,nsite=nsite,nrep=nday,npond=npond,nyear=nyear)

  #Parameters monitored
  jags.params=c("phi","det.mean","K","r",
                "N.total")

  #MCMC settings
  nc=4; nt=1; nb=15000; ni=75000

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

    plot<-ggplot(NQresults,aes(as.factor(Year),Pop_estimate,colour=factor(WetlandPond)))+
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
           width=20,height=10,dpi=300)
  }else if(species=="BSHINER"){
    #Capture and Write results to working directory (R Data Files)
    write.csv(NQresults,"BeautifulShinerPondAbundanceDetection.csv",row.names=F)
    write.csv(GRresults,"BeautifulShinerPondPopGrowthRate.csv",row.names=F)

    plot<-ggplot(NQresults,aes(as.factor(Year),Pop_estimate,colour=factor(WetlandPond)))+
      geom_point(size=4)+
      geom_errorbar(aes(ymin=NLower95,ymax=NUpper95),width=.2,size=1.5)+
      facet_wrap(~WetlandPond,ncol=7,scales="free_y")+
      guides(colour="none")+
      theme_bw()+
      xlab("Year")+
      ylab("Abundance")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=16),
            strip.text.x=element_text(size=12))
    print(plot)
    ggsave("BeautifulShinerWetlandPondAbundanceFigure.tiff",plot=plot,
           width=20,height=10,dpi=300)
  }

  print("Bayesian Population Model Complete",quote=FALSE)





#############################################################################
#############################################################################
#############################################################################


  if(species=="YCHUB"){

  #Reorganize habitat data by Site, Wetland Pond, and Year
  newhab<-data.frame(pond_name=hab$pond_name,rep=hab$day,year=hab$year,site=hab$site,
                     pH=hab$pH,wtemp=hab$wtemp,
                     doxygen=hab$doxygen,wcond=hab$wcond,
                     veg=hab$veg,
                     wdepth=hab$wdepth,include=hab$include)
  #Filter
  habdat1<-newhab %>% filter(include == 1)

  habdat1$pname<-as.numeric(as.factor(habdat1$pond_name))
  habdat1$yr<-as.numeric(as.factor(habdat1$year))

  #Reorder habitat data
  habdata<-habdat1[order(habdat1$rep,habdat1$pname,habdat1$year,habdat1$site),]

  habdata<-data.table(pH=as.numeric(habdata$pH),
                      wtemp=as.numeric(habdata$wtemp),
                      doxygen=as.numeric(habdata$doxygen),
                      wcond=as.numeric(habdata$wcond),
                      veg=as.numeric(habdata$veg),
                      wdepth=as.numeric(habdata$wdepth),
                      yr = as.numeric(habdata$yr),
                      include=as.numeric(habdata$include))

  habdata<-setnafill(habdata, type = "locf")

  nsite=sum(habdata$include)/2


  #Scale habitat variables with mean 0 and input into a 3 dimensional array
  #pH

  pH=scale(habdata$pH)
  pH=array(habdata$pH,c(nsite,2))

  pH2=rowMeans(pH,na.rm=TRUE)

  #Water temperature
  wtemp=scale(habdata$wtemp)
  wtemp=array(wtemp,c(nsite,2))

  wtemp2=rowMeans(wtemp,na.rm=TRUE)


  #Dissolved oxygen
  doxygen=scale(habdata$doxygen)
  doxygen=array(doxygen,c(nsite,2))

  doxygen2=rowMeans(doxygen,na.rm=TRUE)


  #Water conductivity
  wcond=scale(habdata$wcond)
  wcond=array(wcond,c(nsite,2))

  wcond2=rowMeans(wcond,na.rm=TRUE)

  #Turbidity
  #ntu=scale(habdata$ntu)
  #ntu=array(ntu,c(nsite,2))

  #ntu2=rowMeans(ntu,na.rm=TRUE)

  #Algae concentration
  #algal=scale(habdata$algal)
  #algal=array(algal,c(nsite,npond,nyear))

  #Percent submergent aquatic vegetation
  veg=scale(habdata$veg)
  veg=array(veg,c(nsite,2))

  veg2=rowMeans(veg,na.rm=TRUE)

  #Water depth
  wdepth=scale(habdata$wdepth)
  wdepth=array(wdepth,c(nsite,2))

  wdepth2=rowMeans(wdepth, na.rm=TRUE)


  #Add pname field to count, mgmt, and hab .csv files
  #pname is a numeric wetland pond field
  #yr is a numeric year field
  count$pname<-as.numeric(as.factor(count$pond_name))
  count$yr<-as.numeric(as.factor(count$year))


  #Reorganize count data by Site, Wetland Pond, Year, and Species

    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,
                       site=count$site,y=count$YCHUB,include=count$include,
                       pond.name=count$pond_name)
    countdata<-newdat[order(newdat$day,newdat$pname,newdat$year,newdat$site),]

  #Filter
  countdat1<-countdata %>% filter(include == 1)
  countdat1$pname1<-as.numeric(as.factor(countdat1$pond.name))

  #Define 3-dimensional array dimensions
  #nsite = the twenty sites within a wetland pond
  #npond= 11 wetland ponds
  #nday = three survey days

  nrep=max(countdat1$day)

  #Define and create 4 dimensional array
  #Read observed count data into 4-d array
  y=array(as.numeric(countdat1$y),c(nsite,nrep))


  print("Executing JAGS model to assess relationships with select habitat variables. This may take several minutes to hours.",quote=FALSE)

  #Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Habitat.Model.txt"
  cat("
      model{
      phi~dunif(0.1,100)
      #Abundance parameters
      beta~dnorm(0,0.5)I(-2,2)
      beta.veg~dnorm(0,0.5)
      beta.depth~dnorm(0,0.5)I(0.0015,3)
      beta.wtemp~dnorm(0,0.5)

      #Detection parameters
      alpha~dnorm(0,0.01)
      alpha.cond~dnorm(0,0.5)
      alpha.veg~dnorm(0,0.5)

      for(k in 1:nponds){
      delta[k]~dnorm(mu_pond,tau_pond) #Random pond effect
      }
      mu_pond~dnorm(0,0.01)
      tau_pond~dgamma(0.01,0.01)
      sd_pond<-1/sqrt(tau_pond)

      #for(l in 1:nyears){
      #epsilon[l]~dnorm(mu_year,tau_year) #Random year effect
      #}
      #mu_year~dnorm(0,0.01)
      #tau_year~dgamma(0.01,0.01)
      #sd_year<-1/sqrt(tau_pond)

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)

      N[i]~dpois(lambda[i])
      lambda[i]<-muL[i]*eta[i]
      log(muL[i])<- beta + beta.veg*veg2[i] +
      beta.depth*wdepth2[i] + beta.wtemp*wtemp2[i] + delta[ponds[i]]
      }

      for(i in 1:nsite){
      for(j in 1:nrep){

      y[i,j]~dbin(q[i,j],N[i])

      #Detection probabilities
      logit(q[i,j])<-alpha + alpha.cond*wcond[i,j] + alpha.veg*veg[i,j]
      }}
      }",fill=TRUE,file=modelFilename)

  #Initial values
  Nst<-rowSums(y,na.rm=T)+1
  inits=function()list(N=Nst,phi=10)

  ponds<-dplyr::dense_rank(countdat1$pname)
  nponds<-max(ponds)
  years<-dplyr::dense_rank(countdat1$year)
  nyears<-max(years)

  #Bundle data
  ndata=list(y=y,nsite=nsite,nrep=nrep,
             wcond=wcond,veg=veg,
             wdepth2=wdepth2,veg2=veg2,wtemp2=wtemp2,
             ponds=ponds,nponds=nponds)

  #Parameters monitored
  params=c("phi","beta","beta.veg","beta.depth","beta.wtemp",
           "alpha.cond","alpha.veg","sd_pond")
  #MCMC settings
  nc=3; nt=1; nb=15000; ni=75000

  #Call JAGS
  out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,n.chains=nc,
             n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,n.cores=nc,DIC=TRUE)

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdat$pond_name)))))

  #Summarize posteriors for detection (alpha parameters)
  Det.cond<-round(unlist(out2$mean$alpha.cond),3)
  Det.cond<-as.vector(Det.cond)
  Dcond.lower<-round(unlist(out2$q2.5$alpha.cond),3)
  Dcond.lower<-as.vector(Dcond.lower)
  Dcond.upper<-round(unlist(out2$q97.5$alpha.cond),3)
  Dcond.upper<-as.vector(Dcond.upper)

  #Det.temp<-round(unlist(out2$mean$alpha.temp),2)
  #Det.temp<-as.vector(Det.temp)
  #DTemp.lower<-round(unlist(out2$q2.5$alpha.temp),2)
  #DTemp.lower<-as.vector(DTemp.lower)
  #DTemp.upper<-round(unlist(out2$q97.5$alpha.temp),2)
  #DTemp.upper<-as.vector(DTemp.upper)

  Det.veg<-round(unlist(out2$mean$alpha.veg),3)
  Det.veg<-as.vector(Det.veg)
  DVeg.lower<-round(unlist(out2$q2.5$alpha.veg),3)
  DVeg.lower<-as.vector(DVeg.lower)
  DVeg.upper<-round(unlist(out2$q97.5$alpha.veg),3)
  DVeg.upper<-as.vector(DVeg.upper)

  Det.mean<-c(Det.cond,Det.veg)
  Det.lower<-c(Dcond.lower,DVeg.lower)
  Det.upper<-c(Dcond.upper,DVeg.upper)
  Variable<-c("Water Conductivity","Vegetation")
  Parameter<-rep(c("Detection"),2)

  res2.D<-data.frame(Parameter=Parameter,Variable=Variable,Lower=Det.lower,
                     Mean=Det.mean,Upper=Det.upper)

  #Summarize posteriors for abundance (beta parameters)
  N.depth<-round(unlist(out2$mean$beta.depth),3)
  N.depth<-as.vector(N.depth)
  NDepth.lower<-round(unlist(out2$q2.5$beta.depth),3)
  NDepth.lower<-as.vector(NDepth.lower)
  NDepth.upper<-round(unlist(out2$q97.5$beta.depth),3)
  NDepth.upper<-as.vector(NDepth.upper)

  N.wtemp<-round(unlist(out2$mean$beta.wtemp),3)
  N.wtemp<-as.vector(N.wtemp)
  Nwtemp.lower<-round(unlist(out2$q2.5$beta.wtemp),3)
  Nwtemp.lower<-as.vector(Nwtemp.lower)
  Nwtemp.upper<-round(unlist(out2$q97.5$beta.wtemp),3)
  Nwtemp.upper<-as.vector(Nwtemp.upper)

  #N.cond<-round(unlist(out2$mean$beta.cond),2)
  #N.cond<-as.vector(N.cond)
  #Ncond.lower<-round(unlist(out2$q2.5$beta.cond),2)
  #Ncond.lower<-as.vector(Ncond.lower)
  #Ncond.upper<-round(unlist(out2$q97.5$beta.cond),2)
  #Ncond.upper<-as.vector(Ncond.upper)

  N.veg<-round(unlist(out2$mean$beta.veg),3)
  N.veg<-as.vector(N.veg)
  Nveg.lower<-round(unlist(out2$q2.5$beta.veg),3)
  Nveg.lower<-as.vector(Nveg.lower)
  Nveg.upper<-round(unlist(out2$q97.5$beta.veg),3)
  Nveg.upper<-as.vector(Nveg.upper)

  #N.algal<-round(unlist(out2$mean$beta.algal),2)
  #N.algal<-as.vector(N.algal)
  #Nalgal.lower<-round(unlist(out2$q2.5$beta.algal),2)
  #Nalgal.lower<-as.vector(Nalgal.lower)
  #Nalgal.upper<-round(unlist(out2$q97.5$beta.algal),2)
  #Nalgal.upper<-as.vector(Nalgal.upper)

  N.mean<-c(N.depth,N.wtemp,N.veg)
  N.lower<-c(NDepth.lower,Nwtemp.lower,Nveg.lower)
  N.upper<-c(NDepth.upper,Nwtemp.upper,Nveg.upper)
  Variable<-c("Water Depth","Water Temperature",
              "Vegetation")
  Parameter<-rep(c("Abundance"),3)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  phi.mu<-round(unlist(out2$mean$phi),3)
  phi.mu<-as.vector(phi.mu)
  phi.lower<-round(unlist(out2$q2.5$phi),3)
  phi.lower<-as.vector(phi.lower)
  phi.upper<-round(unlist(out2$q97.5$phi),3)
  phi.upper<-as.vector(phi.upper)

  sdpond.mu<-round(unlist(out2$mean$sd_pond),3)
  sdpond.mu<-as.vector(sdpond.mu)
  sdpond.lower<-round(unlist(out2$q2.5$sd_pond),3)
  sdpond.lower<-as.vector(sdpond.lower)
  sdpond.upper<-round(unlist(out2$q97.5$sd_pond),3)
  sdpond.upper<-as.vector(sdpond.upper)

  N.lower=c(phi.lower,sdpond.lower)
  N.upper=c(phi.upper,sdpond.upper)
  N.mean<-c(phi.mu,sdpond.mu)

  Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
  Parameter<-rep(c("Random Effect"),2)
  res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N,res2.D,res2.R)


  #Capture and Write results to working directory (R Data Files)
    write.csv(res2,"YaquiChubSHabitatModelParameters.csv",row.names=F)

    plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
      geom_point(size=4)+
      geom_vline(aes(xintercept=0.0),color="black",size=1)+
      geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,size=1)+
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
    ggsave("YaquiChubWetlandPondModelParameterFigure.tiff",plot=plot,
           width=12,height=7,dpi=300)

  print("Results for the Habitat Model Model are saved and stored in your working directory.",quote=FALSE)

  }else if(species=="BSHINER"){


    #Reorganize habitat data by Site, Wetland Pond, and Year
    newhab<-data.frame(pond_name=hab$pond_name,rep=hab$day,year=hab$year,site=hab$site,
                       pH=hab$pH,wtemp=hab$wtemp,
                       doxygen=hab$doxygen,wcond=hab$wcond,
                       veg=hab$veg,
                       wdepth=hab$wdepth,include=hab$include,onrefuge=hab$onrefuge)
    #Filter
    habdat1<-newhab %>% filter(include == 1 & onrefuge == 1)

    habdat1$pname<-as.numeric(as.factor(habdat1$pond_name))
    habdat1$yr<-as.numeric(as.factor(habdat1$year))

    #Reorder habitat data
    habdata<-habdat1[order(habdat1$rep,habdat1$pname,habdat1$year,habdat1$site),]

    habdata<-data.table(pH=as.numeric(habdata$pH),
                        wtemp=as.numeric(habdata$wtemp),
                        doxygen=as.numeric(habdata$doxygen),
                        wcond=as.numeric(habdata$wcond),
                        veg=as.numeric(habdata$veg),
                        wdepth=as.numeric(habdata$wdepth),
                        yr = as.numeric(habdata$yr),
                        include=as.numeric(habdata$include))

    habdata<-setnafill(habdata, type = "locf")

    nsite=sum(habdata$include)/2


    #Scale habitat variables with mean 0 and input into a 3 dimensional array
    #pH

    pH=scale(habdata$pH)
    pH=array(habdata$pH,c(nsite,2))

    pH2=rowMeans(pH,na.rm=TRUE)

    #Water temperature
    wtemp=scale(habdata$wtemp)
    wtemp=array(wtemp,c(nsite,2))

    wtemp2=rowMeans(wtemp,na.rm=TRUE)


    #Dissolved oxygen
    doxygen=scale(habdata$doxygen)
    doxygen=array(doxygen,c(nsite,2))

    doxygen2=rowMeans(doxygen,na.rm=TRUE)


    #Water conductivity
    wcond=scale(habdata$wcond)
    wcond=array(wcond,c(nsite,2))

    wcond2=rowMeans(wcond,na.rm=TRUE)

    #Turbidity
    #ntu=scale(habdata$ntu)
    #ntu=array(ntu,c(nsite,2))

    #ntu2=rowMeans(ntu,na.rm=TRUE)

    #Algae concentration
    #algal=scale(habdata$algal)
    #algal=array(algal,c(nsite,npond,nyear))

    #Percent submergent aquatic vegetation
    veg=scale(habdata$veg)
    veg=array(veg,c(nsite,2))

    veg2=rowMeans(veg,na.rm=TRUE)

    #Water depth
    wdepth=scale(habdata$wdepth)
    wdepth=array(wdepth,c(nsite,2))

    wdepth2=rowMeans(wdepth, na.rm=TRUE)


    #Add pname field to count, mgmt, and hab .csv files
    #pname is a numeric wetland pond field
    #yr is a numeric year field
    count$pname<-as.numeric(as.factor(count$pond_name))
    count$yr<-as.numeric(as.factor(count$year))


    #Reorganize count data by Site, Wetland Pond, Year, and Species
    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,
                         site=count$site,y=count$BSHINER,include=count$include,
                       onrefuge=count$onrefuge,pond.name=count$pond_name)
      countdata<-newdat[order(newdat$day,newdat$pname,newdat$year,newdat$site),]

    #Filter
    countdat1<-countdata %>% filter(include == 1 & onrefuge == 1)
    countdat1$pname1<-as.numeric(as.factor(countdat1$pond.name))

    #Define 3-dimensional array dimensions
    #nsite = the twenty sites within a wetland pond
    #npond= 11 wetland ponds
    #nday = three survey days

    nrep=max(countdat1$day)

    #Define and create 4 dimensional array
    #Read observed count data into 4-d array
    y=array(as.numeric(countdat1$y),c(nsite,nrep))

    print("Executing JAGS model to assess relationships with select habitat variables. This may take several minutes to hours.",quote=FALSE)

    #Define model file name and create JAGS model to assess relationships with
    #known habitat variables
    modelFilename="Habitat.Model.txt"
    cat("
      model{
      phi~dunif(0.1,100)
      #Abundance parameters
      beta~dnorm(0,0.5)I(-2,2)
      beta.veg~dnorm(0,0.5)
      beta.depth~dnorm(0,0.5)I(0.0015,3)
      beta.wtemp~dnorm(0,0.5)

      #Detection parameters
      alpha~dnorm(0,0.01)
      alpha.cond~dnorm(0,0.5)
      alpha.veg~dnorm(0,0.5)

      for(k in 1:nponds){
      delta[k]~dnorm(mu_pond,tau_pond) #Random pond effect
      }
      mu_pond~dnorm(0,0.01)
      tau_pond~dgamma(0.01,0.01)
      sd_pond<-1/sqrt(tau_pond)

      #for(l in 1:nyears){
      #epsilon[l]~dnorm(mu_year,tau_year) #Random year effect
      #}
      #mu_year~dnorm(0,0.01)
      #tau_year~dgamma(0.01,0.01)
      #sd_year<-1/sqrt(tau_pond)

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)

      N[i]~dpois(lambda[i])
      lambda[i]<-muL[i]*eta[i]
      log(muL[i])<- beta + beta.veg*veg2[i] +
      beta.depth*wdepth2[i] + beta.wtemp*wtemp2[i] + delta[ponds[i]]
      }

      for(i in 1:nsite){
      for(j in 1:nrep){

      y[i,j]~dbin(q[i,j],N[i])

      #Detection probabilities
      logit(q[i,j])<-alpha + alpha.cond*wcond[i,j] + alpha.veg*veg[i,j]
      }}
      }",fill=TRUE,file=modelFilename)

    #Initial values
    Nst<-rowSums(y,na.rm=T)+1
    inits=function()list(N=Nst,phi=10)

    ponds<-dplyr::dense_rank(countdat1$pname)
    nponds<-max(ponds)
    years<-dplyr::dense_rank(countdat1$year)
    nyears<-max(years)

    #Bundle data
    ndata=list(y=y,nsite=nsite,nrep=nrep,
               wcond=wcond,veg=veg,
               wdepth2=wdepth2,veg2=veg2,wtemp2=wtemp2,
               ponds=ponds,nponds=nponds)

    #Parameters monitored
    params=c("phi","beta","beta.veg","beta.depth","beta.wtemp",
             "alpha.cond","alpha.veg","sd_pond")
    #MCMC settings
    nc=3; nt=1; nb=15000; ni=75000

    #Call JAGS
    out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,
               n.chains=nc,n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,
               n.cores=nc,DIC=TRUE)

    #Create Wetland pond labels
    pond.name<-rep(as.character(unique(unlist(sort(countdat$pond_name)))))

    #Summarize posteriors for detection (alpha parameters)
    Det.cond<-round(unlist(out2$mean$alpha.cond),3)
    Det.cond<-as.vector(Det.cond)
    Dcond.lower<-round(unlist(out2$q2.5$alpha.cond),3)
    Dcond.lower<-as.vector(Dcond.lower)
    Dcond.upper<-round(unlist(out2$q97.5$alpha.cond),3)
    Dcond.upper<-as.vector(Dcond.upper)

    #Det.temp<-round(unlist(out2$mean$alpha.temp),2)
    #Det.temp<-as.vector(Det.temp)
    #DTemp.lower<-round(unlist(out2$q2.5$alpha.temp),2)
    #DTemp.lower<-as.vector(DTemp.lower)
    #DTemp.upper<-round(unlist(out2$q97.5$alpha.temp),2)
    #DTemp.upper<-as.vector(DTemp.upper)

    Det.veg<-round(unlist(out2$mean$alpha.veg),3)
    Det.veg<-as.vector(Det.veg)
    DVeg.lower<-round(unlist(out2$q2.5$alpha.veg),3)
    DVeg.lower<-as.vector(DVeg.lower)
    DVeg.upper<-round(unlist(out2$q97.5$alpha.veg),3)
    DVeg.upper<-as.vector(DVeg.upper)

    Det.mean<-c(Det.cond,Det.veg)
    Det.lower<-c(Dcond.lower,DVeg.lower)
    Det.upper<-c(Dcond.upper,DVeg.upper)
    Variable<-c("Water Conductivity","Vegetation")
    Parameter<-rep(c("Detection"),2)

    res2.D<-data.frame(Parameter=Parameter,Variable=Variable,Lower=Det.lower,
                       Mean=Det.mean,Upper=Det.upper)

    #Summarize posteriors for abundance (beta parameters)
    N.depth<-round(unlist(out2$mean$beta.depth),3)
    N.depth<-as.vector(N.depth)
    NDepth.lower<-round(unlist(out2$q2.5$beta.depth),3)
    NDepth.lower<-as.vector(NDepth.lower)
    NDepth.upper<-round(unlist(out2$q97.5$beta.depth),3)
    NDepth.upper<-as.vector(NDepth.upper)

    N.wtemp<-round(unlist(out2$mean$beta.wtemp),3)
    N.wtemp<-as.vector(N.wtemp)
    Nwtemp.lower<-round(unlist(out2$q2.5$beta.wtemp),3)
    Nwtemp.lower<-as.vector(Nwtemp.lower)
    Nwtemp.upper<-round(unlist(out2$q97.5$beta.wtemp),3)
    Nwtemp.upper<-as.vector(Nwtemp.upper)

    #N.cond<-round(unlist(out2$mean$beta.cond),2)
    #N.cond<-as.vector(N.cond)
    #Ncond.lower<-round(unlist(out2$q2.5$beta.cond),2)
    #Ncond.lower<-as.vector(Ncond.lower)
    #Ncond.upper<-round(unlist(out2$q97.5$beta.cond),2)
    #Ncond.upper<-as.vector(Ncond.upper)

    N.veg<-round(unlist(out2$mean$beta.veg),3)
    N.veg<-as.vector(N.veg)
    Nveg.lower<-round(unlist(out2$q2.5$beta.veg),3)
    Nveg.lower<-as.vector(Nveg.lower)
    Nveg.upper<-round(unlist(out2$q97.5$beta.veg),3)
    Nveg.upper<-as.vector(Nveg.upper)

    #N.algal<-round(unlist(out2$mean$beta.algal),2)
    #N.algal<-as.vector(N.algal)
    #Nalgal.lower<-round(unlist(out2$q2.5$beta.algal),2)
    #Nalgal.lower<-as.vector(Nalgal.lower)
    #Nalgal.upper<-round(unlist(out2$q97.5$beta.algal),2)
    #Nalgal.upper<-as.vector(Nalgal.upper)

    N.mean<-c(N.depth,N.wtemp,N.veg)
    N.lower<-c(NDepth.lower,Nwtemp.lower,Nveg.lower)
    N.upper<-c(NDepth.upper,Nwtemp.upper,Nveg.upper)
    Variable<-c("Water Depth","Water Temperature",
                "Vegetation")
    Parameter<-rep(c("Abundance"),3)

    res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                       Mean=N.mean,Upper=N.upper)

    phi.mu<-round(unlist(out2$mean$phi),3)
    phi.mu<-as.vector(phi.mu)
    phi.lower<-round(unlist(out2$q2.5$phi),3)
    phi.lower<-as.vector(phi.lower)
    phi.upper<-round(unlist(out2$q97.5$phi),3)
    phi.upper<-as.vector(phi.upper)

    sdpond.mu<-round(unlist(out2$mean$sd_pond),3)
    sdpond.mu<-as.vector(sdpond.mu)
    sdpond.lower<-round(unlist(out2$q2.5$sd_pond),3)
    sdpond.lower<-as.vector(sdpond.lower)
    sdpond.upper<-round(unlist(out2$q97.5$sd_pond),3)
    sdpond.upper<-as.vector(sdpond.upper)

    N.lower=c(phi.lower,sdpond.lower)
    N.upper=c(phi.upper,sdpond.upper)
    N.mean<-c(phi.mu,sdpond.mu)

    Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
    Parameter<-rep(c("Random Effect"),2)
    res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                       Mean=N.mean,Upper=N.upper)

    res2<-rbind(res2.N,res2.D,res2.R)


    #Capture and Write results to working directory (R Data Files)
      write.csv(res2,"BeautifulShinerHabitatModelParameters.csv",row.names=F)

      plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
        geom_point(size=4)+
        geom_vline(aes(xintercept=0.0),color="black",size=1)+
        geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,size=1)+
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
      ggsave("BeautifulShinerWetlandPondModelParameterFigure.tiff",plot=plot,width=12,
             height=7,dpi=300)

    print("Results for the Habitat Model Model are saved and stored in your working directory.",quote=FALSE)
  }

  print("Executing JAGS model to assess relationships with select Management variables. This may take several minutes to hours.",quote=FALSE)





################################################################################
################################################################################
################################################################################

  if(species=="YCHUB"){
    #Add pname field to count, mgmt, and hab .csv files
    #pname is a numeric wetland pond field
    #yr is a numeric year field
    count$pname<-as.numeric(as.factor(count$pond_name))
    count$yr<-as.numeric(as.factor(count$year))


    #Reorganize count data by Site, Wetland Pond, Year, and Species
    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,
                       site=count$site,y=count$BSHINER,include=count$include,
                       onrefuge=count$onrefuge,pond.name=count$pond_name)
    countdata<-newdat[order(newdat$day,newdat$pname,newdat$year,newdat$site),]

    #Filter
    countdat1<-countdata %>% filter(include == 1)
    countdat1$pname1<-as.numeric(as.factor(countdat1$pond.name))

    #Define 3-dimensional array dimensions
    #nsite = the twenty sites within a wetland pond
    #npond= 11 wetland ponds
    #nday = three survey days

    nrep=max(countdat1$day)

    #Define and create 4 dimensional array
    #Read observed count data into 4-d array
    y=array(as.numeric(countdat1$y),c(nsite,nrep))


  #Reorganize management activity data by Site, Wetland Pond, and Year
  newmgmt<-data.frame(pond_name=mgmt$pond_name,year=mgmt$year,ychubrem=mgmt$YCHUB_removed,
                      bshinerrem=mgmt$BSHINER_removed,ytoprem=mgmt$YTOP_removed,
                      ychubstock=mgmt$YCHUB_stocked,
                      bshinerstock=mgmt$BSHINER_stocked,
                      ytopstock=mgmt$YTOP_stocked,include=mgmt$include)

  #Filter
  mgmtdat1<-newmgmt %>% filter(include == 1)

  mgmtdat1$pname<-as.numeric(as.factor(mgmtdat1$pond_name))
  mgmtdat1$yr<-as.numeric(as.factor(mgmtdat1$year))

  #Reorder management activity data
  mgmtdata<-mgmtdat1[order(mgmtdat1$year,mgmtdat1$pname),]

  #Expand data.frame
  mgmtdat<-mgmtdata %>% group_by(pname,year,ychubrem,bshinerrem,ytoprem,
                                 ychubstock,bshinerstock,ytopstock) %>% expand(sites = 1:10)
  mgmtdat<-mgmtdat[order(mgmtdat$pname),]

  npond=max(unique(mgmtdat$pname))

  #Number of Yaqui chub removed
  ychubrem=mgmtdat$ychubrem


  #Number of beautiful shiner removed
  bshinerem=mgmtdat$bshinerrem
  #bshinerem=array(bshinerem,c(nsite,npond))

  #Number of Yaqui topminnow removed
  ytoprem=mgmtdat$ytoprem
  #ytoprem=array(ytoprem,c(nsite,npond))

  #Number of Yaqui chub stocked
  ychubstock=mgmtdat$ychubstock
  #ychubstock=array(ychubstock,c(nsite,npond))

  #Number of beautiful shiner stocked
  bshinerstock=mgmtdat$bshinerstock
  #bshinestock=array(bshinestock,c(nsite,npond))

  #Number of Yaqui topminnow stocked
  ytopstock=mgmtdat$ytopstock
  #ytopstock=array(ytopstock,c(nsite,npond))


  #Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Management.Model.txt"
  cat("
      model{
      phi~dunif(0.1,100)

      #Abundance parameters
      beta~dnorm(0,0.01)
      beta.ytopstock~dnorm(0,0.01)
      beta.bshinerstock~dnorm(0,0.01)

      #Detection parameters
      alpha~dnorm(0,0.5)

      for(k in 1:nponds){
      delta[k]~dnorm(mu.ponds,tau.ponds)
      }
      mu.ponds~dnorm(0,0.01)
      tau.ponds~dgamma(0.1,0.1)
      sd_pond<-1/tau.ponds

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)

      N[i]~dpois(lambda[i])
      lambda[i]<-muL[i]*eta[i]

      log(muL[i])<-beta + beta.ytopstock*ytopstock[i] +
      beta.bshinerstock*bshinerstock[i] +
      delta[ponds[i]]
      }

      for(i in 1:nsite){
      for(j in 1:nrep){
      y[i,j]~dbin(q[i,j],N[i])
      #Detection probabilities
      logit(q[i,j])<-alpha
      }}

      }",fill=TRUE,file=modelFilename)

  #Initial values
  #Nst<-apply(y,c(1,3,4),sum,na.rm=T)+1
  Nst<-rowSums(y,na.rm=T)+1
  inits=function()list(N=Nst,phi=10)

  #Bundle data
  ndata=list(y=y,nsite=nsite,nrep=nrep,nponds=nponds,ponds=ponds,
             ytopstock=ytopstock,
             bshinerstock=bshinerstock)

  #Parameters monitored
  params=c("phi","beta","beta.ytopstock",
           "beta.bshinerstock","alpha","sd_pond")

  #MCMC settings
  nc=4; nt=1; nb=10000; ni=25000

  #Call JAGS
  out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,n.chains=nc,
             n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,n.cores=nc,DIC=TRUE)

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdat1$pond_name)))))

  #Summarize posteriors for abundance (beta parameters)
  N.ytopstock<-round(unlist(out2$mean$beta.ytopstock),3)
  N.ytopstock<-as.vector(N.ytopstock)
  Nytopstock.lower<-round(unlist(out2$q2.5$beta.ytopstock),3)
  Nytopstock.lower<-as.vector(Nytopstock.lower)
  Nytopstock.upper<-round(unlist(out2$q97.5$beta.ytopstock),3)
  Nytopstock.upper<-as.vector(Nytopstock.upper)

  #N.ychubstock<-round(unlist(out2$mean$beta.ychubstock),2)
  #N.ychubstock<-as.vector(N.ychubstock)
  #Nychubstock.lower<-round(unlist(out2$q2.5$beta.ychubstock),2)
  #Nychubstock.lower<-as.vector(Nychubstock.lower)
  #Nychubstock.upper<-round(unlist(out2$q97.5$beta.ychubstock),2)
  #Nychubstock.upper<-as.vector(Nychubstock.upper)

  N.bshinerstock<-round(unlist(out2$mean$beta.bshinerstock),3)
  N.bshinerstock<-as.vector(N.bshinerstock)
  Nbshinerstock.lower<-round(unlist(out2$q2.5$beta.bshinerstock),3)
  Nbshinerstock.lower<-as.vector(Nbshinerstock.lower)
  Nbshinerstock.upper<-round(unlist(out2$q97.5$beta.bshinerstock),3)
  Nbshinerstock.upper<-as.vector(Nbshinerstock.upper)

  #N.mean<-c(N.ytopstock,N.ytoprem,N.ychubstock,N.ychubrem,N.bshinestock,N.bshinerem)
  #N.lower<-c(Nytopstock.lower,Nytoprem.lower,Nychubstock.lower,Nychubrem.lower,Nbshinestock.lower,Nbshinerem.lower)
  #N.upper<-c(Nytopstock.upper,Nytoprem.upper,Nychubstock.upper,Nychubrem.upper,Nbshinestock.upper,Nbshinerem.upper)

  N.mean<-c(N.ytopstock,N.bshinerstock)
  N.lower<-c(Nytopstock.lower,Nbshinerstock.lower)
  N.upper<-c(Nytopstock.upper,Nbshinerstock.upper)
  Variable<-c("Yaqui Topminnow Stocked",
              "Beautiful Shiner Stocked")
  Parameter<-rep(c("Abundance"),2)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  phi.mu<-round(unlist(out2$mean$phi),3)
  phi.mu<-as.vector(phi.mu)
  phi.lower<-round(unlist(out2$q2.5$phi),3)
  phi.lower<-as.vector(phi.lower)
  phi.upper<-round(unlist(out2$q97.5$phi),3)
  phi.upper<-as.vector(phi.upper)

  sdpond.mu<-round(unlist(out2$mean$sd_pond),3)
  sdpond.mu<-as.vector(sdpond.mu)
  sdpond.lower<-round(unlist(out2$q2.5$sd_pond),3)
  sdpond.lower<-as.vector(sdpond.lower)
  sdpond.upper<-round(unlist(out2$q97.5$sd_pond),3)
  sdpond.upper<-as.vector(sdpond.upper)

  N.lower=c(phi.lower,sdpond.lower)
  N.upper=c(phi.upper,sdpond.upper)
  N.mean<-c(phi.mu,sdpond.mu)

  Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
  Parameter<-rep(c("Random Effect"),2)
  res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N,res2.R)




    #Capture and Write results to working directory (R Data Files)
    write.csv(res2,"YaquiChubManagementStockingParameters.csv",row.names=F)

    plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
      geom_point(size=4)+
      geom_vline(aes(xintercept=0.0),color="black",size=1)+
      geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,size=1)+
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
    ggsave("YaquiChubWetlandManagementStockingParameterFigure.tiff",plot=plot,
           width=12,height=7,dpi=300)

  print("Results of the Management Model (Stocking) are saved and stored in your working directory.",quote=FALSE)


  #Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Management.Model.txt"
  cat("
      model{
      phi~dunif(0.1,100)

      #Abundance parameters
      beta~dnorm(0,0.01)
      beta.ytoprem~dnorm(0,0.01)
      beta.bshinerem~dnorm(0,0.01)


      #Detection parameters
      alpha~dnorm(0,0.5)

      for(k in 1:nponds){
      delta[k]~dnorm(mu.ponds,tau.ponds)
      }
      mu.ponds~dnorm(0,0.01)
      tau.ponds~dgamma(0.1,0.1)
      sd_pond<-1/tau.ponds

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)

      N[i]~dpois(lambda[i])
      lambda[i]<-muL[i]*eta[i]

      log(muL[i])<-beta + beta.ytoprem*ytoprem[i] +
      beta.bshinerem*bshinerem[i] +
      delta[ponds[i]]
      }

      for(i in 1:nsite){
      for(j in 1:nrep){
      y[i,j]~dbin(q[i,j],N[i])
      #Detection probabilities
      logit(q[i,j])<-alpha
      }}

      }",fill=TRUE,file=modelFilename)

  #Initial values
  #Nst<-apply(y,c(1,3,4),sum,na.rm=T)+1
  Nst<-rowSums(y,na.rm=T)+1
  inits=function()list(N=Nst,phi=10)

  #Bundle data
  ndata=list(y=y,nsite=nsite,nrep=nrep,nponds=nponds,ponds=ponds,
             ytoprem=ytoprem,
             bshinerem=bshinerem)

  #Parameters monitored
  params=c("phi","beta","beta.ytoprem",
           "beta.bshinerem","alpha","sd_pond")

  #MCMC settings
  nc=4; nt=1; nb=10000; ni=25000

  #Call JAGS
  out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,n.chains=nc,
             n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,n.cores=nc,DIC=TRUE)

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdat1$pond_name)))))

  N.ytoprem<-round(unlist(out2$mean$beta.ytoprem),3)
  N.ytoprem<-as.vector(N.ytoprem)
  Nytoprem.lower<-round(unlist(out2$q2.5$beta.ytoprem),3)
  Nytoprem.lower<-as.vector(Nytoprem.lower)
  Nytoprem.upper<-round(unlist(out2$q97.5$beta.ytoprem),3)
  Nytoprem.upper<-as.vector(Nytoprem.upper)

  #N.ychubrem<-round(unlist(out2$mean$beta.ychubrem),2)
  #N.ychubrem<-as.vector(N.ychubrem)
  #Nychubrem.lower<-round(unlist(out2$q2.5$beta.ychubrem),2)
  #Nychubrem.lower<-as.vector(Nychubrem.lower)
  #Nychubrem.upper<-round(unlist(out2$q97.5$beta.ychubrem),2)
  #Nychubrem.upper<-as.vector(Nychubrem.upper)

  N.bshinerem<-round(unlist(out2$mean$beta.bshinerem),3)
  N.bshinerem<-as.vector(N.bshinerem)
  Nbshinerem.lower<-round(unlist(out2$q2.5$beta.bshinerem),3)
  Nbshinerem.lower<-as.vector(Nbshinerem.lower)
  Nbshinerem.upper<-round(unlist(out2$q97.5$beta.bshinerem),3)
  Nbshinerem.upper<-as.vector(Nbshinerem.upper)

  #N.mean<-c(N.ytopstock,N.ytoprem,N.ychubstock,N.ychubrem,N.bshinestock,N.bshinerem)
  #N.lower<-c(Nytopstock.lower,Nytoprem.lower,Nychubstock.lower,Nychubrem.lower,Nbshinestock.lower,Nbshinerem.lower)
  #N.upper<-c(Nytopstock.upper,Nytoprem.upper,Nychubstock.upper,Nychubrem.upper,Nbshinestock.upper,Nbshinerem.upper)

  N.mean<-c(N.ytoprem,N.bshinerem)
  N.lower<-c(Nytoprem.lower,Nbshinerem.lower)
  N.upper<-c(Nytoprem.upper,Nbshinerem.upper)
  Variable<-c("Yaqui Topminnow Removed",
              "Beautiful Shiner Removed")
  Parameter<-rep(c("Abundance"),2)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  phi.mu<-round(unlist(out2$mean$phi),3)
  phi.mu<-as.vector(phi.mu)
  phi.lower<-round(unlist(out2$q2.5$phi),3)
  phi.lower<-as.vector(phi.lower)
  phi.upper<-round(unlist(out2$q97.5$phi),3)
  phi.upper<-as.vector(phi.upper)

  sdpond.mu<-round(unlist(out2$mean$sd_pond),3)
  sdpond.mu<-as.vector(sdpond.mu)
  sdpond.lower<-round(unlist(out2$q2.5$sd_pond),3)
  sdpond.lower<-as.vector(sdpond.lower)
  sdpond.upper<-round(unlist(out2$q97.5$sd_pond),3)
  sdpond.upper<-as.vector(sdpond.upper)

  N.lower=c(phi.lower,sdpond.lower)
  N.upper=c(phi.upper,sdpond.upper)
  N.mean<-c(phi.mu,sdpond.mu)

  Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
  Parameter<-rep(c("Random Effect"),2)
  res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N,res2.R)




    #Capture and Write results to working directory (R Data Files)
    write.csv(res2,"YaquiChubManagementRemovedParameters.csv",row.names=F)

    plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
      geom_point(size=4)+
      geom_vline(aes(xintercept=0.0),color="black",size=1)+
      geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,size=1)+
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
    ggsave("YaquiChubWetlandManagementRemovedParameterFigure.tiff",plot=plot,
           width=12,height=7,dpi=300)

  print("Results of the Management Model (Removed) are saved and stored in your working directory.",quote=FALSE)

  }else if(species=="BSHINER"){
    #Add pname field to count, mgmt, and hab .csv files
    #pname is a numeric wetland pond field
    #yr is a numeric year field
    count$pname<-as.numeric(as.factor(count$pond_name))
    count$yr<-as.numeric(as.factor(count$year))


    #Reorganize count data by Site, Wetland Pond, Year, and Species
    newdat<-data.frame(pname=count$pname,year=count$yr,day=count$day,
                       site=count$site,y=count$BSHINER,include=count$include,
                       onrefuge=count$onrefuge,pond.name=count$pond_name)
    countdata<-newdat[order(newdat$day,newdat$pname,newdat$year,newdat$site),]

    #Filter
    countdat1<-countdata %>% filter(include == 1 & onrefuge == 1)
    countdat1$pname1<-as.numeric(as.factor(countdat1$pond.name))

    #Define 3-dimensional array dimensions
    #nsite = the twenty sites within a wetland pond
    #npond= 11 wetland ponds
    #nday = three survey days

    nrep=max(countdat1$day)

    #Define and create 4 dimensional array
    #Read observed count data into 4-d array
    y=array(as.numeric(countdat1$y),c(nsite,nrep))
    #Reorganize management activity data by Site, Wetland Pond, and Year
    newmgmt<-data.frame(pond_name=mgmt$pond_name,year=mgmt$year,ychubrem=mgmt$YCHUB_removed,
                        bshinerrem=mgmt$BSHINER_removed,ytoprem=mgmt$YTOP_removed,
                        ychubstock=mgmt$YCHUB_stocked,
                        bshinerstock=mgmt$BSHINER_stocked,
                        ytopstock=mgmt$YTOP_stocked,include=mgmt$include,
                        onrefuge=mgmt$onrefuge)

    #Filter
    mgmtdat1<-newmgmt %>% filter(include == 1 & onrefuge == 1)

    mgmtdat1$pname<-as.numeric(as.factor(mgmtdat1$pond_name))
    mgmtdat1$yr<-as.numeric(as.factor(mgmtdat1$year))

    #Reorder management activity data
    mgmtdata<-mgmtdat1[order(mgmtdat1$year,mgmtdat1$pname),]

    #Expand data.frame
    mgmtdat<-mgmtdata %>% group_by(pname,year,ychubrem,bshinerrem,ytoprem,
                                   ychubstock,bshinerstock,ytopstock) %>% expand(sites = 1:10)
    mgmtdat<-mgmtdat[order(mgmtdat$pname),]

    npond=max(unique(mgmtdat$pname))

    #Number of Yaqui chub removed
    ychubrem=mgmtdat$ychubrem


    #Number of beautiful shiner removed
    bshinerem=mgmtdat$bshinerrem
    #bshinerem=array(bshinerem,c(nsite,npond))

    #Number of Yaqui topminnow removed
    ytoprem=mgmtdat$ytoprem
    #ytoprem=array(ytoprem,c(nsite,npond))

    #Number of Yaqui chub stocked
    ychubstock=mgmtdat$ychubstock
    #ychubstock=array(ychubstock,c(nsite,npond))

    #Number of beautiful shiner stocked
    bshinerstock=mgmtdat$bshinerstock
    #bshinestock=array(bshinestock,c(nsite,npond))

    #Number of Yaqui topminnow stocked
    ytopstock=mgmtdat$ytopstock
    #ytopstock=array(ytopstock,c(nsite,npond))


    #Define model file name and create JAGS model to assess relationships with
    #known habitat variables
    modelFilename="Management.Model.txt"
    cat("
      model{
      phi~dunif(0.1,100)

      #Abundance parameters
      beta~dnorm(0,0.01)
      beta.ytopstock~dnorm(0,0.01)
      beta.ychubstock~dnorm(0,0.01)

      #Detection parameters
      alpha~dnorm(0,0.5)

      for(k in 1:nponds){
      delta[k]~dnorm(mu.ponds,tau.ponds)
      }
      mu.ponds~dnorm(0,0.01)
      tau.ponds~dgamma(0.1,0.1)
      sd_pond<-1/tau.ponds

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)

      N[i]~dpois(lambda[i])
      lambda[i]<-muL[i]*eta[i]

      log(muL[i])<-beta + beta.ytopstock*ytopstock[i] +
      beta.ychubstock*ychubstock[i] +
      delta[ponds[i]]
      }

      for(i in 1:nsite){
      for(j in 1:nrep){
      y[i,j]~dbin(q[i,j],N[i])
      #Detection probabilities
      logit(q[i,j])<-alpha
      }}

      }",fill=TRUE,file=modelFilename)

    #Initial values
    #Nst<-apply(y,c(1,3,4),sum,na.rm=T)+1
    Nst<-rowSums(y,na.rm=T)+1
    inits=function()list(N=Nst,phi=10)

    #Bundle data
    ndata=list(y=y,nsite=nsite,nrep=nrep,nponds=nponds,ponds=ponds,
               ytopstock=ytopstock,ychubstock=ychubstock)

    #Parameters monitored
    params=c("phi","beta","beta.ytopstock","beta.ychubstock",
             "alpha","sd_pond")

    #MCMC settings
    nc=4; nt=1; nb=10000; ni=25000

    #Call JAGS
    out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,n.chains=nc,
               n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,n.cores=nc,DIC=TRUE)

    #Create Wetland pond labels
    pond.name<-rep(as.character(unique(unlist(sort(countdat1$pond_name)))))

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

    #N.bshinerstock<-round(unlist(out2$mean$beta.bshinerstock),2)
    #N.bshinerstock<-as.vector(N.bshinerstock)
    #Nbshinerstock.lower<-round(unlist(out2$q2.5$beta.bshinerstock),2)
    #Nbshinerstock.lower<-as.vector(Nbshinerstock.lower)
    #Nbshinerstock.upper<-round(unlist(out2$q97.5$beta.bshinerstock),2)
    #Nbshinerstock.upper<-as.vector(Nbshinerstock.upper)

    #N.mean<-c(N.ytopstock,N.ytoprem,N.ychubstock,N.ychubrem,N.bshinestock,N.bshinerem)
    #N.lower<-c(Nytopstock.lower,Nytoprem.lower,Nychubstock.lower,Nychubrem.lower,Nbshinestock.lower,Nbshinerem.lower)
    #N.upper<-c(Nytopstock.upper,Nytoprem.upper,Nychubstock.upper,Nychubrem.upper,Nbshinestock.upper,Nbshinerem.upper)

    N.mean<-c(N.ytopstock,N.ychubstock)
    N.lower<-c(Nytopstock.lower,Nychubstock.lower)
    N.upper<-c(Nytopstock.upper,Nychubstock.upper)
    Variable<-c("Yaqui Topminnow Stocked",
                "Yaqui Chub Stocked")
    Parameter<-rep(c("Abundance"),2)

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


    #Capture and Write results to working directory (R Data Files)
      write.csv(res2,"BeautifulShinerManagementStockingParameters.csv",row.names=F)

      plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
        geom_point(size=4)+
        geom_vline(aes(xintercept=0.0),color="black",size=1)+
        geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,size=1)+
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
      ggsave("BeautifulShinerWetlandManagementStockingParameterFigure.tiff",plot=plot,width=12,
             height=7,dpi=300)

    print("Results of the Management Model (Stocking) are saved and stored in your working directory.",quote=FALSE)


    #Define model file name and create JAGS model to assess relationships with
    #known habitat variables
    modelFilename="Management.Model.txt"
    cat("
      model{
      phi~dunif(0.1,100)

      #Abundance parameters
      beta~dnorm(0,0.01)
      beta.ytoprem~dnorm(0,0.01)
      beta.ychubrem~dnorm(0,0.01)


      #Detection parameters
      alpha~dnorm(0,0.5)

      for(k in 1:nponds){
      delta[k]~dnorm(mu.ponds,tau.ponds)
      }
      mu.ponds~dnorm(0,0.01)
      tau.ponds~dgamma(0.1,0.1)
      sd_pond<-1/tau.ponds

      for(i in 1:nsite){
      eta[i]~dgamma(phi,phi)

      N[i]~dpois(lambda[i])
      lambda[i]<-muL[i]*eta[i]

      log(muL[i])<-beta + beta.ytoprem*ytoprem[i] +
      beta.ychubrem*ychubrem[i] +
      delta[ponds[i]]
      }

      for(i in 1:nsite){
      for(j in 1:nrep){
      y[i,j]~dbin(q[i,j],N[i])
      #Detection probabilities
      logit(q[i,j])<-alpha
      }}

      }",fill=TRUE,file=modelFilename)

    #Initial values
    #Nst<-apply(y,c(1,3,4),sum,na.rm=T)+1
    Nst<-rowSums(y,na.rm=T)+1
    inits=function()list(N=Nst,phi=10)

    #Bundle data
    ndata=list(y=y,nsite=nsite,nrep=nrep,nponds=nponds,ponds=ponds,
               ytoprem=ytoprem,ychubrem=ychubrem)

    #Parameters monitored
    params=c("phi","beta","beta.ytoprem","beta.ychubrem",
             "alpha","sd_pond")

    #MCMC settings
    nc=4; nt=1; nb=15000; ni=75000

    #Call JAGS
    out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,n.chains=nc,
               n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,n.cores=nc,DIC=TRUE)

    #Create Wetland pond labels
    pond.name<-rep(as.character(unique(unlist(sort(countdat1$pond_name)))))

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

    #N.bshinerem<-round(unlist(out2$mean$beta.bshinerem),2)
    #N.bshinerem<-as.vector(N.bshinerem)
    #Nbshinerem.lower<-round(unlist(out2$q2.5$beta.bshinerem),2)
    #Nbshinerem.lower<-as.vector(Nbshinerem.lower)
    #Nbshinerem.upper<-round(unlist(out2$q97.5$beta.bshinerem),2)
    #Nbshinerem.upper<-as.vector(Nbshinerem.upper)

    #N.mean<-c(N.ytopstock,N.ytoprem,N.ychubstock,N.ychubrem,N.bshinestock,N.bshinerem)
    #N.lower<-c(Nytopstock.lower,Nytoprem.lower,Nychubstock.lower,Nychubrem.lower,Nbshinestock.lower,Nbshinerem.lower)
    #N.upper<-c(Nytopstock.upper,Nytoprem.upper,Nychubstock.upper,Nychubrem.upper,Nbshinestock.upper,Nbshinerem.upper)

    N.mean<-c(N.ytoprem,N.ychubrem)
    N.lower<-c(Nytoprem.lower,Nychubrem.lower)
    N.upper<-c(Nytoprem.upper,Nychubrem.upper)
    Variable<-c("Yaqui Topminnow Removed",
                "Yaqui Chub Removed")
    Parameter<-rep(c("Abundance"),2)

    res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                       Mean=N.mean,Upper=N.upper)

    phi.mu<-round(unlist(out2$mean$phi),3)
    phi.mu<-as.vector(phi.mu)
    phi.lower<-round(unlist(out2$q2.5$phi),3)
    phi.lower<-as.vector(phi.lower)
    phi.upper<-round(unlist(out2$q97.5$phi),3)
    phi.upper<-as.vector(phi.upper)

    sdpond.mu<-round(unlist(out2$mean$sd_pond),3)
    sdpond.mu<-as.vector(sdpond.mu)
    sdpond.lower<-round(unlist(out2$q2.5$sd_pond),3)
    sdpond.lower<-as.vector(sdpond.lower)
    sdpond.upper<-round(unlist(out2$q97.5$sd_pond),3)
    sdpond.upper<-as.vector(sdpond.upper)

    N.lower=c(phi.lower,sdpond.lower)
    N.upper=c(phi.upper,sdpond.upper)
    N.mean<-c(phi.mu,sdpond.mu)

    Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
    Parameter<-rep(c("Random Effect"),2)
    res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                       Mean=N.mean,Upper=N.upper)

    res2<-rbind(res2.N,res2.R)


      #Capture and Write results to working directory (R Data Files)
      write.csv(res2,"BeautifulShinerManagementRemovedParameters.csv",row.names=F)

      plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
        geom_point(size=4)+
        geom_vline(aes(xintercept=0.0),color="black",size=1)+
        geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,size=1)+
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
      ggsave("BeautifulShinerWetlandManagementRemovedParameterFigure.tiff",plot=plot,width=12,
             height=7,dpi=300)

    print("Results of the Management Model (Removed) are saved and stored in your working directory.",quote=FALSE)



  }


  print("If needed, please consult with the Regional Statistician (Dr. David R. Stewart) or the Regional Data Management Team once complete and if you have any concerns.")
}
