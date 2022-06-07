install.packages("jagsUI",repos = "http://cran.us.r-project.org")
install.packages("rjags",repos = "http://cran.us.r-project.org")
install.packages("ggplot2",repos = "http://cran.us.r-project.org")
install.packages("devtools",repos = "http://cran.us.r-project.org")
install.packages("data.table",repos = "http://cran.us.r-project.org")
install.packages("dplyr",repos = "http://cran.us.r-project.org")
install.packages("tidyr",repos = "http://cran.us.r-project.org")
install.packages("rstan",repos = "http://cran.us.r-project.org")

library(jagsUI)
library(rjags)
library(ggplot2)
library(devtools)
library(data.table)
#install_github("drstewart11/ONMixture")
#library(ONMixture)
library(dplyr)
library(tidyr)
library(rstan)

setwd("C:/Users/dstewart/Documents/R Data Files/Yaqui Fish")
countdat=read.table("RYaqui_pondsurvey_2021_2021.csv",header=T,sep=",",na.strings="")
mgmtdat=read.table("RYaqui_management_2021_2021.csv",header=T,sep=",",na.strings="")
habdat=read.table("RYaqui_habitat_2021_2021.csv",header=T,sep=",",na.strings="")

#Abundance, detection, trend, species:habitat relationships
#countmix(count=countdat,mgmt=mgmtdat,hab=habdat,species="BSHINER")

#Percent small: habitat relationships
#betamod(countdat,mgmtdat,habdat,species=c("BSHINER"))





