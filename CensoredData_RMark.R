########################
# Code to run Cormack-Jolly-Seber models to estimate survival of post-metamorphic lowland 
# leopard frogs in southern Arizona based on 3 years of capture-recapture data.
#
# Models run using the RMark package, which calls program MARK (must be installed prior to running this code)
#
# Reference:
# Zylstra, E.R., D.E. Swann, and R.J. Steidl. Surface-water availability governs survival of 
# an amphibian in arid mountain streams. Freshwater Biology 64:164-174.
########################

#--- Load libraries -----------------------------------------------------------------#

#install.packages(c("RMark","plyr"))
library(RMark)
library(plyr) 

#--- Load data ----------------------------------------------------------------------#

#set working directory
setwd("C:/...")

rm(list=ls(all=TRUE))

intdata <- read.table("CensoredData_Covariates.csv",sep=',',header=TRUE)

#--- Format and standardize survey-level covariates ---------------------------------#

#In general, variables that start with "end" have been created to use as covariates for recapture (p) because they reflect 
#conditions associated with the end of each survival interval (i.e., recapture occasion)
#Variables that have "mid" are to be used as covariates for survival (phi) because they reflect mean conditions or conditions at 
#the midpoint of each survival interval

intdata$midyr <- as.factor(intdata$midyr)
intdata$midwateryr <- as.factor(intdata$midwateryr)
intdata$endyr <- as.factor(intdata$endyr)
intdata$endwateryr <- as.factor(intdata$endwateryr)

intdata$end.waterz <- (intdata$end.water - mean(intdata$end.water))/sd(intdata$end.water)
intdata$end.humidz <- (intdata$end.humid - mean(intdata$end.humid))/sd(intdata$end.humid)
intdata$end.windz <- (intdata$end.wind - mean(intdata$end.wind))/sd(intdata$end.wind)
intdata$end.tempminz <- (intdata$end.tempmin - mean(intdata$end.tempmin))/sd(intdata$end.tempmin)
intdata$end.perimz <- (intdata$end.perim - mean(intdata$end.perim))/sd(intdata$end.perim)

intdata$water.meanz <- (intdata$water.mean - mean(intdata$water.mean))/sd(intdata$water.mean)
intdata$tminz <- (intdata$tmin - mean(intdata$tmin))/sd(intdata$tmin)
intdata$dewz <- (intdata$dew - mean(intdata$dew))/sd(intdata$dew)
intdata$mid.perimz <- (intdata$mid.perim - mean(intdata$mid.perim))/sd(intdata$mid.perim)

#--- Format things for RMark --------------------------------------------------------#

#Input file contains captures from May 2013 through 8 Nov 2015 (248 capture occasions)

#Converting MARK input file, with group IDs (codes for each stream reach)
llf.inp <- convert.inp("CensoredData_InputFile.inp",group.df=data.frame(site=c('CH','MA','RN','WHL','WHM','WHU')),use.comments=TRUE)
summary(llf.inp)

#Make processed dataframe and design data:
  #Create vector with interval lengths (in months if I divide by 30)
  #We'll need every 6th value since it's a time*site dataframe
  intervals <- intdata$duration[seq(1,nrow(intdata),by=6)]/30 
llf.process <- process.data(data=llf.inp,model="CJS",begin.time=1,time.intervals=intervals,groups="site")
llf.ddl <- make.design.data(llf.process)

#making occasion-based covariates:
time <- cumsum(c(1,intervals))  #if we assume time=1 at first occasion, time at other 247 occasions

p.covar <- intdata[,c('group','endyr','endwateryr','end.humidz','end.windz','end.tempminz','end.waterz','end.perimz')]
p.covar$time <- rep(time[2:length(time)],each=6)
names(p.covar) <- c('group','yr','wateryr','humid','wind','temp.min','water','veg.perim','time')
llf.ddl$p <- merge_design.covariates(llf.ddl$p,p.covar,bygroup=TRUE,bytime=TRUE)

phi.covar <- intdata[,c('group','midyr','midwateryr','tminz','water.meanz','dewz','mid.perimz')]
phi.covar$time <- rep(time[1:(length(time)-1)],each=6)
names(phi.covar) <- c('group','yr','wateryr','t.min','water.mn','dew','veg.perim','time')
llf.ddl$Phi <- merge_design.covariates(llf.ddl$Phi,phi.covar,bygroup=TRUE,bytime=TRUE)

#Numerous ways of viewing/checking the data:
# str(llf.ddl)
# llf.ddl$Phi[1:20,]
# llf.ddl$p[1:20,]
# summary(llf.ddl$Phi)
# summary(llf.ddl$p)

#--- Correlations among covariates? -------------------------------------------------#
#Covariates for p include:
    #wateryr (2013-2015) -- as factor
    #water (mean pool volume with water, interpolated values for each day within reach)
    #humidity (field measurement averaged across sites on particular date: mean)
    #wind (field measurement averaged across sites on particular date: max)
    #temperature (field measurement averaged across sites on particular date: min)
    #perimeter vegetation (reach- and wateryr-specific, percent perimeter with vegetation) 

#Covariates for Phi include: 
    #wateryr (2013-2015) -- as factor
    #water (mean pool volume with water, interpolated values for each day within reach)
    #dewpoint (PRISM[north/south]: mean of daily values)
    #temperature (PRISM[north/south]: mean of daily minimums)
    #perimeter vegetation (reach- and wateryr-specific, percent perimeter with vegetation) 

#Check potential correlations among covariates in recapture model:
    round(cor(intdata[,c('end.water','end.humid','end.wind','end.tempmin','end.perim')]),2)

#Check potential correlations among covariates in survival model:
    round(cor(intdata[,c('water.mean','dew','tmin','mid.perim')]),2)

#--- Running CJS models -------------------------------------------------------------#

models.comb <- function(){

  p.full  <- list(formula= ~wateryr+water+I(water^2)+humid+wind+temp.min+veg.perim)
  p.null  <- list(formula= ~1) 

  Phi.full    <- list(formula= ~wateryr+water.mn+I(water.mn^2)+t.min+dew+veg.perim)
  Phi.null     <- list(formula= ~1)
  Phi.noyear  <- list(formula= ~water.mn+I(water.mn^2)+t.min+dew+veg.perim)

  cml <- create.model.list("CJS")
  results <- mark.wrapper(cml,data=llf.process,ddl=llf.ddl,output=F)
  return(results)
}
results <- models.comb()

results$model.table[,c('model','npar','AICc','DeltaAICc','weight')]
  #model 1 = phi(full)p(full)
  #model 2 = phi(full)p(.)
  #model 3 = phi(everything except water year)p(full)
  #model 4 = phi(everything except water year)p(.)
  #model 5 = phi(.)p(full)
  #model 6 = phi(.)p(.)

#--- Overall/Mean estimates ---------------------------------------------------------#

#Estimates of survival from null model: phi(.)p(full)
(monphi <- results[[5]]$results$real[1,1])                #monthly = 0.83
(monphi.se <- results[[5]]$results$real[1,2])             #se = 0.013
(annphi <- monphi^12)                                     #annual = 0.11
(annphi.se <- sqrt((12^2)*(monphi^(11*2))*(monphi.se^2))) #se = 0.020
cbind(annphi,annphi-1.96*annphi.se,annphi+1.96*annphi.se)

#Estimates of survival from null model: phi(.)p(.)
results[[6]]$results$real                                         #monthly = 0.83
results[[6]]$results$real[1,2]/results[[6]]$results$real[1,1]*100 #CV = 1.59%

#Estimates of recapture probability from null model: phi(full)p(.)
results[[2]]$results$real[nrow(results[[2]]$results$real),]  #0.26 (0.011)

#Summarizing a few results from model that doesn't include year in survival submodel
phi3 <- coef(m3)[1:6,]
phi3 <- as.matrix(phi3)
nphi3 <- nrow(phi3)
p3   <- coef(m3)[7:15,]
p3   <- as.matrix(p3)
np3  <- nrow(p3)

  #Mean survival when water >=XX%?
  waterq <- seq(20,100,by=1)
  ll <- length(waterq)
  waterqq <- (waterq-mean(intdata$water.mean))/sd(intdata$water.mean)
  X.qq <- matrix(c(rep(1,ll),waterqq,waterqq^2,rep(0,ll),rep(0,ll),rep(0,ll)),byrow=F,ncol=nphi3)
  predsqq.logit <- X.qq %*% phi3[,1]
  predsqq <- plogis(predsqq.logit)
  std.errorsqq <- sqrt(diag(X.qq %*% m3$results$beta.vcv[1:nphi3,1:nphi3] %*% t(X.qq)))
  lcl.logitqq <- predsqq.logit - 1.96*std.errorsqq
  ucl.logitqq <- predsqq.logit + 1.96*std.errorsqq 
  cbind(waterq,predsqq,plogis(lcl.logitqq),plogis(ucl.logitqq))

#--- Inferences based on full model -------------------------------------------------#

m <- results[[1]]
phi <- coef(m)[1:8,] 
phi <- as.matrix(phi)
nphi <- nrow(phi)
p   <- coef(m)[9:17,]
p   <- as.matrix(p)
np  <- nrow(p)

#z-tests for covariate effects:
covars <- m$results$beta
covars$z <- covars$estimate/covars$se
covars$p <- 2*(1-pnorm(abs(covars$z),0,1))
covars

#---- Effects of surface water availability on survival and recapture probability ---#

## Survival: Water availability by year
  tapply(intdata$water.mean,intdata$midwateryr,summary)
  #minimums: 2013 = 0.16, 2014 = 0.04, 2015 = 0.30

  water13 <- seq(min(phi.covar$water.mn[phi.covar$wateryr==2013]),
                 max(phi.covar$water.mn[phi.covar$wateryr==2013]),length=100)
  water13p <- (water13*sd(intdata$water.mean) + mean(intdata$water.mean))/100
  water14 <- seq(min(phi.covar$water.mn[phi.covar$wateryr==2014]),
                 max(phi.covar$water.mn[phi.covar$wateryr==2014]),length=100)
  water14p <- (water14*sd(intdata$water.mean) + mean(intdata$water.mean))/100
  water15 <- seq(min(phi.covar$water.mn[phi.covar$wateryr==2015]),
                 max(phi.covar$water.mn[phi.covar$wateryr==2015]),length=100)
  water15p <- (water15*sd(intdata$water.mean) + mean(intdata$water.mean))/100
  #set weather variables at their mean:
  X13 <- matrix(c(rep(c(1,0,0),each=100),water13,water13^2,rep(0,300)),byrow=F,ncol=nphi)
  X14 <- matrix(c(rep(c(1,1,0),each=100),water14,water14^2,rep(0,300)),byrow=F,ncol=nphi)
  X15 <- matrix(c(rep(c(1,0,1),each=100),water15,water15^2,rep(0,300)),byrow=F,ncol=nphi)
  X <- rbind(X13,X14,X15)
  preds.logit <- X %*% phi[,1]
  preds <- plogis(preds.logit)
  std.errors <- sqrt(diag(X %*% m$results$beta.vcv[1:nphi,1:nphi] %*% t(X)))
  lcl.logit <- preds.logit - 1.96*std.errors
  ucl.logit <- preds.logit + 1.96*std.errors

## Recapture: Water availability by year
  tapply(intdata$end.water,intdata$endwateryr,summary)
  #minimums: 2013 = 0.16, 2014 = 0.03, 2015 = 0.30

  waterr13 <- seq(min(p.covar$water[p.covar$wateryr==2013]),
                 max(p.covar$water[p.covar$wateryr==2013]),length=100)
  waterr13p <- (water13*sd(intdata$end.water) + mean(intdata$end.water))/100
  waterr14 <- seq(min(p.covar$water[p.covar$wateryr==2014]),
                 max(p.covar$water[p.covar$wateryr==2014]),length=100)
  waterr14p <- (water14*sd(intdata$end.water) + mean(intdata$end.water))/100
  waterr15 <- seq(min(p.covar$water[p.covar$wateryr==2015]),
                 max(p.covar$water[p.covar$wateryr==2015]),length=100)
  waterr15p <- (water15*sd(intdata$end.water) + mean(intdata$end.water))/100  
  
  #set other variables at their mean:
  Xr13 <- matrix(c(rep(c(1,0,0),each=100),waterr13,waterr13^2,rep(0,400)),byrow=F,ncol=np)
  Xr14 <- matrix(c(rep(c(1,1,0),each=100),waterr14,waterr14^2,rep(0,400)),byrow=F,ncol=np)
  Xr15 <- matrix(c(rep(c(1,0,1),each=100),waterr15,waterr15^2,rep(0,400)),byrow=F,ncol=np)
  Xr <- rbind(Xr13,Xr14,Xr15)
  predsr.logit <- Xr %*% p[,1]
  predsr <- plogis(predsr.logit)
  std.errorsr <- sqrt(diag(Xr %*% m$results$beta.vcv[(nphi+1):(nphi+np),(nphi+1):(nphi+np)] %*% t(Xr)))
  lclr.logit <- predsr.logit - 1.96*std.errorsr
  uclr.logit <- predsr.logit + 1.96*std.errorsr

## Stacked plots
  
  cexfix <- 0.9  #will apply for axes.  Setting mtext = 0.9 of this value.  Legends 0.90 of this value
  #jpeg('....jpg',width=80,height=140,units='mm',res=600)
  #tiff('....tif',width=80,height=140,units='mm',res=600)
  par(mfrow=c(2,1),mar=c(0.5,3.0,0.5,0.7)+0.1,oma=c(2.0,0,0,0),cex=cexfix)
  plot(preds[1:100]~water13p,type='l',lty=2,col='black',xaxt='n',yaxt='n',xlab='',ylab='',
       ylim=c(-0.04,1),xlim=c(-0.04,1),bty='n',xaxs="i",yaxs='i',lwd=1)  
    usr <- par('usr')  #these are plotting limits (incl extra bit)
    axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
    axis(1,at=seq(0,1,by=0.2),labels=F,tcl=-0.25,mgp=c(1.5,0.5,0))
    axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
    axis(2,at=seq(0,1,by=0.2),labels=c('0.0','0.2','0.4','0.6','0.8','1.0'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0),cex.axis=0.9)
    lines(preds[101:200]~water14p,lty=1,col='black',lwd=1)
    lines(preds[201:300]~water15p,lty=3,col='black',lwd=1)
    polygon(x=c(water13p,rev(water13p)),y=c(plogis(lcl.logit[1:100]),rev(plogis(ucl.logit[1:100]))),border=NA,col=rgb(0,0,0,0.1))
    polygon(x=c(water14p,rev(water14p)),y=c(plogis(lcl.logit[101:200]),rev(plogis(ucl.logit[101:200]))),border=NA,col=rgb(0,0,0,0.1))
    polygon(x=c(water15p,rev(water15p)),y=c(plogis(lcl.logit[201:300]),rev(plogis(ucl.logit[201:300]))),border=NA,col=rgb(0,0,0,0.1))
    mtext('Monthly survival (95% CI)',side=2,las=0,line=2.1,cex=cexfix)
    text(x=0,y=0.95*usr[4],adj=c(0,0),labels='(a)',cex=0.95)
    legend(x=0.63,y=0.25,c('2013','2014','2015'),lty=c(2,1,3),col='black',bty='n',cex=cexfix)
  plot(predsr[1:100]~waterr13p,type='l',lty=2,col='black',xaxt='n',yaxt='n',xlab='',ylab='',
       ylim=c(-0.04,1),xlim=c(-0.04,1),bty='n',xaxs="i",yaxs='i',lwd=1)  
    usr <- par('usr')  #these are plotting limits (incl extra bit)
    axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
    axis(1,at=seq(0,1,by=0.2),labels=c('0.0','0.2','0.4','0.6','0.8','1.0'),
         tcl=-0.25,mgp=c(1.5,0.4,0),cex.axis=0.9)
    axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
    axis(2,at=seq(0,1,by=0.2),labels=c('0.0','0.2','0.4','0.6','0.8','1.0'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0),cex.axis=0.9)
    lines(predsr[101:200]~waterr14p,lty=1,col='black',lwd=1)
    lines(predsr[201:300]~waterr15p,lty=3,col='black',lwd=1)
    polygon(x=c(waterr13p,rev(waterr13p)),y=c(plogis(lclr.logit[1:100]),rev(plogis(uclr.logit[1:100]))),border=NA,col=rgb(0,0,0,0.1))
    polygon(x=c(waterr14p,rev(waterr14p)),y=c(plogis(lclr.logit[101:200]),rev(plogis(uclr.logit[101:200]))),border=NA,col=rgb(0,0,0,0.1))
    polygon(x=c(waterr15p,rev(waterr15p)),y=c(plogis(lclr.logit[201:300]),rev(plogis(uclr.logit[201:300]))),border=NA,col=rgb(0,0,0,0.1))
    mtext('Recapture probability (95% CI)',side=2,las=0,line=2.1,cex=cexfix)
    mtext('Surface-water availability',side=1,line=1.5,cex=cexfix)
    text(x=0,y=0.95*usr[4],adj=c(0,0),labels='(b)',cex=0.95)
  #dev.off()

  #mean survival when water >=60%? (looking at plot, survival will only increase with water >60%)
    water60 <- (60-mean(intdata$water.mean))/sd(intdata$water.mean)
    #set weather variables at their mean:
    X13.60 <- matrix(c(1,0,0,water60,water60^2,0,0,0),byrow=F,ncol=nphi)
    X14.60 <- matrix(c(1,1,0,water60,water60^2,0,0,0),byrow=F,ncol=nphi)
    X15.60 <- matrix(c(1,0,1,water60,water60^2,0,0,0),byrow=F,ncol=nphi)
    X.60 <- rbind(X13.60,X14.60,X15.60)
    preds60.logit <- X.60 %*% phi[,1]
    preds60 <- plogis(preds60.logit)
    std.errors60 <- sqrt(diag(X.60 %*% m$results$beta.vcv[1:nphi,1:nphi] %*% t(X.60)))
    lcl.logit60 <- preds60.logit - 1.96*std.errors60
    ucl.logit60 <- preds60.logit + 1.96*std.errors60 
    cbind(preds60,plogis(lcl.logit60),plogis(ucl.logit60))
    
  #mean survival when pools full?
    water100 <- (100-mean(intdata$water.mean))/sd(intdata$water.mean)
    #set weather variables at their mean:
    X13.100 <- matrix(c(1,0,0,water100,water100^2,0,0,0),byrow=F,ncol=nphi)
    X14.100 <- matrix(c(1,1,0,water100,water100^2,0,0,0),byrow=F,ncol=nphi)
    X15.100 <- matrix(c(1,0,1,water100,water100^2,0,0,0),byrow=F,ncol=nphi)
    X.100 <- rbind(X13.100,X14.100,X15.100)
    preds100.logit <- X.100 %*% phi[,1]
    preds100 <- plogis(preds100.logit)
    std.errors100 <- sqrt(diag(X.100 %*% m$results$beta.vcv[1:nphi,1:nphi] %*% t(X.100)))
    lcl.logit100 <- preds100.logit - 1.96*std.errors100
    ucl.logit100 <- preds100.logit + 1.96*std.errors100 
    cbind(preds100,plogis(lcl.logit100),plogis(ucl.logit100)) 
  
  #mean survival when water >=50%? (looking at plot, survival will only increase with water >50%)
    water50 <- (50-mean(intdata$water.mean))/sd(intdata$water.mean)
    #set weather variables at their mean:
    X13.50 <- matrix(c(1,0,0,water50,water50^2,0,0,0),byrow=F,ncol=nphi)
    X14.50 <- matrix(c(1,1,0,water50,water50^2,0,0,0),byrow=F,ncol=nphi)
    X15.50 <- matrix(c(1,0,1,water50,water50^2,0,0,0),byrow=F,ncol=nphi)
    X.50 <- rbind(X13.50,X14.50,X15.50)
    preds50.logit <- X.50 %*% phi[,1]
    preds50 <- plogis(preds50.logit)
    std.errors50 <- sqrt(diag(X.50 %*% m$results$beta.vcv[1:nphi,1:nphi] %*% t(X.50)))
    lcl.logit50 <- preds50.logit - 1.96*std.errors50
    ucl.logit50 <- preds50.logit + 1.96*std.errors50 
    cbind(preds50,plogis(lcl.logit50),plogis(ucl.logit50))  
  
  #mean survival when water at yearly min?
    #set weather variables at their mean:
    X13.min <- matrix(c(1,0,0,water13[1],water13[1]^2,0,0,0),byrow=F,ncol=nphi)
    X14.min <- matrix(c(1,1,0,water14[1],water14[1]^2,0,0,0),byrow=F,ncol=nphi)
    X15.min <- matrix(c(1,0,1,water15[1],water15[1]^2,0,0,0),byrow=F,ncol=nphi)
    X.min <- rbind(X13.min,X14.min,X15.min)
    predsmin.logit <- X.min %*% phi[,1]
    predsmin <- plogis(predsmin.logit)
    std.errorsmin <- sqrt(diag(X.min %*% m$results$beta.vcv[1:nphi,1:nphi] %*% t(X.min)))
    lcl.logitmin <- predsmin.logit - 1.96*std.errorsmin
    ucl.logitmin <- predsmin.logit + 1.96*std.errorsmin 
    cbind(predsmin,plogis(lcl.logitmin),plogis(ucl.logitmin))  

#Clean up MARK files:
rm(results)
cleanup(ask=F)
list.files()





