########################
# Code to run Cormack-Jolly-Seber models to estimate survival of post-metamorphic lowland 
# leopard frogs in southern Arizona based on 3 years of capture-recapture data.
#
# Models run using the RMark package, which calls program MARK (must be installed prior to running this code)
#
# This R code was used to estimate survival based on UNCENSORED capture-recapture data from 
# only those frogs with complete spot maps (information reported in the appendix of the following paper)
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

intdata <- read.table("UncensoredData_Covariates.csv",sep=',',header=TRUE)

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

#Input file contains captures from May 2013 through 8 Nov 2015 (259 capture occasions)

#Converting MARK input file, with group IDs (codes for each stream reach)
llf.inp <- convert.inp("UncensoredData_CompleteSpotMaps_InputFile.inp",group.df=data.frame(site=c('CH','MA','RN','WHL','WHM','WHU')),use.comments=TRUE)
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
results.completemaps.nocensor <- models.comb()

results.completemaps.nocensor$model.table[,c('model','npar','AICc','DeltaAICc','weight')]
  #model 1 = phi(full)p(full)
  #model 2 = phi(full)p(.)
  #model 3 = phi(everything except water year)p(full)
  #model 4 = phi(everything except water year)p(.)
  #model 5 = phi(.)p(full)
  #model 6 = phi(.)p(.)

#--- Overall/Mean estimates ---------------------------------------------------------#

#Estimates of survival from null model: phi(.)p(full)
(monphi <- results.completemaps.nocensor[[5]]$results$real[1,1])
(monphi.se <- results.completemaps.nocensor[[5]]$results$real[1,2])
(annphi <- monphi^12)
(annphi.se <- sqrt((12^2)*(monphi^(11*2))*(monphi.se^2)))
cbind(annphi,annphi-1.96*annphi.se,annphi+1.96*annphi.se) #(0.0589, 0.1155)

#Estimates of survival from null model: phi(.)p(.)
results.completemaps.nocensor[[6]]$results$real               #phi = 0.86
results.completemaps.nocensor[[6]]$results$real[1,2]/results.completemaps.nocensor[[6]]$results$real[1,1]*100 #CV = 1.25%

#Clean up MARK files:
rm(results.completemaps.nocensor)
cleanup(ask=F)
list.files()
