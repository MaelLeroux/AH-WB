########################################################
# Script Playback AH-WB
#
# Leroux et al., 2022
########################################################



#-------------------------------------------------
# House keeping
#-------------------------------------------------

rm(list=ls())
setwd("~/Documents/Documents - Mael’s MacBook Pro/AH-WB/R")

library(brms)
library(mice)
library(lme4)
library(nlme)
library(glmmTMB)
library(DHARMa)
library(car)
library(PerformanceAnalytics)
library(RColorBrewer)
cls <- brewer.pal(9,"Set1")
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


#-------------------------------------------------
# 1) Upload/explore data
#-------------------------------------------------

dat <- read.csv("ahwb-data-playback.csv")
head(dat)
str(dat)
dat$cdt <- factor(dat$cdt,levels=c("AH","WB","AH-WB"))  
dat$cdt.num <- as.numeric(dat$cdt)
dat$ID <- as.factor(dat$ID)

# reorder
dat <- dat[order(dat$ID,partial=dat$cdt),]
dat$duration_look_speak[dat$duration_look_speak==0] <- 1
dat$latency_first_look[dat$latency_first_look==0] <- 1
dat$nb_look[dat$nb_look==0] <- 1

#-------------------------------------------------
# 2) GLMMs Analyses 
#-------------------------------------------------

##############
## Duration ##
##############

m1 <- glmmTMB(duration_look_speak ~ relevel(cdt,ref="AH-WB") + (1|ID),data=dat,family=Gamma(link="log"))
summary(m1)
simulateResiduals(m1, plot = T)
m1bis <- glmmTMB(duration_look_speak ~ relevel(cdt,ref="WB") + (1|ID),data=dat,family=Gamma(link="log"))
summary(m1bis)

#############
## Latency ##
#############

m2 <- glmmTMB(latency_first_look ~ relevel(cdt,ref="AH-WB") + (1|ID),data=dat,family=Gamma(link="log"))
summary(m2)
simulateResiduals(m2, plot = T)
m2bis <- glmmTMB(latency_first_look ~ relevel(cdt,ref="WB") + (1|ID),data=dat,family=Gamma(link="log"))
summary(m2bis)

#####################
## Number of looks ##
#####################

m3 <- glmmTMB(nb.look ~ relevel(cdt,ref="AH-WB") + (1|ID),data=dat,family=nbinom1)
summary(m3)
simulateResiduals(m3, plot = T)
m3bis <- glmmTMB(nb.look ~ relevel(cdt,ref="WB") + (1|ID),data=dat,family=nbinom1)
summary(m3bis)

#-------------------------------------------------
# 3) Bayesian analyses with multiple imputation 
#-------------------------------------------------

##############
## Duration ##
##############

# create data frame with missing conditions:
tapply(dat$ID,dat[,c('ID','cdt')],length)

dim(dat)
head(dat)

table(dat$cdt,dat$ID)

IDmiss <- c("KT", "KU", "MS")
cdtmiss <- c("WB", "WB", "AH")

new.dat <- data.frame(ID=c( "KT", "KU", "MS"),
                      cdt= c("WB", "WB", "AH"),
                      combination=rep("no",3),
                      duration_look_speak=rep(NA,3))

new.dat <- rbind(dat[,c("ID","cdt","combination","duration_look_speak")],new.dat)
new.dat$combination <- as.factor(new.dat$combination)
str(new.dat)

# Imputation with mice (using multiple mean matching)
imp <- mice(new.dat, m = 50, print = FALSE) # impute 50 new data sets with default method (pmm: predictive mean matching)
plot(density(as.numeric(imp$imp$duration_look_speak[2,]))) # values of missing data created

# store/extract imputed data
plotdat <- complete(imp,1) # show the 4th complete data set (out of 50)
plot(imp)
str(imp)

impdat <- complete(imp, "long")
plotdat$mean <- tapply(impdat$duration_look_speak, impdat$.id,mean)
plotdat$sd <- tapply(impdat$duration_look_speak, impdat$.id,sd)
plotdat$cdt.num <- as.numeric(plotdat$cdt)

# Analyse multiple imputation using brm_multiple Bayesian model 
m1_imp <- brm_multiple(duration_look_speak ~ relevel(cdt,ref="AH-WB") + (1|ID),family=Gamma(link="log"), data = imp)
summary(m1_imp)
m1_impbis <- brm_multiple(duration_look_speak ~ relevel(cdt,ref="WB")  + (1|ID),family=Gamma(link="log"), data = imp)
summary(m1_impbis)

#############
## Latency ##
#############

# create data frame with missing conditions:
tapply(dat$ID,dat[,c('ID','cdt')],length)

dim(dat)
head(dat)

table(dat$cdt,dat$ID)

IDmiss <- c("KT", "KU", "MS")
cdtmiss <- c("WB", "WB", "AH")

new.dat <- data.frame(ID=c( "KT", "KU", "MS"),
                      cdt= c("WB", "WB", "AH"),
                      combination=rep("no",3),
                      latency_first_look=rep(NA,3))

new.dat <- rbind(dat[,c("ID","cdt","combination","latency_first_look")],new.dat)
new.dat$combination <- as.factor(new.dat$combination)
str(new.dat)

# Imputation with mice (using multiple mean matching)
imp <- mice(new.dat, m = 50, print = FALSE) # impute 50 new data sets with default method (pmm: predictive mean matching)
plot(density(as.numeric(imp$imp$latency_first_look[2,]))) # values of missing data created

# store/extract imputed data
plotdat <- complete(imp,1) # show the 4th complete data set (out of 50)
plot(imp)
str(imp)

impdat <- complete(imp, "long")
plotdat$mean <- tapply(impdat$latency_first_look, impdat$.id,mean)
plotdat$sd <- tapply(impdat$latency_first_look, impdat$.id,sd)
plotdat$cdt.num <- as.numeric(plotdat$cdt)

# Analyse multiple imputation using brm_multiple Bayesian model
m2_imp <- brm_multiple(latency_first_look ~ relevel(cdt,ref="AH-WB")  + (1|ID),family=Gamma(link="log"), data = imp)
summary(m2_imp)
m2_impbis <- brm_multiple(latency_first_look ~ relevel(cdt,ref="WB")  + (1|ID),family=Gamma(link="log"), data = imp)
summary(m2_impbis)

#####################
## Number of looks ##
#####################

# create data frame with missing conditions:
tapply(dat$ID,dat[,c('ID','cdt')],length)

dim(dat)
head(dat)

table(dat$cdt,dat$ID)

IDmiss <- c("KT", "KU", "MS")
cdtmiss <- c("WB", "WB", "AH")

new.dat <- data.frame(ID=c( "KT", "KU", "MS"),
                      cdt= c("WB", "WB", "AH"),
                      combination=rep("no",3),
                      nb.look=rep(NA,3))

new.dat <- rbind(dat[,c("ID","cdt","combination","nb.look")],new.dat)
new.dat$combination <- as.factor(new.dat$combination)
str(new.dat)

# Imputation with mice (using multiple mean matching)
imp <- mice(new.dat, m = 50, print = FALSE) # impute 50 new data sets with default method (pmm: predictive mean matching)
plot(density(as.numeric(imp$imp$nb.look[2,]))) # values of missing data created

# store/extract imputed data
plotdat <- complete(imp,1) # show the 4th complete data set (out of 50)
plot(imp)
str(imp)

impdat <- complete(imp, "long")
plotdat$mean <- tapply(impdat$nb.look, impdat$.id,mean)
plotdat$sd <- tapply(impdat$nb.look, impdat$.id,sd)
plotdat$cdt.num <- as.numeric(plotdat$cdt)

# Analyse multiple imputation using brm_multiple Bayesian model
m3_imp <- brm_multiple(nb.look ~ relevel(cdt,ref="AH-WB")  + (1|ID),family=negbinomial, data = imp)
summary(m3_imp)
m3_impbis <- brm_multiple(nb.look ~ relevel(cdt,ref="WB")  + (1|ID),family=negbinomial, data = imp)
summary(m3_impbis)
