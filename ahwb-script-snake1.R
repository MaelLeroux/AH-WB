################################################################################
# Script snake presentation AH-WB
#
# Leroux et al., 2022                                                                                                                                                       
################################################################################


#-------------------------------------------------
# House keeping
#-------------------------------------------------

rm(list= ls())

setwd("~/Documents/Documents - Mael’s MacBook Pro/AH-WB/R")

library(reshape)
library(ggplot2)
library(pastecs)
library(psych)
library(car)
library(ggm)
library(polycor)
library(Hmisc)
library(boot)
library(QuantPsyc)
library(vegan)
library(irr)
library(DescTools)
library(nlme)
library(multcomp)
library(predictmeans)
library(ggthemes) 
library(gridExtra)
library(glmm)
library(trust)
library(digest)
library(Matrix)
library(mvtnorm)
library(lme4)
library(DHARMa)

#-------------------------------------------------
# 1) Upload/explore data
#-------------------------------------------------

dat<- read.csv("ahwb-data-snake.csv")

str(dat)

overview3$combination <- as.factor(dat$combination)
dat$recruit <- as.factor(dat$recruit)
dat$recruitable_id <- as.factor(dat$recruitable_id)

#-------------------------------------------------
# 2) GLMMs Analysis 
#-------------------------------------------------


####################################
## Number of individuals recuited ##
####################################

recruitedglmm<- glmer(nb_recruited ~ combination + (1|cdt) + (1|ID), data= dat, family = poisson)

Anova(recruitedglmm, type= "III")   
summary(recruitedglmm) 

simulateResiduals(recruitedglmm)

