
# ———————————————————————————————
# Prior predictive checks for the SSPM and ASPM by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# ———————————————————————————————


## South Atlantic Albacore Tuna  from Polacheck et al., 1993
## Fitting Surplus Production Models: Comparing Methods and Measuring Uncertainty
years=1967:1989

## Catch (MT)
Ct=c(15.9, 25.7, 28.5, 23.7, 25, 33.3, 28.2, 19.7, 17.5, 19.3, 21.6, 23.1, 22.5, 22.5, 23.6, 29.1, 14.4, 
     13.2, 28.4, 34.6, 37.5, 25.9, 25.3)

## cpue ( kg/100hooks)
It=c(61.89, 78.98, 55.59, 44.61, 56.89, 38.27, 33.84, 36.13, 41.95, 36.63, 36.33, 38.82, 34.32, 37.64, 34.01, 32.16, 26.88, 36.61, 30.07, 30.75, 23.36, 22.36, 21.91)

Et=Ct/It


ntimes=length(years)

names(Ct)<-years
names(It)<-years



######## biological info (for ASPM) from Punt et al., 1995 #########
##  Stock assessment and risk analysis for the South Atlantic population of albacore Thunnus alalunga using an age-structured production model ##

########################
###### Input pars ######
########################

a_mat=5
a_sel=3.5
nages=12

sig_sel=0.2
sig_mat=1e-10


Linf=124.74
kappa=0.2284
a0=-0.9892
lwa=1.3718*1e-5*1e-6
lwb=3.0973
sigR=sqrt(log(0.4^2 + 1)) ## CV=0.4
FemaleProp=0.5

## assumed median values of M and steepness parameters
M=0.3
h=0.85

tau2=log(0.1^2 + 1) ## CV=0.1


## R0
B0PriorFromPunt=c(80, 300)

