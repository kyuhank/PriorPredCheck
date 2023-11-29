cat('\n\n Setup begins \n\n')  

#############
## library ##
#############

library(posterior)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(GGally)
library(tidyverse)
library(cmdstanr)
library(posterior)


###############################
### Load data and functions ###
###############################

source("../../Data_and_functions/Albacore_Data.R")
source("../../Data_and_functions/Global_functions.R")
source("../../Data_and_functions/ASPM_functions.R")

##################################
## read environmental variables ##
##################################

if (!interactive()) {
  
  nSimBoot <- as.integer(Sys.getenv("nSimBoot", 200))
  nSimSBC <- as.integer(Sys.getenv("nSimSBC", 200))
  
  ## mcmc sampling option ##
  nChains <- as.integer(Sys.getenv("nChains", 8))
  nWarmup <- as.integer(Sys.getenv("nWarmup", 1000))
  nSample <- as.integer(Sys.getenv("nSample", 1000))
  
  nWarmupDemon <- as.integer(Sys.getenv("nWarmupDemon", 5000))
  nSampleDemon <- as.integer(Sys.getenv("nSampleDemon", 5000))
  
}


cat('\n\n Setup finished \n\n') 

#############
#############
#############

cat('\n\n Run ASPM_PriorChecks.R \n\n')  
source("ASPM_PriorChecks.R")
cat('\n\n ASPM_PriorChecks.R finished \n\n')  

cat('\n\n Run ASPM_DemonAlbacore.R \n\n')  
source("ASPM_DemonAlbacore.R")
cat('\n\n ASPM_DemonAlbacore.R finished \n\n')  

cat('\n\n Run ASPM_SBC.R \n\n')  
source("ASPM_SBC.R")
cat('\n\n ASPM_SBC.R finished \n\n')  

cat('\n\n Run ASPM_Boot.R \n\n')  
source("ASPM_Boot.R")
cat('\n\n ASPM_Boot.R finished \n\n')  

