cat('\n\n Setup begins \n\n')  

#############
## library ##
#############

capture.output(if (!interactive()) {
  
  install.packages('../../cmdstanr-0.5.2.tar.gz', repos=NULL)
  source('../../installcmdstan.R')
  require(cmdstanr)
  install_cmdst(dir = '../../', cores = 8)
  cmdstanr::set_cmdstan_path("./cmdstan-2.29.2")
  
  LocalRun=0
  
}, file=nullfile())

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
source("../../Data_and_functions/SSPM_functions.R")


##################################
## read environmental variables ##
##################################

if (!interactive()) {
  
  nSimBoot <- as.integer(Sys.getenv("nSimBoot", 10))
  nSimSBC <- as.integer(Sys.getenv("nSimSBC", 10))
  
  ## mcmc sampling option ##
  nChains <- as.integer(Sys.getenv("nChains", 5))
  nWarmup <- as.integer(Sys.getenv("nWarmup", 500))
  nSample <- as.integer(Sys.getenv("nSample", 500))
  
  nWarmupDemon <- as.integer(Sys.getenv("nWarmupDemon", 1000))
  nSampleDemon <- as.integer(Sys.getenv("nSampleDemon", 1000))
  
  }


cat('\n\n Setup finished \n\n')  

#############
#############
#############

cat('\n\n Run SSPM_PriorChecks.R \n\n')  
source("SSPM_PriorChecks.R")
cat('\n\n SSPM_PriorChecks.R finished \n\n')  

cat('\n\n Run SSPM_DemonAlbacore.R \n\n')  
source("SSPM_DemonAlbacore.R")
cat('\n\n SSPM_DemonAlbacore.R finished \n\n')  

cat('\n\n Run SSPM_SBC.R \n\n')  
source("SSPM_SBC.R")
cat('\n\n SSPM_SBC.R finished \n\n')  

cat('\n\n Run SSPM_Boot.R \n\n')  
source("SSPM_Boot.R")
cat('\n\n SSPM_Boot.R finished \n\n')  

