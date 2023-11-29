##############
## Preamble ##
##############

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(bakeR)
require(dplyr)

###############
### options ###
###############

nSimBoot=500
nSimSBC=1500
nChains=8
nWarmup=1000
nSample=1000
nWarmupDemon=5000
nSampleDemon=5000

#######################
#### pre-processing ###
#######################

sspm_albacore <- list(list(pars = list("nSimBoot"=nSimBoot,
                                       "nSimSBC"=nSimSBC,
                                       "nChains"=nChains,
                                       "nWarmup"=nWarmup,
                                       "nSample"=nSample,
                                       "nWarmupDemon"=nWarmupDemon,
                                       "nSampleDemon"=nSampleDemon)))

names(sspm_albacore) <-"SSPM_albacore"

#######################
##### run the job #####
#######################

gateaux_job_runner(sspm_albacore,
                   server = "gateaux.io/api",
                   JWT=JWT, 
                   report_name = "Prior-checks-SSPM",
                   log_jobs = FALSE)

