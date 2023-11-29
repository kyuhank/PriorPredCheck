##############
## Preamble ##
##############

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(bakeR)
require(dplyr)

###############
### options ###
###############

nSimBoot=200
nSimSBC=200
nChains=8
nWarmup=1000
nSample=1000
nWarmupDemon=5000
nSampleDemon=5000

#######################
#### pre-processing ###
#######################

aspm_albacore <- list(list(pars = list("nSimBoot"=nSimBoot,
                                       "nSimSBC"=nSimSBC,
                                       "nChains"=nChains,
                                       "nWarmup"=nWarmup,
                                       "nSample"=nSample,
                                       "nWarmupDemon"=nWarmupDemon,
                                       "nSampleDemon"=nSampleDemon)))

names(aspm_albacore) <-"ASPM_albacore"

#######################
##### run the job #####
#######################

gateaux_job_runner(aspm_albacore,
                   server = "gateaux.io/api",
                   JWT=JWT, 
                   report_name = "Prior-checks-ASPM",
                   log_jobs = FALSE)
