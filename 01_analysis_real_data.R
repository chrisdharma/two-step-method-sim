library(survey) 
library(tableone)
library(svrep)
library(dplyr)
library(tidyverse)

library(haven)
library(boot)
library(rsample)

setwd("/Users/christofferdharma/Documents/ThesisDat")
source("/Users/christofferdharma/Dropbox/PhD Stuff/Thesis_Work/Obj1-Simulation/AJE Paper/Codes_To_Submit/real_data/00_functions_real_data.R")

load("allbsw.Rdata")

M = 1
B = 200

#We want three estimates:
# 1. unweighted
# 2. ALP
# 3. two step

#Data prep:

allbsw <- allbsw %>% drop_na(age_grp,income,education,orient,msm,employ,rural,race,immigration,disclosure,survey,msm,WTS_M)
allbsw$existing_mh <- as.numeric(as.character(allbsw$existing_mh))
allbsw$depress_binary <- as.numeric(as.character(allbsw$depress_binary))
allbsw$mhsu_harmonized <- as.numeric(as.character(allbsw$mhsu_harmonized))
allbsw$consult_mh <- as.numeric(as.character(allbsw$consult_mh))
allbsw$poor_srmh <- as.numeric(as.character(allbsw$poor_srmh))
allbsw$life_stress_bin <- as.numeric(as.character(allbsw$life_stress_bin))
allbsw$local_com_belong_bin <- as.numeric(as.character(allbsw$local_com_belong_bin))
allbsw$loneliness_ucla_bin <- as.numeric(as.character(allbsw$loneliness_ucla_bin))
allbsw$community_connection_msm <- as.numeric(as.character(allbsw$community_connection_msm))
allbsw$WTS_M <- as.numeric(as.character(allbsw$WTS_M))

allmsm <- allbsw %>% filter(msm == 1)

#Subset the data
SNall <- allmsm %>% dplyr::filter(survey == 1)
cchsall <- allbsw %>% dplyr::filter(survey == 0) 

SNall1<-drop_outcome(dat=SNall,outcomevar=loneliness_ucla_bin)

x.sn<-getboot(dat=SNall1,B=B)

#Bootseq needs to be fixed just once
bootseq <- sample(seq(1:1000),B)

main_function(datc=SNall1,datp=cchsall,outcomevar=loneliness_ucla_bin)
main_function(datc=SNall1,datp=cchsall,outcomevar=community_connection_msm)

#Calculate from CCHS side

#1. CCHS weights
cchsall_msm <- cchsall %>% dplyr::filter(msm == 1)
main_function_cchs(outcomevar=local_com_belong_bin,datpmsm=cchsall_msm,datc=SNall,datp=cchsall)



# CCHS = all, B = 200
# > main_function(datc=SNall1,datp=cchsall,outcomevar=loneliness_ucla_bin)
# lower.2.5%       est upper.97.5%
#   uw       0.5020783 0.5125217   0.5231652
# alp      0.4558841 0.4744832   0.4929023
# twostep  0.4562196 0.4709097   0.4930615

#CCHS = all, B = 200
# > main_function(datc=SNall1,datp=cchsall,outcomevar=community_connection_msm)
# lower.2.5%       est upper.97.5%
#   uw       0.3119944 0.3232252   0.3350313
# alp      0.3367403 0.3501012   0.3711421
# twostep  0.3327314 0.3462935   0.3690493

#CCHS = all, B = 200,
# do we want to try M = 10?
# > main_function_cchs(outcomevar=local_com_belong_bin,datpmsm=cchsall_msm,datc=SNall,datp=cchsall)
# lower.2.5%       est upper.97.5%
#   cchs_uw  0.5619766 0.5994196   0.6408071
# twostep  0.5760722 0.6019891   0.6456141

# logistic_regression_fn <- function(data, indices) {
#   # Resample the data
#   sample_data <- data[indices, ]
#   
#   # Fit logistic regression model
#   model <- glm(Formula_fit.msm, data = sample_data, family = binomial)
#   
#   return(coef(model))  # Return model coefficients
# }
# 
# set.seed(456)  # For reproducibility
# B <- 100  # Number of bootstrap samples
# boot_results <- boot(data = cchsall, statistic = logistic_regression_fn, R = B)
# 
