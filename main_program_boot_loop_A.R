#Main program, this is the only thing you need to run. 

library(survey) 
library(sampling)
library(dplyr)
library(stats)
library(sps)
library(boot)
library(rsample)
library(MASS)

#This is the only parameter you need changed for every misclassification scenario as shown in the paper, p1, p2

sens_low = 0.9
sens_high = 0.5

#Change the location where you put these programs
source("/Users/christofferdharma/Dropbox/PhD Stuff/Thesis_Work/Obj1-Simulation/AJE Paper/Codes_To_Submit/simulation/simulation_new/00_functions_A.R")
source("/Users/christofferdharma/Dropbox/PhD Stuff/Thesis_Work/Obj1-Simulation/AJE Paper/Codes_To_Submit/simulation/simulation_new/00_data_generation_A.R")

Nsimr = 500
B = 200

estot = 5
estimate.boot <- replicate(estot, matrix(nrow = B, ncol = Nsimr), simplify = FALSE)
mean.res <- matrix(, nrow = Nsimr, ncol = estot)

start_time <- Sys.time()

for (j in 1:Nsimr) {
  set.seed(j)
  
  datp <- as.data.frame(data_FP[poisson_sampling_VL(N,pi_p),])
  sps = sps_repweights(datp$di_p, replicates = B, dist = \(x) rexp(x) - 1)
  colnames(sps) <- paste0("BSW", 1:ncol(sps))
  datp$survey = 0
  
  datp_bsw <- cbind(datp,sps)
  
  datpA <- datp_bsw %>% dplyr::filter(A == 1)
  datpAp <- datp_bsw %>% dplyr::filter(Ap == 1)
  
  datc <- as.data.frame(FP1[poisson_sampling_VL(nrow(FP1),pi_c),])
  datc$survey = 1

  datctemp <- datc
  spsc = sps_repweights(datctemp$di_c, replicates = B, dist = \(x) rexp(x) - 1)
  colnames(spsc) <- paste0("BSW", 1:ncol(spsc))
  datctemp_bsw <- cbind(datctemp,spsc)
  
  for (i in 1:B) {
    datc[[as.symbol(paste0('BSW', i))]] <- 1
  }
  
  bt_samples_sc <- bootstraps(datc, times = B)
  bt_samples_sp <- bootstraps(datp_bsw, times = B)
  
  x.sc <- lapply(1:B, function(i) {
    data<-analysis(bt_samples_sc$splits[[i]])
    return(data)
  })
  
  x.sp <- lapply(1:B, function(i) {
    data<-analysis(bt_samples_sp$splits[[i]])
    return(data)
  })
  
  datALP <- ALP_only(datp=datp_bsw,datc=datc,xvarp=A)
  datALPAp <- ALP_only(datp=datp_bsw,datc=datc,xvarp=Ap)
  
  #First, collect the point estimates
  mean.res[j,1] = mean(datc$y) #unweighted Sc
  mean.res[j,2] =  sum(datctemp_bsw$y*datctemp_bsw$di_c)/sum(datctemp_bsw$di_c)
  mean.res[j,3] = sum(datALP$y*datALP$ALP)/sum(datALP$ALP) #ALP_weight true
  mean.res[j,4] = sum(datALPAp$y*datALPAp$ALP)/sum(datALPAp$ALP) #ALP_weight fake

  ######
  datc1 <- datALPAp %>% dplyr::filter(survey == 1)
  datp_all <- datp_bsw 
  
  datp_all_mike <- Mikeformula(datc1=datc1,datp_all=datp_all)
  result.two.step <- Imputation(M=1,datp_all=datp_all_mike,datc=datc1)
  mean.res[j,5] = mean(result.two.step)

  for (i in 1:B) {
    datALPA <- ALP_only_boot(datp=datp_bsw,datc=x.sc[[i]],xvarp=A)
    estimate.boot[[3]][i,j] = sum(datALPA$y*datALPA$ALP_weight)/sum(datALPA$ALP_weight)
    
    datALPAp <- ALP_only_boot(datp=datp_bsw,datc=x.sc[[i]],xvarp=Ap)
    estimate.boot[[4]][i,j] = sum(datALPAp$y*datALPAp$ALP_weight)/sum(datALPAp$ALP_weight)
    
    datc1 <- datALPAp %>% dplyr::filter(survey == 1)

    datp_all_mike <- Mikeformula(datc1=datc1,datp_all=datp_bsw,datpboot=x.sp[[i]],boot="T")
    result.two.step <- Imputation(M=1,datp_all=datp_all_mike,datc=datc1,boot="T")
    estimate.boot[[5]][i,j] = mean(result.two.step)
  }
  
  estimate.boot[[1]][,j] = colMeans(unlist(sapply(x.sc, function(x) x[, "y"]))) #unweighted nonprobability
  estimate.boot[[2]][,j] = colSums(datctemp_bsw$y*datctemp_bsw[,(ncol(datctemp_bsw)-B+1):ncol(datctemp_bsw)])/colSums(datctemp_bsw[,(ncol(datctemp_bsw)-B+1):ncol(datctemp_bsw)]) #true weight Sc
}
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
# 7.491591 hours

estimate.boot[[1]]
estimate.boot[[2]]
estimate.boot[[3]]
estimate.boot[[4]]
estimate.boot[[5]]

mean(mean.res[,1]) 
mean(mean.res[,2]) 
mean(mean.res[,3]) 
mean(mean.res[,4]) 
mean(mean.res[,5]) 

#This is the summary you need, it gives you the mean, RB and MSE
summary_function(mean.res[,1])
summary_function(mean.res[,2])
summary_function(mean.res[,3])
summary_function(mean.res[,4])
summary_function(mean.res[,5])

#1. unweighted Sc - corresponds to 1
#2. TW - corresponds to 2
#3. ALP weight true - corrsponds to 3
#4. ALP weight fake - corresponds to 4
#5. TSM - corresponds to 5 with 1 Sc boot

# apply 1 = row, 2 = column
#This gives you the variance * 10^2
mean(apply(estimate.boot[[1]],2,var))*100
mean(apply(estimate.boot[[2]],2,var))*100
mean(apply(estimate.boot[[3]],2,var))*100
mean(apply(estimate.boot[[4]],2,var))*100
mean(apply(estimate.boot[[5]],2,var))*100

#These give you the variance ratio
mean(apply(estimate.boot[[1]],2,var))/summary_function(mean.res[,1])[3]
mean(apply(estimate.boot[[2]],2,var))/summary_function(mean.res[,2])[3]
mean(apply(estimate.boot[[3]],2,var))/summary_function(mean.res[,3])[3]
mean(apply(estimate.boot[[4]],2,var))/summary_function(mean.res[,4])[3]
mean(apply(estimate.boot[[5]],2,var))/summary_function(mean.res[,5])[3]

#These give you the coverage probability (CP)
coverage_prob(matrixinp = estimate.boot[[1]],targetp=targetp,Nsim=Nsimr)
coverage_prob(matrixinp = estimate.boot[[2]],targetp=targetp,Nsim=Nsimr)
coverage_prob(matrixinp = estimate.boot[[3]],targetp=targetp,Nsim=Nsimr)
coverage_prob(matrixinp = estimate.boot[[4]],targetp=targetp,Nsim=Nsimr)
coverage_prob(matrixinp = estimate.boot[[5]],targetp=targetp,Nsim=Nsimr)
