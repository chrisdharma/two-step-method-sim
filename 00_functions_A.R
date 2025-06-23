#These are all the functions used in the simulation, will be called in the main program

library(survey) 
library(sampling)
library(dplyr)
library(stats)
library(MASS)

#Poisson sampling from VL
poisson_sampling_VL <- function(N, pi_i) {
  idx <- rbinom(N, size = 1, prob = pi_i)
  ord <- 1:N
  ord[idx == 1]
}

#ALP only - paper 1.
ALP_only <- function(datp,datc,xvarp) {

  datcm <- datc
  datcm$pi_c <- NULL
  datcm$di_c <- NULL
  
  datcm$di_p = 1
  
  datp1 <- datp %>% dplyr::filter({{xvarp}} == 1)
  
  datpc <- rbind(datp1,datcm)

  svyds.wei  = svydesign(ids=~1, weight = ~ di_p, data = datpc) ### 
  Formula_fit.survey = as.formula("survey ~ x1 + x2 + x3 + x4") #For ALP we need all 4
  lgtreg.w = svyglm(Formula_fit.survey, family= quasibinomial, design = svyds.wei)    

  beta.w = summary(lgtreg.w)$coeff[,1]   
  datpc$ps.w = lgtreg.w$fitted.values 

  datpc1 <- datpc %>% dplyr::filter(survey == 1) %>%
  mutate(ALP_weight = (1-ps.w) /ps.w)
  
  return(datpc1)
}

ALP_only_boot <- function(datp,datc,xvarp) {
  
  datcm <- datc
  datcm$pi_c <- NULL
  datcm$di_c <- NULL
  
  datcm$di_p = 1
  
  datp1 <- datp %>% dplyr::filter({{xvarp}} == 1)
  
  datpc <- rbind(datp1,datcm)
  
  variable = paste0("datpc$BSW",i)
  svyds.wei  = svydesign(ids=~1, weight = ~ eval(parse(text=variable)), data = datpc) ### 
  Formula_fit.survey = as.formula("survey ~ x1 + x2 + x3 + x4") 
  lgtreg.w = svyglm(Formula_fit.survey, family= quasibinomial, design = svyds.wei)    
  
  beta.w = summary(lgtreg.w)$coeff[,1]   
  datpc$ps.w = lgtreg.w$fitted.values 
  
  datpc1 <- datpc %>% dplyr::filter(survey == 1) %>%
    mutate(ALP_weight = (1-ps.w) /ps.w)
  
  return(datpc1)
}


#Simulate ALP_only
# simulation_ALPonly <- function(Nsim=Nsimr,xvarp) {
#   mean.res <- matrix(, nrow = Nsim, ncol = 2)
#   for (i in 1:Nsim) {
#     set.seed(i)
#     dat<-Poisson_sampling(pi_p=pi_p,pi_c=pi_c,popp=data_FP,popc=E1)
#     
#     datp<-dat[[1]]
#     datc<-dat[[2]]
#     datc$di = 1
#     
#     datp1 <- datp %>% dplyr::filter({{xvarp}} == 1)
#     
#     datALP <- ALP_only(datp=datp1,datc=datc)
#     mean.res[i,1] = weighted.mean(datALP$y,datALP$ALP_weight)
#     mean.res[i,2] = mean(datALP$y)
#   }
#   return(mean.res)
# }

Mikeformula <- function(datc1,datp_all,datpboot,boot="F") {

  #generate Atilde
  #Left hand side of formula
  m_sc = glm(Ap ~ x1 + x2 + x3 + x4, family = quasibinomial(link = logit), data = datc1)  ### Ap ~ x1 + x2 + x3 + x4 in sc
  #Apply this to Sp
  p_sc = expit(coef(m_sc)[1] + coef(m_sc)[2]*datp_all$x1 + coef(m_sc)[3]*datp_all$x2 + coef(m_sc)[4]*datp_all$x3 + coef(m_sc)[5]*datp_all$x4)
  
  #Right hand side of formula
  if (boot == "F") {m_sp = glm(Ap ~ x3 + x4, family = quasibinomial(link = logit), data = datp_all)}
  if (boot == "T") {m_sp = glm(Ap ~ x3 + x4, family = quasibinomial(link = logit), data = datpboot)}
  p_sp = expit(coef(m_sp)[1] + coef(m_sp)[2]*datp_all$x3 + coef(m_sp)[3]*datp_all$x4) 
  
  datp_all$pnewmsm <- ((1-p_sc)/p_sc)*(p_sp/(1-p_sp))

  return(datp_all)
}

# Imputation Method
Imputation <- function(M = 10,datp_all,datc,boot="F") {
imp <- c()
for (i in 1:M) {
  datp_all$newmsm<-ifelse(datp_all$Ap==1,1,rbinom(size=1,n=nrow(datp_all),prob=datp_all$pnewmsm))
  
  rm.datp <- "pnewmsm"
  datp_all_new <- datp_all %>% dplyr::filter(newmsm == 1) %>%
    dplyr::select(-any_of(rm.datp))
  
  rm.datc2 <- c("ps.w","ALP_weight")
  
  datc2 <- datc %>% dplyr::select(-any_of(rm.datc2)) %>%
    mutate(newmsm = 1)
  
  if (boot == "F") {newmsmdat <- ALP_only(datp=datp_all_new,datc=datc2,xvarp=newmsm)}
  if (boot == "T") {newmsmdat <- ALP_only_boot(datp=datp_all_new,datc=datc2,xvarp=newmsm)}
  
  datcallnewALP <- newmsmdat %>% dplyr::filter(survey == 1)
  outcometemp <- datcallnewALP %>% dplyr::select(y)
  outcome2 <- as.numeric(outcometemp[[1]])
  weight <- datcallnewALP %>% dplyr::select(c(ALP_weight))
  imp[i]<-sum(outcome2*weight[[1]])/sum(weight[[1]])
}
return(imp)
}


#Calculate the three evaluation criteria

relative_bias <- function(matrixinp,targetp,Nsim=Nsimr) {
  (sum((matrixinp-targetp) / targetp) * 100) / Nsim
}

emp_variance <- function(matrixinp,Nsim=Nsimr) {
  sum ((matrixinp - (sum(matrixinp)/Nsim))^2) / (Nsim-1)
}

mse <- function(matrixinp,targetp,Nsim=Nsimr) {
  sum((matrixinp-targetp)^2) / Nsim
}

coverage_prob <- function(matrixinp,targetp,Nsim) {
  capture_count <- c()
  for (i in 1:Nsimr) {
    # Calculate the 0.025 and 0.975 quantiles
    quantiles <- quantile(matrixinp[, i], probs = c(0.025, 0.975))
    lower_bound <- quantiles[1]
    upper_bound <- quantiles[2]
    capture_count[i] <- ifelse((targetp >= lower_bound) & (targetp <= upper_bound),1,0)
  }
  return(sum(capture_count)/Nsim)
}

summary_function <- function(inputmat) {
  relative_bias<-relative_bias(matrixinp=inputmat,targetp=targetp)
  emp_variance<-emp_variance(matrixinp=inputmat)
  mse<-mse(matrixinp=inputmat,targetp=targetp)
  mean_estimate<-mean(inputmat)
  output<-cbind(mean_estimate,relative_bias,emp_variance,mse)
  return(output)
}


#expit function
expit <- function(x) {
  exp(x)/(1+exp(x))
}

# cnst_sc_alp <- sapply(1:2, function(i) {
#   uniroot(f = function(const) {
#     pi_i <- exp(const + X%*%beta)
#     sum(pi_i) - nc[i]
#   }, interval = c(-100, 100))$root
# })

