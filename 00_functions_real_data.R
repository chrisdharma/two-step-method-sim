#These are all the functions used in the simulation, will be called in the main program

library(survey) 
library(sampling)
library(dplyr)
library(stats)
library(MASS)

#Formulas:
Formula_fit.survey = as.formula("survey ~ rural + age_grp + (age_grp**2) + income + education + race + employ + immigration + age_grp*employ") 
Formula_fit.msm = as.formula("msm ~ rural + age_grp + (age_grp**2) + income + education + race + employ + immigration + age_grp*employ") 
Formula_fit.disc = as.formula("disclosure ~ rural + age_grp + (age_grp**2) + income + education + race + employ + immigration + age_grp*employ") 

getboot <- function(dat,B) {
  
  bt_samples <- bootstraps(dat, times = B)
  
  x.spc <- lapply(1:B, function(i) {
    data<-analysis(bt_samples$splits[[i]])
    return(data)
  })
  
  return(x.spc)
}

drop_outcome <- function(dat,outcomevar) {  
  newdat <- dat %>% drop_na({{outcomevar}})
  return(newdat)
}

boot.sn_unadj <- function(outcomevar) {
  newdat <- SNall %>% drop_na({{outcomevar}})
  outcomeall <- newdat %>% dplyr::select({{outcomevar}})
  estimate<- prop.table(table(outcomeall[[1]]))[[2]]
  bs_estimates<-c()
  for (i in 1:B) {
    bt_samples_sn <- bootstraps(newdat, times = B)
    x.sn<-analysis(bt_samples_sn$splits[[i]]) %>% as_tibble()
    outcometemp <- x.sn %>% dplyr::select({{outcomevar}})
    bs_estimates[i]<-prop.table(table(outcometemp[[1]]))[[2]]
  }
  return(c(lower=quantile(bs_estimates,.025),est=estimate,upper=quantile(bs_estimates,.975)))
}


ALP_only <- function(datp,datc,boot="F",xvarp) {
  
  datp1 <- datp %>% dplyr::filter({{xvarp}} == 1)
  
  datpc <- rbind(datp1,datc)
  
  if (boot == "T") {
    variable = paste0("datpc$BSW",bootseq[i])
    svyds.wei  = svydesign(ids=~1, weight = ~ eval(parse(text=variable)), data = datpc)
    }
  if (boot == "F") {svyds.wei  = svydesign(ids=~1, weight = datpc$WTS_M, data = datpc)}
  
  lgtreg.w = svyglm(Formula_fit.survey, family= quasibinomial, design = svyds.wei)    
  
  beta.w = summary(lgtreg.w)$coeff[,1]   
  datpc$ps.w = lgtreg.w$fitted.values 
  
  datpc1 <- datpc %>% dplyr::filter(survey == 1) %>%
    mutate(ALP_weight = (1-ps.w) /ps.w)
}


Mikeformula <- function(datc1,datp_all,datpboot,boot="F") {
  
  m_sc = glm(Formula_fit.disc, family = quasibinomial(link = logit), data = datc1)  
  p_sc = predict(m_sc, newdata=datp_all, type="response")
  
  #Right hand side of formula
  if (boot == "F") {m_sp = glm(Formula_fit.msm, family = quasibinomial(link = logit), data = datp_all)}
  if (boot == "T") {m_sp = glm(Formula_fit.msm, family = quasibinomial(link = logit), data = datpboot)}
  p_sp = predict(m_sp, newdata=datp_all, type="response")
  
  pnewmsm <- ((1-p_sc)/p_sc)*(p_sp/(1-p_sp))
  
  return(pnewmsm)
}

Imputation <- function(M = 10,datp_all,datc,outcomevar,boot="F",pnewmsm) {
  imp <- c()
  for (i in 1:M) {
    datp_all2 <- datp_all
    datp_all2$newmsm<-ifelse(datp_all2$msm == 1,1,rbinom(size=1,n=nrow(datp_all2),prob=pnewmsm))
    
    rm.datp <- "pnewmsm"
    datp_all_new <- datp_all2 %>% dplyr::filter(newmsm == 1) %>%
      dplyr::select(-any_of(rm.datp))
    
    rm.datc2 <- c("ps.w","ALP_weight")
    
    datc2 <- datc %>% dplyr::select(-any_of(rm.datc2)) %>%
      mutate(newmsm = 1)
    
    if (boot == "F") {newmsmdat <- ALP_only(datp=datp_all_new,datc=datc2,xvarp=newmsm)}
    if (boot == "T") {newmsmdat <- ALP_only(datp=datp_all_new,datc=datc2,xvarp=newmsm,boot="T")}
    
    datcallnewALP <- newmsmdat %>% dplyr::filter(survey == 1)
    datcallnewALP1 <- datcallnewALP %>% drop_na({{outcomevar}})
    outcometemp <- datcallnewALP1 %>% dplyr::select({{outcomevar}})
    outcome2 <- as.numeric(outcometemp[[1]])
    weight <- datcallnewALP1 %>% dplyr::select(c(ALP_weight))
    imp[i]<-sum(outcome2*weight[[1]])/sum(weight[[1]])
  }
  return(imp)
}

Imputation_cchs <- function(M = 10,datp_all,datc,outcomevar,boot="F",pnewmsm,index=i) {
  imp <- c()
  for (j in 1:M) {
    datp_all2 <- datp_all
    datp_all2$newmsm<-ifelse(datp_all2$msm == 1,1,rbinom(size=1,n=nrow(datp_all2),prob=pnewmsm))
    
    newmsmdat <- datp_all2 %>% dplyr::filter(newmsm == 1) %>% drop_na({{outcomevar}})
    outcometemp <- newmsmdat %>%
      dplyr::select({{outcomevar}}) 
    outcome2 <- as.numeric(outcometemp[[1]])

    if (boot == "F") {weight <- newmsmdat$WTS_M }
    if (boot == "T") {
      variable = paste0("newmsmdat$BSW",bootseq[index])
      weight <- eval(parse(text=variable))
      }
    imp[j]<-sum(outcome2*weight)/sum(weight)
  }
  return(mean(imp))
}

main_function_cchs <- function(outcomevar,datpmsm,datc,datp) {
  #1. CCHS survey weights
  newdat <- datpmsm %>% drop_na({{outcomevar}})
  outcomeall <- newdat %>% dplyr::select({{outcomevar}})
  estimate1 <- sum(outcomeall[[1]]*newdat$WTS_M)/sum(newdat$WTS_M)
  
  #2. two step method
  pnewmsm <- Mikeformula(datc1=datc,datp_all=datp)
  result.two.step <- Imputation_cchs(M=1,datp_all=cchsall,outcomevar={{outcomevar}},datc=datc,pnewmsm=pnewmsm)
  estimate2 <- mean(result.two.step)
  
  bs_estimates<-c()
  result.two.step<-c()
  
  for (i in 1:B) {
    #1. CCHS survey weights
    variable = paste0("newdat$BSW",bootseq[i])
    weight <- eval(parse(text=variable))
    bs_estimates[i]<- sum(outcomeall[[1]]*weight)/sum(weight)
    
    #2. two step method
    set.seed(i)
    x.sp<-getboot(dat=cchsall,B=1)
    pnewmsm <- Mikeformula(datc1=datc,datp_all=cchsall,datpboot=x.sp[[1]],boot="T")
    result.two.step[i] <- Imputation_cchs(M=1,datp_all=cchsall,datc=datc,outcomevar={{outcomevar}},boot="T",pnewmsm=pnewmsm,index=i)
  }
  cchs_uw<-c(lower=quantile(bs_estimates,.025),est=estimate1,upper=quantile(bs_estimates,.975))
  twostep<-c(lower=quantile(result.two.step,.025),est=estimate2,upper=quantile(result.two.step,.975))
  return(rbind(cchs_uw,twostep))
}

main_function <- function(datc,datp,outcomevar) {
  #Point estimates:
  
  #1. Unweighted 
  outcomeall <- datc %>% dplyr::select({{outcomevar}})
  estimate1<- prop.table(table(outcomeall[[1]]))[[2]]
  
  #2. ALP only
  datALP1 <- ALP_only(datp=datp,datc=datc,xvarp=msm)
  temp <- datALP1 %>% drop_na({{outcomevar}})
  outcometemp <- temp %>% dplyr::select({{outcomevar}}) 
  estimate2 = sum(outcometemp*temp$ALP_weight,na.rm=TRUE)/sum(temp$ALP_weight,na.rm=TRUE)
  
  #3. Two-step
  pnewmsm <- Mikeformula(datc1=datc,datp_all=datp)
  result.two.step <- Imputation(M=1,datp_all=cchsall,outcomevar={{outcomevar}},datc=datc,pnewmsm=pnewmsm)
  estimate3 = mean(result.two.step)
  
  bs_estimates<-c()
  imputed_estimate<-c()
  result.two.step<-c()
  
  #Bootstraps:
  
  for (i in 1:B) {
    #1. unweighted
    outcometemp <- x.sn[[i]] %>% dplyr::select({{outcomevar}})
    bs_estimates[i]<-prop.table(table(outcometemp[[1]]))[[2]]
    
    #2. ALP only
    datALP <- ALP_only(datp=datp,datc=x.sn[[i]],boot="T",xvarp=msm)
    temp <- datALP %>% drop_na({{outcomevar}})
    outcometemp <- temp %>% dplyr::select({{outcomevar}}) 
    imputed_estimate[i] = sum(outcometemp*temp$ALP_weight)/sum(temp$ALP_weight)
    
    #3. Two-step approach
    #There is a memory problem saving all cchs bootstrap data, 
    #so jsut take x.sp within this loop once
    set.seed(i)
    x.sp<-getboot(dat=cchsall,B=1)
    pnewmsm <- Mikeformula(datc1=datALP,datp_all=cchsall,datpboot=x.sp[[1]],boot="T")
    result.two.step[i] <- Imputation(M=1,datp_all=cchsall,datc=datALP,outcomevar={{outcomevar}},boot="T",pnewmsm=pnewmsm)
  }
  uw<-c(lower=quantile(bs_estimates,.025),est=estimate1,upper=quantile(bs_estimates,.975))
  alp<-c(lower=quantile(imputed_estimate,.025),est=estimate2,upper=quantile(imputed_estimate,.975))
  twostep<-c(lower=quantile(result.two.step,.025),est=estimate3,upper=quantile(result.two.step,.975))
  return(rbind(uw,alp,twostep))
}

bootstrap_ci <- function(boot.matrix,point.est) {
  ##Want to calculate M standard errors and M estimates from each column
  var_impute_boot<-apply(boot.matrix,2,var)
  mean_impute_boot<-apply(boot.matrix,2,mean)
  #var is just the variance formula in R with n-1 (also B-1 for bootstrap, so this is correct) 
  #margin = 2 for columns
  
  #Then calculate standard MI CI var(MI)
  var_impute <- sum(var_impute_boot)/M + (M+1)*sum((mean_impute_boot-point.est)^2)/(M*(M-1))
  # And then take SE
  sd_impute <- sqrt(var_impute)
  
  #W is just the average variance from all imputations
  W = mean(var_impute_boot)
  V = sum((mean_impute_boot-point.est)^2) / (M-1)
  
  #With t dist and R degrees of freedom
  R <- (M-1)*(1 + ((M*W)/ ((M+1)*V)))^2
  alpha = 0.05
  t.score = qt(p=alpha/2, df=R,lower.tail=F)
  
  lower = point.est - t.score * sd_impute
  upper = point.est + t.score * sd_impute
  return(c(lower=lower,est=point.est,upper=upper))
}

