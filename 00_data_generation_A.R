#Program to generate the simulated data, FP and FP1, it will be called in the main program

library(survey) 
library(sampling)
library(dplyr)
library(stats)
library(MASS)
set.seed(123456)
N = 1000000
Nsimr = 100
M = 10
prev = 0.10

A <-  rbinom(N, size = 1, prob = prev)  #### status in the paper 

u1 = rbinom(N, size = 1, prob = 0.5)
u2 <- runif(N, min = 0, max = 2)
u3 <- rexp(N, rate = 1)
u4 <- rchisq(N, df = 4)

x1 = u1 
x2 = u2 + (0.3*x1)
x3 = u3 + 0.2 * (x1 + x2)
x4 = u4 + 0.1 * (x1 + x2 + x3)

mu_i = -x1 -x2 + x3 + x4

y = rnorm(n=N,mean=mu_i,sd=1)
mu <- 1/N*sum(y) 

# sens_low = 0.9
# sens_high = 0.5
  
#Try 3 scenarios
# sens_low = 0.9; sens_high = 0.5 --> real, medium. about 16% misses
# sens_low = 0.5; sens_high = 0.1 --> very high bias, extreme. 56% miss
# sens_low = 0.95; sens_high = 0.8 --> very low bias, possible. 7% misses

Ap <- ifelse(A == 1 & x4 < 7, rbinom(length(A), 1, sens_low), ifelse(A == 1 & x4 >= 7, rbinom(length(A), 1, sens_high),A)) ### reported Ap; 

table(A,Ap)

# cor(x1,x4)
# cor(Ap,x4)
# cor(Ap,y)
# cor(A,y)
# cor(x4,y)

#Create population based sample s_p
#probability sample Sp target sample size np = 12,500, with inclusion prob pi = np * qi / sum (qi)
np1 = 12500

cnst_sp1 <- uniroot(f = function(x) {
  q_i <- x + x3                                                #### this is the x_p in the paper. 
  max(q_i)/min(q_i) - 20
}, interval = c(0, 1))$root    
#### cnst_sp1 = 0.87

q_i = cnst_sp1 + x3 
pi_p = (np1*q_i/sum(q_i))
di_p = 1 / pi_p #survey weights
which(pi_p < 0)

#Combine data, this is the finite population we're drawing from
data_FP = cbind(x1,x2,x3,x4,y,pi_p,di_p,Ap,A)

#This is the subset of population we're interested in, and the parameter we're estimating
FP1 <- data_FP %>% subset(A == 1)
mean(FP1[,"y"])  
targetp = mean(FP1[,"y"])  
# E(Y | A = 1) 
# 3.967401

#Non probability survey Sc:
nc = 2500
beta = c(0.18, 0.18, -0.27, -0.10) ### Use weights from LXW 

all_covars = c("x1", "x2", "x3", "x4")
sps.vars = list(c("x1", "x2", "x3", "x4"), c("x1", "x2", "x3", "x4"))

N1 = nrow(FP1)
n_all_covars = length(all_covars) ### in this case = 4 (full set of available covariates)
X_sub = FP1[,all_covars]

#### constants for weights of sc 
const = uniroot(f = function(const) {
  pi_i <- exp(const + as.matrix(X_sub)%*%beta)     #### x_c in the paper = (x1, x2, x3, x4). 
  sum(pi_i) - nc
}, interval = c(-100, 100))$root   
# const = -3.27

pi_c = exp(const + as.matrix(X_sub)%*%beta)
di_c = 1 / pi_c
FP1 <- cbind(FP1,"pi_c"=pi_c,"di_c"=di_c)
colnames(FP1)[10:11] <- c("pi_c", "di_c")

FP1[,"pi_c"] = exp(const + as.matrix(X_sub)%*%beta)
FP1[,"di_c"] = 1 / pi_c
