# SurvDNN

R package for paper: Enhanced Deep Neural Network and Feature Importance Test for Complex Survival Data

## Installation 
- System requirement: Rtools (Windows); None (MacOS/Linux)
- In R: ```devtools::install_github("ShiyuWAN3/SurvDNN)" ```
- For MAC user, if you  can not find *libquadmath.dylib* or similar files, you need to install gcc and gfortran first, and then copy and move the folder ```/usr/local/gfortran/lib``` to ```/usr/local/lib```

## Check Package Dependancy

Check package dependancy. Please use the following code to check whether dependencies are successfully loaded.

Specifically, for DeepSurv and DeepHit, we use methods implemented in the existing package ```survivalmodels```. To use this package, appropriate version of Python should be installed.

``` r
check_dependency()
```

## A Toy Simulation for survival outcome with Geompertz Baseline Hazard

Firstly, using the following code to simulate a survival dataset:

``` r
rm(list = ls())
suppressMessages(library(SurvDNN))
library(MASS)
library(survival)
library(tidyverse)
library(simsurv)

# Sample Size, Number of Features and correlation among features
N = 1200
N_train = 1000
p_block = 10
p = 50
rho = 0.5

# Survival Simulation Parameter
alpha.t = 3
lambda.t = 0.00002
alpha.c = 1.5

# Simulate Continuous Feartures

X_Scov = data.frame(mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho))
  colnames(X_Scov) = paste("X",1:(4*p_block),sep = "")

# Simulate Binary Features
for (z in 1:p_block) {
    if (z==1){
      Z_Scov = data.frame(rbinom(N,1,0.4))
    }else{
      Z_random = data.frame(rbinom(N,1,0.4))
      Z_Scov = cbind(Z_Scov,Z_random)
    }
}
  
colnames(Z_Scov) = paste("Z",1:(p_block),sep = "")

# Simulate Full Feature Matrix
C_Scov = cbind(X_Scov,Z_Scov)

# Significant variable: Z1,Z2,X1,X11,X21,X31
  
sims_beta = c(1,2,0.5,1,rep(0,p-1))
  
names(sims_beta) = c("Z1","X1Z2",
                    paste("X",p_block+1,"Square",sep = ""),
                    paste("X",2*p_block+1,"X",3*p_block+1,sep = ""),
                    paste("X",1:(4*p_block),sep = ""),paste("Z",2:p_block,sep = ""))
  
# Generate Non-Linear Term

C_Scov[,"X1Z2"] = C_Scov[,"X1"]*C_Scov[,"Z2"]
C_Scov[,paste("X",p_block+1,"Square",sep = "")] = (C_Scov[,paste("X",p_block+1,sep = "")])^2
C_Scov[,paste("X",2*p_block+1,"X",3*p_block+1,sep = "")] = C_Scov[,paste("X",2*p_block+1,sep = "")]*C_Scov[,paste("X",3*p_block+1,sep = "")]

# Generate Survival Time Using Package simsurv

Y_Scov= simsurv(dist = "gompertz",lambdas = lambda.t,gammas = alpha.t,
                betas = sims_beta,
                x= C_Scov)

C_Scov$id=1:nrow(C_Scov)

# This is the Dataset Containing True Survival Time and Features
Full_Scov = left_join(Y_Scov,C_Scov,by = "id")

# Generate Censor Distribution

# Censor Distribution

theta = 3.5

set.seed(1234)

Full_Scov$censor_time = rweibull(N,shape = alpha.c,scale = theta)

# Generate Survival Status and Observed Survival Time

Full_Scov$eventtime1 = apply(Full_Scov,1,function(x){
    Cens = x["censor_time"]
    Event = x["eventtime"]
    return(min(c(Cens,Event)))
})

Full_Scov$status1 = apply(Full_Scov,1,function(x){
    Cens = x["censor_time"]
    Event = x["eventtime"]
    if (Event >= Cens){
        return(0)
    } else {
        return(1)
    }
})

# Generate the dnnetSurvInput object

Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Scov$eventtime1[1:N_train],
                             e = Full_Scov$status1[1:N_train])
Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Scov$eventtime1[(N_train+1):N],
                             e = Full_Scov$status1[(N_train+1):N])
```

## Define the Hyperparameters Used for SurvDNN

``` r
esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
               l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
               n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE)
```

## PermFIT for SurvDNN

``` r
PermSurvDNN = permfit_survival(train = Train_Scov,n_perm =100,method = "ensemble_dnnet",
            k_fold = 5,
            n.ensemble = 100, esCtrl = esCtrl) %>% try()
View(PermSurvDNN@importance)
```

## PermFIT for XGboost

``` r
PermXGboost = permfit_survival(train = Train_Scov,n_perm =100,method = "Xgboost",
                                   k_fold = 5,nrounds = 50) %>% try()
```

## PermFIT for RSF

``` r
PermRSF = permfit_survival(train = Train_Scov,n_perm =100,method = "random_forest",
                 k_fold = 5,ntrees = 500) %>% try()
```

## PermFIT for Cox

``` r
Perm_cox_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "survival_cox",
                k_fold = 5) %>% try()
```

## PermFIT for DeepHit

``` r
PermDeepHit = permfit_survival(train = Train_Scov,n_perm =100,method = "DeepHit",
                            k_fold = 5,activation = "relu",
                            frac = 0.2,early_stopping = T,
                            num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                            batch_size = 50) %>% try()
```

## PermFIT for DeepSurv


``` r 
PermDeepSurv = permfit_survival(train = Train_Scov,n_perm =100,method = "DeepSurv",
                            k_fold = 5,activation = "relu",
                            frac = 0.2,early_stopping = T,
                            num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                            batch_size = 50) %>% try()
```

## Fit SurvDNN Model and Predict Relative Risk Using This Model

``` r 
SurvDNN = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                            object = Train_Scov,
                            n.ensemble = 100, esCtrl = esCtrl) %>% try()

centralized_model = model.centralize(SurvDNN,Train_Scov) %>% try()

# Predict Relative Risk Using the Test Dataset

Predict_SurvDNN = predict(centralized_model[[1]],Valid_Scov@x) %>% try() 

# Calculate the C-index

Cindex= try(1 - rcorr.cens(Predict_SurvDNN,Surv(Valid_Scov@y,Valid_Scov@e))[1])
```

## Predict Survival Probability

``` r
# Predicted Survival Probability of Participants in the Test Dataset
surv_df = predict_surv_df(centralized_model,Valid_Scov)
View(surv_df)
```