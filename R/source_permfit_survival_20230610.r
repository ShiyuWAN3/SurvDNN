library(survival)
library(tidyverse)
library(simsurv)
library(MASS)
library(deepTL)
library(survivalmodels)
library(randomForestSRC)
library(survivalsvm)
library(reticulate)
library(Hmisc)
library(xgboost)
library(survivalsvm)
library(pheatmap)

##### 0. Args Setting #####


#' Cox-PH log-likelihood calculation
#'
#' @param t_threshold Calculate the risk set at time t_throshold
#' @param times Patients' event times
#'
#' @return Returns the IDs of patients whose event time is longer then t_threshold
#' @export
#'
#' @examples
#'
#'
risk.set <- function(t_threshold,times) {
  return(which(times >= t_threshold))
}



#' Calculated the Cox-PH's log-likelihood based on event status, event time and risk prediction
#'
#' @param Status Patients' event statuses
#' @param Times Patients' event times
#' @param f_hat_y X \%*\% Beta
#'
#' @return Returns the log-likelihood of a Cox-PH Model
#' @export
#'
#' @examples
#'
loglik_coxph = function(Status,Times,f_hat_y){
  if (length(Status) != length(Times)){
    return("Check your input!")
  } else {
    loglik = 0
    Xbeta = f_hat_y %>% as.numeric() # 被替换成为 predict value的地方
    df= data.frame(status=Status,
                   times=Times,
                   yhat=f_hat_y)
    df$yhat[which(is.infinite(df$yhat))] = NA
    df$loglik = apply(df,1,function(x){
      status = x["status"]
      time = x["times"]
      riskset= risk.set(t_threshold = time,times = df$times)
      yhat_i = x["yhat"]
      yhat_riskset= df$yhat[riskset]
      return(as.numeric(status)*(yhat_i-log(sum(exp(yhat_riskset),na.rm = T))))
    })
    df$loglik[which(is.infinite(df$loglik))]=NA
    loglik = sum(df$loglik,na.rm = T)
  }
  if (is.infinite(loglik)){
    return(NA)
  }else{
    return(loglik)
  }
}

mod_permfit <- function(method, model.type, object, ...) {
  if (model.type == "survival"){
    if (method == "random_forest"){
      mod <- do.call(randomForestSRC::rfsrc,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "survival_aft"){
      df_train=cbind(object@y,object@e,object@x) %>% as.data.frame()
      colnames(df_train)=c("EventTime","EventStatus",
                           colnames(object@x))
      mod <- survival::survreg(Surv(EventTime,EventStatus)~.,data = df_train,dist = "logistic")
    } else if (method == "survival_cox") {
      df_train=cbind(object@y,object@e,object@x) %>% as.data.frame()
      colnames(df_train)=c("EventTime","EventStatus",
                           colnames(object@x))
      mod <- survival::coxph(Surv(EventTime,EventStatus)~.,data = df_train)
    } else if (method == "DeepSurv") {
      mod <- do.call(survivalmodels::deepsurv,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "DeepHit") {
      mod <- do.call(survivalmodels::deephit,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "Xgboost"){
      label_train = ifelse(object@e == 1,object@y,-object@y)
      df_train_xg = xgboost::xgb.DMatrix(data = as.matrix(object@x),label = label_train)
      mod =do.call(xgboost::xgboost,
                   appendArg(appendArg(list(...),"data",df_train_xg,TRUE),
                             "objective","survival:cox",TRUE))
    } else if (method == "Survival_SVM"){
      mod <- do.call(survivalsvm::survivalsvm,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "DnnSurv") {
      mod <- do.call(survivalmodels::dnnsurv,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "lasso") {
      lasso_family <- "cox"
      cv_lasso_mod <- glmnet::cv.glmnet(object@x, cbind(time = object@y,status = object@e), family = lasso_family)
      mod <- glmnet::glmnet(object@x, cbind(time = object@y,status = object@e), family = lasso_family,
                            lambda = cv_lasso_mod$lambda[which.min(cv_lasso_mod$cvm)])
    } else if (method == "ensemble_dnnet") {
      mod <- do.call(deepTL::ensemble_dnnet, appendArg(list(...), "object", object, TRUE))
    } else if (method == "dnnet") {
      spli_obj = splitDnnet(object, 0.8)
      mod <- do.call(dnnet, appendArg(appendArg(list(...), "train", spli_obj$train, TRUE),
                                      "validate", spli_obj$valid, TRUE))
    } else {
      return("Survival Analysis Not Applicable")
    }
  } else {
    if(method == "ensemble_dnnet") {
      mod <- do.call(ensemble_dnnet, appendArg(list(...), "object", object, TRUE))
    } else if (method == "random_forest") {
      mod <- do.call(randomForest::randomForest,
                     appendArg(appendArg(list(...), "x", object@x, TRUE), "y", object@y, TRUE))
    } else if (method == "lasso") {
      lasso_family <- ifelse(model.type == "regression", "gaussian",
                             ifelse(model.type == "binary-classification", "binomial", "cox"))
      cv_lasso_mod <- glmnet::cv.glmnet(object@x, object@y, family = lasso_family)
      mod <- glmnet::glmnet(object@x, object@y, family = lasso_family,
                            lambda = cv_lasso_mod$lambda[which.min(cv_lasso_mod$cvm)])
    } else if (method == "linear") {
      if(model.type == "regression") {
        mod <- stats::lm(y ~ ., data.frame(x = object@x, y = object@y))
      } else if(model.type == "binary-classification") {
        mod <- stats::glm(y ~ ., family = "binomial", data = data.frame(x = object@x, y = object@y))
      } else {
        mod <- survival::coxph(survival::Surv(y, e) ~ ., data = data.frame(x = object@x, y = object@y, e = object@e))
      }
    } else if (method == "svm") {
      if(model.type == "regression") {
        mod <- e1071::tune.svm(object@x, object@y, gamma = 10**(-(0:4)), cost = 10**(0:4/2),
                               tunecontrol = e1071::tune.control(cross = 5))
        mod <- mod$best.model
      } else if(model.type == "binary-classification") {
        mod <- e1071::tune.svm(object@x, object@y, gamma = 10**(-(0:4)), cost = 10**(0:4/2),
                               tunecontrol = e1071::tune.control(cross = 5))
        mod <- svm(object@x, object@y, gamma = mod$best.parameters$gamma, cost = mod$best.parameters$cost, probability = TRUE)
      } else {
        return("Not Applicable")
      }
    } else if (method == "dnnet") {
      spli_obj <- splitDnnet(object, 0.8)
      mod <- do.call(dnnet, appendArg(appendArg(list(...), "train", spli_obj$train, TRUE),
                                      "validate", spli_obj$valid, TRUE))
    } else {
      return("Not Applicable")
    }
  }
  return(mod)
}

predict_mod_permfit <- function(mod, object, method, model.type) {
  if(model.type == "regression") {
    if(!method %in% c("linear", "lasso")) {
      return(predict(mod, object@x))
    } else if(method == "linear") {
      return(predict(mod, object@x)[, "s0"])
    } else {
      return(predict(mod, object@x)[, "s0"])
    }
  } else if(model.type == "binary-classification") {
    if(method %in% c("dnnet")) {
      return(predict(mod, object@x)[, mod@label[1]])
    } else if(method == "ensemble_dnnet") {
      return(predict(mod, object@x)[, mod@model.list[[1]]@label[1]])
    } else if(method == "random_forest") {
      return(predict(mod, object@x, type = "prob")[, 1])
    } else if (method == "lasso") {
      return(1 - predict(mod, object@x, type = "response")[, "s0"])
    } else if (method == "linear") {
      return(1 - predict(mod, data.frame(x = object@x, y = object@y), type = "response"))
    } else if (method == "svm") {
      return(attr(predict(mod, object@x, decision.values = TRUE, probability = TRUE),
                  "probabilities")[, levels(object@y)[1]])
    }
  } else if(model.type == "survival") {
    if(method %in% c("ensemble_dnnet", "dnnet")) {
      pred = predict(mod, object@x)
      return(list(pred,pred))
    } else if(method == "random_forest") {
      pred = as.numeric(predict(mod, data.frame(x = object@x))$predicted)
      return(list(pred,pred))
    } else if (method == "lasso") {
      pred = predict(mod, object@x)[, "s0"] %>% as.numeric() %>% exp()
      return(list(pred,log(pred)))
    } else if (method == "survival_cox") {
      pred = predict(mod, data.frame(object@x), type = "risk")
      return(list(pred,log(pred)))
    } else if (method == "survival_aft") {
      pred = -as.numeric(predict(mod, data.frame(object@x), type = "lp"))
      return(list(pred,pred))
    } else if (method == "Survival_SVM") {
      pred = -as.numeric(predict(mod, data.frame(x = object@x))$predicted)
      return(list(pred,pred))
    } else if (method == "DeepSurv" | method == "DeepHit" | method == "DnnSurv") {
      pred = predict(mod, data.frame(x = object@x),type = "risk")
      return(list(pred,log(pred)))
    } else if (method == "Xgboost"){
      pred = predict(mod,newdata = object@x,type = "risk")
      return(list(pred,log(pred)))
    }
  } else {
    return("Not Applicable")
  }
}


Cindex = function(Status,Times,f_hat_y){
  Inf_fhaty = which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)
  if (sum(which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)) == 0){
    rcorr_cens = rcorr.cens(f_hat_y,Surv(Times,Status))[1]
    Cindex = 1 - as.numeric(rcorr_cens)
    return(Cindex)
  }else{
    rcorr_cens = rcorr.cens(f_hat_y[-Inf_fhaty],Surv(Times[-Inf_fhaty],Status[-Inf_fhaty]))[1]
    Cindex = 1 - as.numeric(rcorr_cens)
    return(Cindex)
  }
}

timedAUC = function(Status,Times,f_hat_y){
  Inf_fhaty = which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)
  if (sum(which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)) == 0){
    death_case = which(Status == 1)
    choose_time = Times[death_case]
    time_range = quantile(choose_time,c(0.25,0.5,0.75))
    tdAUC = timeROC(T = Times, delta = Status,
                    marker = f_hat_y,weighting = "marginal",cause = 1,
                    times = time_range)
    return(mean(tdAUC$AUC))
  }else{
    new_time = Times[-Inf_fhaty]
    new_status = Status[-Inf_fhaty]
    new_f_hat_y = f_hat_y[-Inf_fhaty]
    death_case = which(new_status == 1)
    choose_time = new_time[death_case]
    time_range = quantile(choose_time,c(0.25,0.5,0.75))
    tdAUC = timeROC(T = new_time, delta = new_status,
                    marker = new_f_hat_y,weighting = "marginal",cause = 1,
                    times = time_range)
    return(mean(tdAUC$AUC))
  }
}

log_lik_diff <- function(model.type, y_hat, y_hat0, object,
                         y_max = 1-10**-10, y_min = 10**-10,
                         y_hatcoxl,y_hat0coxl) {
  if(model.type == "regression") {
    return((object@y - y_hat)**2 - (object@y - y_hat0)**2)
  } else if(model.type == "binary-classification") {
    y_hat <- ifelse(y_hat < y_min, y_min, ifelse(y_hat > y_max, y_max, y_hat))
    y_hat0 <- ifelse(y_hat0 < y_min, y_min, ifelse(y_hat0 > y_max, y_max, y_hat0))
    return(-(object@y == levels(object@y)[1])*log(y_hat) - (object@y != levels(object@y)[1])*log(1-y_hat) +
             (object@y == levels(object@y)[1])*log(y_hat0) + (object@y != levels(object@y)[1])*log(1-y_hat0))
  } else if(model.type == "survival") {
    return(c(Cindex(Status = object@e,Times = object@y,f_hat_y = y_hat)-
               Cindex(Status = object@e,Times = object@y,f_hat_y = y_hat0),
             loglik_coxph(Status = object@e,Times = object@y,f_hat_y = y_hatcoxl)-
               loglik_coxph(Status = object@e,Times = object@y,f_hat_y = y_hat0coxl),
             timedAUC(Status = object@e,Times = object@y,f_hat_y = y_hat)-
               timedAUC(Status = object@e,Times = object@y,f_hat_y = y_hat0)))
  } else {
    return("Not Applicable")
  }
}

##### 4. Permfit Survival #####

permfit_survival <- function(train, validate = NULL, k_fold = 5,
                             n_perm = 100, pathway_list = list(),
                             method = c("ensemble_dnnet", "random_forest",
                                        "lasso", "linear", "svm", "dnnet",
                                        "survival_cox","survival_aft","DeepSurv",
                                        "DeepHit","Xgboost","Survival_SVM",
                                        "DnnSurv")[8],
                             shuffle = NULL,
                             ...){
  n_pathway <- length(pathway_list)
  # Number of observations
  n <- dim(train@x)[1]
  # Number of Covariates
  p <- dim(train@x)[2]

  # 根据导入数据的类型，判断应该使用的模型的类型
  if(class(train) == "dnnetSurvInput") {
    model.type <- "survival"
  } else if(class(train) == "dnnetInput") {
    if(is.factor(train@y)) {
      model.type <- "binary-classification"
    } else {
      model.type <- "regression"
    }
  } else {
    stop("'train' has to be either a dnnetInput or dnnetSurvInput object.")
  }

  # 不进行交叉验证，需要一个测试集
  if(k_fold == 0) {
    if(is.null(validate)){
      stop("A validation set is required when k = 0. ")
    }
    # 测试集的样本量
    n_valid <- dim(validate@x)[1]
    # 在训练集上拟合模型
    mod = mod_permfit(method = method,model.type = model.type,object = train,...)
    # Risk Score
    f_hat_x <- predict_mod_permfit(mod = mod,object = validate,
                                   method = method,model.type = model.type)
    valid_ind <- list(1:length(validate@y))
    y_pred <- f_hat_x

    # n_pathway: 对于多个分组的变量，会生成一系列的dummy
    # 这个时候要同时进行permfit

    # To Be Done
    if(n_pathway >= 1) {
      p_score <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorel <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorea <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      for(i in 1:n_pathway) {
        for(l in 1:n_perm){
          # Validate Dataset
          set.seed(i*1000+l)
          x_i <- validate@x
          # 只permutation一个Block
          x_i[, pathway_list[[i]]] <- x_i[, pathway_list[[i]]][sample(n_valid), ]
          Perm_Validate = importDnnetSurv(x = x_i,
                                          y = validate@y,
                                          e = validate@e)
          pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                       object = Perm_Validate,method = method)
          diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                            y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                            y_hat0coxl = pred_i[[2]],object = train_spl$train)
          p_score[l, , i] = diff[1]
          p_scorel[l, , i] = diff[2]
          p_scorea[l, , i] = diff[3]
        }
      }
    }

    # PERMFIT:Core
    p_score2 <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2l <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2a <- array(NA, dim = c(n_perm, n_valid, p))

    for (i in 1:p) {
      for (l in 1:n_perm) {
        set.seed(i*1000+l)
        x_i <- validate@x
        x_i[, i] <- x_i[, i][sample(dim(x_i)[1])]
        Perm_Validate = importDnnetSurv(x = x_i,
                                        y = validate@y,
                                        e = validate@e)
        pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                     object = Perm_Validate,method = method)
        diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                          y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                          y_hat0coxl = pred_i[[2]],object = train_spl$train)
        p_score2[l, , i] = diff[1]
        p_score2l[l, , i] = diff[2]
        p_score2a[l, , i] = diff[3]
      }
    }
  } else { # K-Fold 交叉验证
    valid_ind <- list()
    # 对 n 进行不放回抽样：重新排序
    if(is.null(shuffle)) {
      shuffle <- sample(n)
    }
    # train 的观测数
    n_valid <- n
    # 生成了很多个 0
    y_pred <- numeric(length(train@y))

    ## n_pathway：针对dummy的
    if(n_pathway >= 1){
      p_score <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorel <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorea <- array(NA, dim = c(n_perm, n_valid, n_pathway))
    }

    # 初始化permutation score矩阵
    p_score2 <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2l <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2a <- array(NA, dim = c(n_perm, n_valid, p))

    # 验证集误差初始化：生成k-折个0
    valid_error <- numeric(k_fold)

    # K-Fold:: Start
    for(k in 1:k_fold){
      # 先将n个样本生成一个排列，然后按顺序取
      train_spl <- splitDnnet(train, shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)])
      # 记录所用的验证集的观测编号
      valid_ind[[k]] <- shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)]
      # Typo: Valid 是其中的 4折
      mod = mod_permfit(method = method,model.type = model.type,object = train_spl$valid,...)
      #print(mod)
      # Risk Score
      f_hat_x <- predict_mod_permfit(mod = mod,object = train_spl$train,
                                     method = method,model.type = model.type)

      #print(f_hat_x)
      # valid_error[k] <- sum(log_lik_diff(model.type, f_hat_x, f_hat_x, train_spl$train))

      y_pred[valid_ind[[k]]] <- f_hat_x[1]

      if(k == 1) {
        final_model <- mod
      } else if(method == "ensemble_dnnet") {
        final_model@model.list <- c(final_model@model.list, mod@model.list)
        final_model@loss <- c(final_model@loss, mod@loss)
        final_model@keep <- c(final_model@keep, mod@keep)
      }

      # To Be Done
      if(n_pathway >= 1){
        for(i in 1:n_pathway) {
          for(l in 1:n_perm) {
            x_i <- train_spl$train@x
            x_i[, pathway_list[[i]]] <- x_i[, pathway_list[[i]]][sample(dim(x_i)[1]), ]
            Perm_Validate = importDnnetSurv(x = x_i,
                                            y = train_spl$train@y,
                                            e = train_spl$train@e)
            pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                         object = Perm_Validate,method = method)
            diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                              y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                              y_hat0coxl = pred_i[[2]],object = train_spl$train)
            p_score[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[1]
            p_scorel[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[2]
            p_scorea[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[3]
          }
        }
      }
      # Permfit::Core
      # p: 第p个协变量

      for (i in 1:p) {
        for (l in 1:n_perm) {
          set.seed(i*1000+l)
          x_i <- train_spl$train@x
          x_i[, i] <- x_i[, i][sample(dim(x_i)[1])]
          Perm_Validate = importDnnetSurv(x = x_i,
                                          y = train_spl$train@y,
                                          e = train_spl$train@e)
          pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                       object = Perm_Validate,method = method)
          diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                            y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                            y_hat0coxl = pred_i[[2]],object = train_spl$train)
          p_score2[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[1]
          p_score2l[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[2]
          p_score2a[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[3]
        }
      }
      print(k)
    }
    mod <- final_model
    # valid_error <- sum(valid_error)/n_valid
  }

  # Importance Score: C-Index
  if(is.null(colnames(train@x))) {
    imp <- data.frame(var_name = paste0("V", 1:p))
  } else  {
    imp <- data.frame(var_name = colnames(train@x))
  }
  imp$importance <- apply(apply(p_score2, 2:3, mean), 2, mean, na.rm = TRUE)
  imp$importance_sd <- sqrt(apply(apply(p_score2, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
  imp$importance_pval <- 1 - stats::pnorm(imp$importance/imp$importance_sd)
  if(n_perm > 1) {
    imp$importance_sd_x <- apply(apply(p_score2, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
    imp$importance_pval_x <- 1 - stats::pnorm(imp$importance/imp$importance_sd_x)
  }

  # Importance Score: likelihood
  if(is.null(colnames(train@x))) {
    impl <- data.frame(var_name = paste0("V", 1:p))
  } else  {
    impl <- data.frame(var_name = colnames(train@x))
  }
  impl$importance <- apply(apply(p_score2l, 2:3, mean), 2, mean, na.rm = TRUE)
  impl$importance_sd <- sqrt(apply(apply(p_score2l, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
  impl$importance_pval <- 1 - stats::pnorm(impl$importance/impl$importance_sd)
  if(n_perm > 1) {
    impl$importance_sd_x <- apply(apply(p_score2l, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
    impl$importance_pval_x <- 1 - stats::pnorm(impl$importance/impl$importance_sd_x)
  }

  # Importance Score: tdAUC
  if(is.null(colnames(train@x))) {
    impa <- data.frame(var_name = paste0("V", 1:p))
  } else  {
    impa <- data.frame(var_name = colnames(train@x))
  }
  impa$importance <- apply(apply(p_score2a, 2:3, mean), 2, mean, na.rm = TRUE)
  impa$importance_sd <- sqrt(apply(apply(p_score2a, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
  impa$importance_pval <- 1 - stats::pnorm(impa$importance/impa$importance_sd)
  if(n_perm > 1) {
    impa$importance_sd_x <- apply(apply(p_score2a, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
    impa$importance_pval_x <- 1 - stats::pnorm(impa$importance/impa$importance_sd_x)
  }

  imp_block <- data.frame()
  if(n_pathway >= 1) {
    if(is.null(names(pathway_list))) {
      imp_block <- data.frame(block = paste0("P", 1:n_pathway))
    } else {
      imp_block <- data.frame(block = names(pathway_list))
    }
    imp_block$importance <- apply(apply(p_score, 2:3, mean), 2, mean, na.rm = TRUE)
    imp_block$importance_sd <- sqrt(apply(apply(p_score, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
    imp_block$importance_pval <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd)
    if(n_perm > 1) {
      imp_block$importance_sd_x <- apply(apply(p_score, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
      imp_block$importance_pval_x <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd_x)
    }
  }

  return(new("PermFIT", model = list(impl,impa), importance = imp, block_importance = imp_block,
             validation_index = valid_ind, y_hat = y_pred))
}
