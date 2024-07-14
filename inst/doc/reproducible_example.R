## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,message=FALSE ,warning=FALSE, tidy.opts=list(width.cutoff=60),tidy=TRUE)

## -----------------------------------------------------------------------------
library(sgs)
groups = c(rep(1:20, each=3),
           rep(21:40, each=4),
           rep(41:60, each=5),
           rep(61:80, each=6),
           rep(81:100, each=7))

data = gen_toy_data(p=500, n=400, groups = groups, seed_id=3)

## -----------------------------------------------------------------------------
model = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", lambda = 0.5, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE, screen=TRUE)

## -----------------------------------------------------------------------------
model_path = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", lambda = "path", path_length = 5, min_frac = 0.1, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE, screen=TRUE)

## -----------------------------------------------------------------------------
model$beta[model$selected_var+1] # the +1 is to account for the intercept
model$group_effects[model$selected_grp]
model$selected_var
model$selected_grp

## -----------------------------------------------------------------------------
fdr_sensitivity = function(fitted_ids, true_ids,num_coef){
  # calculates FDR, FPR, and sensitivity
  num_true = length(intersect(fitted_ids,true_ids))
  num_false = length(fitted_ids) - num_true
  num_missed = length(true_ids) - num_true
  num_true_negatives = num_coef - length(true_ids)
  out=c()
  out$fdr = num_false / (num_true + num_false)
  if (is.nan(out$fdr)){out$fdr = 0}
  out$sensitivity = num_true / length(true_ids)
  if (length(true_ids) == 0){
    out$sensitivity = 1
  }
  out$fpr = num_false / num_true_negatives
  out$f1 = (2*num_true)/(2*num_true + num_false + num_missed)
  if (is.nan(out$f1)){out$f1 = 1}
  return(out)
}

## -----------------------------------------------------------------------------
fdr_sensitivity(fitted_ids = model$selected_var, true_ids = data$true_var_id, num_coef = 500)
fdr_sensitivity(fitted_ids = model$selected_grp, true_ids = data$true_grp_id, num_coef = 100)

## -----------------------------------------------------------------------------
cv_model = fit_sgs_cv(X = data$X, y = data$y, groups=groups, type = "linear", path_length = 20, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=TRUE, screen = TRUE)

## -----------------------------------------------------------------------------
print(cv_model)

## -----------------------------------------------------------------------------
cv_model$best_lambda_id

## -----------------------------------------------------------------------------
fdr_sensitivity(fitted_ids = cv_model$fit$selected_var, true_ids = data$true_var_id, num_coef = 500)
fdr_sensitivity(fitted_ids = cv_model$fit$selected_grp, true_ids = data$true_grp_id, num_coef = 100)

## -----------------------------------------------------------------------------
plot(cv_model,how_many = 10)

## -----------------------------------------------------------------------------
predict(model,data$X)[1:5]
dim(predict(cv_model,data$X))

## -----------------------------------------------------------------------------
sigmoid = function(x) {
  1 / (1 + exp(-x))
}
y = ifelse(sigmoid(data$X %*% data$true_beta + rnorm(400))>0.5,1,0)
train_y = y[1:350] 
test_y = y[351:400]
train_X = data$X[1:350,] 
test_X = data$X[351:400,]

## -----------------------------------------------------------------------------
cv_model = fit_sgs_cv(X = train_X, y = train_y, groups=groups, type = "logistic", path_length = 20, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, min_frac = 0.05, standardise="l2",intercept=FALSE,verbose=TRUE, screen = TRUE)

## -----------------------------------------------------------------------------
predictions = predict(cv_model,test_X)

## -----------------------------------------------------------------------------
predictions$response[1:5,cv_model$best_lambda_id]
predictions$class[1:5,cv_model$best_lambda_id]
sum(predictions$class[,cv_model$best_lambda_id] == test_y)/length(test_y)

## -----------------------------------------------------------------------------
groups = c(rep(1:20, each=3),
           rep(21:40, each=4),
           rep(41:60, each=5),
           rep(61:80, each=6),
           rep(81:100, each=7))

data = gen_toy_data(p=500, n=400, groups = groups, seed_id=3)

## -----------------------------------------------------------------------------
model = fit_gslope(X = data$X, y = data$y, groups = groups, type="linear", lambda = 0.5, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE, screen=TRUE)

## -----------------------------------------------------------------------------
screen_time = system.time(model_screen <- fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", path_length = 100, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE, screen=TRUE))
no_screen_time = system.time(model_no_screen <- fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", path_length = 100, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE, screen=FALSE))
screen_time
no_screen_time

## -----------------------------------------------------------------------------
screen_time = system.time(model_screen <- fit_gslope(X = data$X, y = data$y, groups = groups, type="linear", path_length = 100, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE, screen=TRUE))
no_screen_time = system.time(model_no_screen <- fit_gslope(X = data$X, y = data$y, groups = groups, type="linear", path_length = 100, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE, screen=FALSE))
screen_time
no_screen_time

