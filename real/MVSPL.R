##----------------------------------
##  Code Version 1.0
##  The real data experiment code for Multi-View Self-Paced Learning.
##  MVSPL
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 27, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
## Setting path ##
WhereAmI <- "/path"
source(paste0(WhereAmI,"Function.R"))

## load libraries
library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)

##### Input X and Y #####
X_view1 <- readRDS(file = paste0(WhereAmI,"data/Breast_meta_view1.rds"))
Gene <- X_view1[,1]
X_view1 <- data.matrix(t(X_view1[,-1]))

X_view2 <- readRDS(file = paste0(WhereAmI,"data/Breast_meta_view2.rds"))
Gene <- X_view2[,1]
X_view2 <- data.matrix(t(X_view2[,-1]))

X <- list()
X[[1]] <- X_view1
X[[2]] <- X_view2
View_num <- 2

Y <- read.table(file = paste0(WhereAmI,"/data/Breast_meta_labels.txt"), header = TRUE, sep = "\t")
Samp <- Y[,1]
Y <- data.matrix(Y[,-1])

##-------------------------------------------------------
## Random select samples and Setting training and test data
##-------------------------------------------------------
randidx <- sample(c(1:length(Y)),size=length(Y))
splits <- vector(mode = "numeric", length = length(Y))
splits_trainidx <- randidx[1:(0.7*length(randidx))]
splits_testidx <- randidx[(0.7*length(randidx)+1):length(randidx)]
splits[splits_trainidx] <- 0
splits[splits_testidx] <- 1
splits <- as.matrix(splits)

trainidx <- which(splits[,1]==0)
testidx <- which(splits[,1]==1)

X_train <- list()
X_test <- list()
for(i in 1:View_num){
  X_train[[i]] <- X[[i]][trainidx,]
  X_test[[i]] <- X[[i]][testidx,]
}
Y_train <- Y[trainidx]
Y_test <- Y[testidx]

##------------------------------------------------------
## Training the classifier
##------------------------------------------------------
#The list storing the result for each iteration for classifier1
iter_num <- 100
gamma <- 0.2
lambda <- c(0.05, 0.05)      # age parameter
uplambda <- 0.05
num_add <- 4            
num_up <- 2

valpredmatrix <- list()
evlpredmatrix <- list()
coefmatrix <- list()
nonzerocoefmatrix <- list()
coef_idx <- list()
coef_value <- list()
selectedidx <- list()

loss <- matrix(0, nrow = length(trainidx), ncol = View_num)
v_iter <- matrix(0, nrow = length(trainidx), ncol = View_num)

for(iter in 1:iter_num) {
  valpredmatrix[[iter]] <- matrix(0, nrow = length(trainidx), ncol = View_num)
  evlpredmatrix[[iter]] <- matrix(0, nrow = length(testidx), ncol = View_num)
  coefmatrix[[iter]] <- list()
  nonzerocoefmatrix[[iter]] <- matrix(0, nrow = 1, ncol = View_num)
}

val_labels <- matrix(rep(Y_train, each = View_num),ncol = View_num, byrow = T)
evl_labels <- matrix(rep(Y_test, each = View_num),ncol = View_num, byrow = T)
valmaps <- replicate(iter_num,0)
evlmaps <- replicate(iter_num,0)

##-------------------------------------
## Step 2.1: Initialization classifier
##-------------------------------------
for(i in 1:View_num){
  cvfit <- cv.glmnet(x = X_train[[i]],
                     y = Y_train,
                     alpha = 1,
                     family = "binomial",
                     type.measure = "class")
  valpredmatrix[[1]][,i] <- predict(cvfit, X_train[[i]], type = "response", s = "lambda.min")
}

##-----------------------
## Step 2.2: Optimization
##-----------------------
for (iter in 1:iter_num){
  
  if(length(unlist(selectedidx)) == (length(Y_train)*View_num)){break}
  cat("\nStarting the ",iter,"-th iteration.\n", sep = "")
  
  for(j in 1:View_num){
    # update v_view
    dev_decval = valpredmatrix[[iter]]
    v_iter = MVSPL.rank(dev_decval = dev_decval, 
                        dev_labels = Y_train, 
                        v_iter = v_iter, View_id = j, 
                        lambda = lambda, gamma = gamma, 
                        View_num = View_num,num_add = num_add)
    
    for(i in 1:View_num){
      selectedidx[[i]] = which(v_iter[,i]==1)
    }
    
    # update w_view Logistic with Lasso or Elasitc net
    train.idx = selectedidx[[j]]
    cvfit<-cv.glmnet(x = data.matrix(X_train[[j]][train.idx,]),
                     y = data.matrix(Y_train[train.idx]),
                     alpha = 1,
                     family = "binomial",
                     type.measure = "class") 
    valprediction <- predict(cvfit, X_train[[j]], type = "response", s = "lambda.min")
    tsprediction <- predict(cvfit, X_test[[j]], type = "response", s = "lambda.min")
    coefprediction <- as.vector(coef(cvfit, s = "lambda.min")[-1])
    numbernonzerocoef <- length(which(coefprediction!=0))
    
    valpredmatrix[[iter]][,j] <- as.numeric(valprediction)
    evlpredmatrix[[iter]][,j] <- tsprediction
    coefmatrix[[iter]][[j]] <- coefprediction
    nonzerocoefmatrix[[iter]][,j] <- numbernonzerocoef
  }
  
  #evaluate the training and test error
  val_loss <- sum((valpredmatrix[[iter]] - val_labels)^2)
  evl_loss <- sum((evlpredmatrix[[iter]] - evl_labels)^2)
  valmaps[iter] <- val_loss
  evlmaps[iter] <- evl_loss
  
  
  # update lambda and valpredmatrix for next iteriation
  lambda <- uplambda + lambda
  num_add <- num_add + num_up
  
  valpredmatrix[[iter+1]] <- valpredmatrix[[iter]]
  
}

##----------------------------------------------------
# Step 3: Find the run with the best valudation map
##----------------------------------------------------
## According to the training loss value to determine the best result
best.iter <- which(valmaps == min(valmaps[1:length(which(valmaps!=0))]))
best_valperf <- valpredmatrix[[best.iter]]
best_evlperf <- evlpredmatrix[[best.iter]]
best_coef <- coefmatrix[[best.iter]]

best_coefidx <- NULL
for(i in 1:View_num){
  if(i == 1){
    best_coefidx <- which(best_coef[[i]]!=0)
  }else{
    best_coefidx <- sort(union(best_coefidx, which(best_coef[[i]]!=0)))
  }
}
coef.gene <- as.character(Gene[best_coefidx])
coef.value <- best_coef[[1]][best_coefidx] + best_coef[[2]][best_coefidx]
coef.info <- cbind(coef.gene,coef.value)

## The final prediction label
valpred <- predict.label(predictlabels = best_valperf, truelabels = Y_train, View_num)
evlpred <- predict.label(predictlabels = best_evlperf, truelabels = Y_test, View_num)

MVSPL.result <- list("best.valperf" = valpred, 
                     "best.evlperf" = evlpred,
                     "coef.info" = coef.info)
