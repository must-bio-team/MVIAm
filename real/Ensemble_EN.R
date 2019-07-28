##----------------------------------
##  Code Version 1.0
##  The real data experiment code for the logistic regression model with Elastic Net embeded to the ensemble framework.
##  Ensemble_EN
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 27, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
## Setting path ##
WhereAmI <- "D:/path"
source(paste0(WhereAmI,"Function.R"))


##---------------------
## load libraries
##---------------------
library(amritr)
library(glmnet)
library(Matrix)
library(tseries)
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
randidx <- sample(c(1:length(Y)),size = length(Y))
splits <- vector(mode = "numeric", length = length(Y))
splits_trainidx <- randidx[1:(0.7*length(randidx))]                     # 70% samples are treated as training samples
splits_testidx <- randidx[(0.7*length(randidx)+1):length(randidx)]      # 30% samples are treated as test samples
splits[splits_trainidx] <- 0
splits[splits_testidx] <- 1
splits <- as.matrix(splits)

trainidx <- which(splits[,1]==0)
testidx <- which(splits[,1]==1)

## Initalize ##
data.train <- list()
data.test <- list()
y_train <- NULL
y_test <- NULL

for(i in 1:View_num){
  data.train[[i]] <- X[[i]][trainidx,]
  data.test[[i]] <- X[[i]][testidx,]
  y_train <- Y[trainidx]
  y_test <- Y[testidx]
}


##-------------------------
## 2.Ensemble_EN classifier
##-------------------------
ensemblePanels <- lapply(data.train, function(i){
  cvfit_ensem <- cv.glmnet(x = i, 
                           y = y_train, 
                           alpha = 0.9,
                           family = "binomial",
                           type.measure = "class")
})

ensembleValiPredction <- mapply(function(cvfit,x){
  valprediction <- predict(cvfit, x, type = "response", s = "lambda.min")
}, cvfit = ensemblePanels, x = data.train)

ensembleTestPredction <- mapply(function(cvfit,x){
  tsprediction <- predict(cvfit, x, type = "response", s = "lambda.min")
}, cvfit = ensemblePanels, x = data.test)

ensembleCoef <- lapply(ensemblePanels, function(i){
  coefprediction <- as.vector(coef(i, s = "lambda.min")[-1])
})

Coef.idx <- NULL
for(i in 1:View_num){
  if(i == 1){
    Coef.idx <- which(ensembleCoef[[i]]!=0)
  }else{
    Coef.idx <- sort(union(Coef.idx, which(ensembleCoef[[i]]!=0)))
  }
}
coef.gene <- as.character(Gene[Coef.idx])
coef.value <- ensembleCoef[[1]][Coef.idx] + ensembleCoef[[2]][Coef.idx]
coef.inf <- cbind(coef.gene,coef.value)
coef.info <- coef.inf[order(abs(as.numeric(coef.inf[,2])),decreasing=T),]

## The final prediction label
valpred <- predict.label(predictlabels = ensembleValiPredction, truelabels = y_train, View_num)
evlpred <- predict.label(predictlabels = ensembleTestPredction, truelabels = y_test, View_num)


results <- list("valpredmatrix" = valpred, 
                "evlpredmatrix" = evlpred, 
                "coef.info" = coef.info)

