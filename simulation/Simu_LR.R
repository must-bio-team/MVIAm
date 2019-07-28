##----------------------------------
##  Code Version 1.0
##  This simulation code for logistic regression with Elastic net penalty.
##  L_EN
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 27, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
## Setting path ##
WhereAmI <- "D:/path"

## Load package ##
library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)

##---------------
## Input X and Y
##---------------
X <- read.table(paste0(WhereAmI,"x_meta1.txt"), header = T, sep = "\t") 
Gene <- X[,1]
X <- data.matrix(t(X))

Y <- read.table(paste0(WhereAmI,"Dataset_label.txt"), header = F, sep = "\t")
Y <- data.matrix(Y[,1])

##----------------------------------------------------
## random select samples and Setting training and test
##----------------------------------------------------
randidx <- sample(c(1:length(Y)),size=length(Y))
splits <- vector(mode = "numeric", length = length(Y))
splits_trainidx <- randidx[1:(0.7*length(randidx))]          # 70% samples are treated as training samples
splits_testidx <- randidx[(0.7*length(randidx)+1):length(randidx)]      # 30% samples are treated as test samples
splits[splits_trainidx] <- 0
splits[splits_testidx] <- 1
splits <- as.matrix(splits)

trainidx <- which(splits[,1]==0)
testidx <- which(splits[,1]==1)

##---------------------
## L_EN classifier
##---------------------
cvfit <- cv.glmnet(x = X[trainidx,], 
                 y = Y[trainidx], 
                 alpha = 0.9,                  # lasso: alpha = 1, L2: alpha = 0
                 family = "binomial",  
                 type.measure = "class")   
valprediction <- predict(cvfit, X[trainidx,], type = "class", s = "lambda.min")
tsprediction <- predict(cvfit, X[testidx,], type = "class", s = "lambda.min")
coefprediction <- as.vector(coef(cvfit, s = "lambda.min")[-1])
numbernonzerocoef <- length(which(coefprediction!=0))

results <- list("valpredmatrix" = valprediction, 
                "evlpredmatrix" = tsprediction, 
                "Coef" = coefprediction, 
                "Coef.num" = numbernonzerocoef)
