##----------------------------------
##  Code Version 1.0
##  The real data experiment code for logistic regression model with L1/Elastic net penalty.
##  L_1/LEN
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 27, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
## Setting path ##
WhereAmI <- "D:/path"

## Load packages
library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)

##------------------------
## Input X and Y
##------------------------
X <- readRDS(file = paste0(WhereAmI,"data/Breast_meta_view1.rds"))
Gene <- X[,1]
X <- data.matrix(t(X[,-1]))

Y <- read.table(file = paste0(WhereAmI,"/data/Breast_meta_labels.txt"), header = TRUE, sep = "\t")
Samp <- Y[,1]
Y <- data.matrix(Y[,-1])

##-------------------------------------------------------
## Random select samples and Setting training and test data
##-------------------------------------------------------
randidx <- sample(c(1:length(Y)), size = length(Y))
splits <- vector(mode = "numeric", length = length(Y))
splits_trainidx <- randidx[1:(0.7*length(randidx))]
splits_testidx <- randidx[(0.7*length(randidx)+1):length(randidx)]
splits[splits_trainidx] <- 0
splits[splits_testidx] <- 1
splits <- as.matrix(splits)

trainidx <- which(splits[,1]==0)
testidx <- which(splits[,1]==1)

##------------------------------------------------------
## Training the classifier
##------------------------------------------------------
cvfit <- cv.glmnet(x = X[trainidx,],
                   y = Y[trainidx],
                   alpha = 0.9,            # Lasso: alpha = 1; L2: alpha = 0
                   family = "binomial",
                   type.measure = "class")
valprediction <- predict(cvfit,X[trainidx,], type = "class", s = "lambda.min")
tsprediction <- predict(cvfit,X[testidx,],type = "class", s = "lambda.min")
coef <- as.vector(coef(cvfit, s = "lambda.min")[-1])
coef.idx <- which(coef!=0)
coef.gene <- as.character(Gene[coef.idx])
coef.value <- coef[coef.idx]
coef.inf <- cbind(coef.gene,coef.value)
coef.info <- coef.inf[order(abs(as.numeric(coef.inf[,2])),decreasing=T),]

results <- list("valpredmatrix" = valprediction, 
                "evlpredmatrix" = tsprediction, 
                "Coef.info" =coef.info)
