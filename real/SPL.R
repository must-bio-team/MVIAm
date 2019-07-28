##----------------------------------
##  Code Version 1.0
##  The real data experiment code for self-paced learning with L1 penalty.
##  SPL
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 27, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
## Setting path ##
WhereAmI <- "/path"
source(paste0(WhereAmI,"Function.R"))

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
randidx <- sample(c(1:length(Y)),size = length(Y))
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
# tuning the parameters
lambda <- 4          # age parameter
uplambda <- 2
Iter_num <- 100

results <- selfpaced.single(splitid = 1, 
                           lambda = lambda, 
                           uplambda = uplambda, 
                           Iter_num = Iter_num)

## According to the training loss value to select best results
best.iter <- which(results$valmaps == min(results$valmaps[1:length(results$itervidx)]))
best_valperf <- results$valpredmatrix[[best.iter]]
best_evlperf <- results$evlpredmatrix[[best.iter]]
best_coef <- results$Coef[[best.iter]]
coef.idx <- which(best_coef!=0)
coef.gene <- as.character(Gene[coef.idx])
coef.value <- best_coef[coef.idx]
coef.info <- cbind(coef.gene,coef.value)

## Record best result
results.best <- list("valpredmatrix" = best_valperf, 
                     "evlpredmatrix" = best_evlperf,
                     "Coef.info" =coef.info)