##----------------------------------
##  Code Version 1.0
##  This simulation code for self-paced learning with L1 penalty.
##  SPL
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on July 27, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
## Setting path ##
WhereAmI <- "D:/path"
source(paste0(WhereAmI,"Function.R"))

## Load package ##
library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)

# Input X and Y
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

##--------------------------------------
## Self-paced Learning with L1 penalty
##--------------------------------------
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
best_numcoef <- results$Coef.num[[best.iter]]

## Record best result
results.best <- list("valpredmatrix" = best_valperf, 
                     "evlpredmatrix" = best_evlperf,
                     "Coef" = best_coef, 
                     "Coef.num" = best_numcoef)


