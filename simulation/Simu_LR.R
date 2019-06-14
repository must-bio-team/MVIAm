# This code is using logistic regression model with penalty function (Lasso, Elastic net, SCAD and MCP).

library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)


# Input X and Y
X = read.table('simulation/simu_data_meta/x_meta1.txt',header = T, sep = "\t") # import X
Gene = X[,1]
X = data.matrix(t(X[,-1]))

Y = read.table('simulation/simu_data_meta/y.txt',header = F, sep = "\t") # import y
Y = data.matrix(Y[,1])

# random select samples and Setting training and testing
randidx = sample(c(1:length(Y)),size=length(Y))
splits = vector(mode = "numeric", length = length(Y))
splits_trainidx = randidx[1:(0.7*length(randidx))]
splits_testidx = randidx[(0.7*length(randidx)+1):length(randidx)]
splits[splits_trainidx] = 0
splits[splits_testidx] = 1
splits = as.matrix(splits)

trainidx = which(splits[,1]==0)
testidx = which(splits[,1]==1)

cvfit<-cv.glmnet(X[trainidx,], Y[trainidx], alpha = 0.9, family = "binomial", type.measure = "class")
valprediction <- predict(cvfit,X[trainidx,],type="response",s="lambda.min")
tsprediction <- predict(cvfit,X[testidx,],type="response",s="lambda.min")
coefprediction = as.vector(coef(cvfit,s="lambda.min"))
numbernonzerocoef = length(which(coefprediction!=0))


results <- list("valpredmatrix" = valprediction, 
                "evlpredmatrix" = tsprediction, 
                "Coef" = coefprediction, 
                "NumbernonzeroCoef" = numbernonzerocoef)
