# This code is based on Multi-view Self-Paced Learning 
# using logistic regression model with penalty function.

library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)

selfpace1.rank <- function(dev_decval, dev_labels, v_iter, lambda, View_id, gamma) {
  
  # initialize the loss function
  loss = matrix(0, nrow = length(dev_labels), ncol = View_num)
  
  #calculate the loss
  for(m in 1: View_num){
    if(m != View_id){
      loss[,m] = (dev_decval[,m] - dev_labels)^2    #squared error
    }else{
      next;
    }
  }
  
  # Update v(View_num-j)
  for(m in 1:View_num){
    if(m != View_id){
      for(i in 1:length(dev_labels)){
        if(loss[i,m] < lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m])){
          v_iter[i,m] = 1
        }else{
          v_iter[i,m] = 0
        }
      }
    }
  }
  
  # Update vj
  loss[,View_id] = (dev_decval[,View_id] - dev_labels)^2
  for(i in 1:length(dev_labels)){
    if(loss[i,View_id] < lambda[View_id] + gamma * (sum(v_iter[i,])-v_iter[i,View_id])){
      v_iter[i,View_id] = 1
    }else{
      v_iter[i,View_id] = 0
    }
  }
  #return the result
  return(v_iter)
}


# Input X and Y
X_view1 = read.table('simulation/simu_data_meta/x_meta1.txt',header = F, sep = "\t")
X_view1 = data.matrix(t(X_view1[,-1]))

X_view2 = read.table('simulation/simu_data_meta/x_meta2.txt',header = F, sep = "\t")
X_view2 = data.matrix(t(X_view2[,-1]))

X = list()
X[[1]] = X_view1
X[[2]] = X_view2

Y = read.table('simulation/simu_data_meta/y.txt',header = F, sep = "\t")
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

#The list storing the result for each iteration for classifier1
View_num = 2
iter_num = 100
gamma = 0.2
lambda = c(0.05, 0.05)  # age parameter
uplambda = 0.05

valpredmatrix = list()
evlpredmatrix = list()
coefmatrix = list()
nonzerocoefmatrix = list()
iters.V_idx = list()
loss = matrix(0, nrow = length(trainidx), ncol = View_num)

#Initial v
v_iter = matrix(0, nrow = length(trainidx), ncol = View_num)

for(iter in 1:iter_num) {
  valpredmatrix[[iter]] = matrix(0, nrow = length(trainidx), ncol = View_num)
  evlpredmatrix[[iter]] = matrix(0, nrow = length(testidx), ncol = View_num)
  coefmatrix[[iter]] =  matrix(0, nrow = length(X[[1]][1,])+1, ncol = View_num)
  nonzerocoefmatrix[[iter]] = matrix(0, nrow = 1, ncol = View_num)
  iters.V_idx[[iter]] = list()
}

# Step1: Initialization
for(i in 1:View_num){
  cvfit<-cv.glmnet(X[[i]][trainidx,],
                   Y[trainidx],
                   alpha = 1,
                   family = "binomial",
                   type.measure = "class")
  valprediction <- predict(cvfit,X[[i]][trainidx,],type="response",s="lambda.min")
  coefprediction = as.vector(coef(cvfit,s="lambda.min"))
  numbernonzerocoef = length(which(coefprediction!=0))
  valpredmatrix[[1]][,i] = as.numeric(valprediction)
  coefmatrix[[1]][,i] = coefprediction
  nonzerocoefmatrix[[1]][,i] = numbernonzerocoef
}

dev_x = list()
dev_labels = Y[trainidx]
for(q in 1:View_num){
  dev_x[[q]] = X[[q]][trainidx,]
}

# Step2: Optimization
for (iter in 1:iter_num){
  for(j in 1:View_num){
    # update v_view
    dev_decval = valpredmatrix[[iter]]
    v_iter = selfpace1.rank(dev_decval, dev_labels = Y[trainidx], v_iter, lambda, View_id = j,gamma = gamma)
    
    selectedidx = list()
    for(i in 1:View_num){
      selectedidx[[i]] = which(v_iter[,i]==1)
    }
    
    # update w_view Logistic with Lasso or Elasitc net
    trainidx_iter = selectedidx[[j]]
    cvfit<-cv.glmnet(dev_x[[j]][trainidx_iter,],
                     dev_labels[trainidx_iter],
                     alpha = 1,
                     family = "binomial",
                     type.measure = "class") 
    valprediction <- predict(cvfit,X[[j]][trainidx,],type="response",s="lambda.min")
    tsprediction <- predict(cvfit,X[[j]][testidx,],type="response",s="lambda.min")
    coefprediction = as.vector(coef(cvfit,s="lambda.min"))
    numbernonzerocoef = length(which(coefprediction!=0))
    
    valpredmatrix[[iter]][,j] = as.numeric(valprediction)
    evlpredmatrix[[iter]][,j] = tsprediction
    coefmatrix[[iter]][,j] = coefprediction
    nonzerocoefmatrix[[iter]][,j] = numbernonzerocoef
    iters.V_idx[[iter]][[j]] = trainidx_iter
  }
  
  # update lambda and valpredmatrix for next iteriation
  lambda= uplambda + lambda
  
  if(iter < iter_num){
    valpredmatrix[[iter+1]]=valpredmatrix[[iter]]
  }
}

results <- list("valpredmatrix" = valpredmatrix, 
                "evlpredmatrix" = evlpredmatrix, 
                "itervidx" = iters.V_idx, 
                "Coef"=coefmatrix, 
                "NumbernonzeroCoef"=nonzerocoefmatrix)
