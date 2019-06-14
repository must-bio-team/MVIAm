# This code is based on Self-Paced Learning using logistic regression model with penalty function (Lasso, Elastic net, SCAD and MCP).

library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)

selfpace1.rank <- function(dev_decval, dev_labels, lambda) {
  #calculate the loss
  loss = (dev_decval-dev_labels)^2	#squared error
  #loss = 1/(1+e^(-1*loss))			#logistic
  
  selecedidx = which(loss <= lambda)
  
  #return the result
  result = selecedidx
  return(result)
}

selfpaced.single <- function(splitid, lambda, uplambda, Iter_num) {
  trainidx = which(splits[,splitid]==0)
  testidx = which(splits[,splitid]==1)
  
  #the list storing the result for each iteration
  valpredmatrix = list()
  evlpredmatrix = list()
  coefmatrix = list()
  nonzerocoefmatrix = list()
  iters.V_idx = list()
  
  #the starting values
  cvfit<-cv.glmnet(X[trainidx,],Y[trainidx],alpha=1,family="binomial",type.measure = "class")
  valpred <- predict(cvfit,X[trainidx,],type="response",s="lambda.min")
  
  v.idx = selfpace1.rank(valpred, dev_labels = Y[trainidx], lambda = lambda)
  this.training.vidx = v.idx
  
  #validation set is the same as the training set
  dev_labels = Y[trainidx]
  val_labels = dev_labels
  dev_x = X[trainidx,]				
  
  #store starting values
  initial.vidx = list()
  initial.vidx[[1]] = this.training.vidx
  iters.vidx = list()	
  
  for(iter in 1:Iter_num) {
    iters.vidx[[iter]] = this.training.vidx
    
    cvfit<-cv.glmnet(dev_x[this.training.vidx,],
                     dev_labels[this.training.vidx],
                     alpha=1,family="binomial",
                     type.measure = "class")  
    valprediction <- predict(cvfit,X[trainidx,],type="response",s="lambda.min")
    coefprediction = as.vector(coef(cvfit,s="lambda.min"))
    numbernonzerocoef = length(which(coefprediction!=0))
    trprediction = valprediction
    tsprediction <- predict(cvfit,X[testidx,],type="response",s="lambda.min")
    
    #self-paced learning
    selectedidx = selfpace1.rank(trprediction, dev_labels, lambda)
    this.training.vidx = selectedidx
    
    #change the parameter accoding to the step size
    lambda = uplambda + lambda

    #store the prediction for this iteration
    coefmatrix[[iter]]= coefprediction
    nonzerocoefmatrix[[iter]] = numbernonzerocoef
    valpredmatrix[[iter]] = as.numeric(valprediction)
    evlpredmatrix[[iter]] = tsprediction
    iters.V_idx[[iter]] = selectedidx
    
  }
  
  results <- list("valpredmatrix" = valpredmatrix, 
                  "evlpredmatrix" = evlpredmatrix, 
                  "initialvals" = initial.vidx, 
                  "itervidx" = iters.V_idx, 
                  "Coef"=coefmatrix, 
                  "NumbernonzeroCoef"=nonzerocoefmatrix)
  return(results)
}

# Input X and Y
X = read.table('simulation/simu_data_meta/x_meta1.txt',header = F, sep = "\t") # import X
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

# tuning the parameters
lambda = 0.05  # age parameter
uplambda = 0.05
Iter_num = 100

results = selfpaced.single(1, lambda = lambda, uplambda = uplambda, Iter_num = Iter_num)
