selfpace1.rank <- function(dev_decval, dev_labels, lambda) {
  #calculate the loss
  loss = (dev_decval-dev_labels)^2	#squared error
  #loss = 1/(1+e^(-1*loss))			#logistic
  
  ## lambda as threshold
  # selecedidx = which(loss <= lambda)
  # 
  # #return the result
  # result = selecedidx
  # return(result)
  
  ## lambda as threshold by select sample numbers
  posidx = which(dev_labels==1)	#postive id mapping
  negidx = which(dev_labels==0)	#negative id mapping
  
  #calculate pos_lambda1 and neg_lambda2 according to the rank
  pos_lambda = sort(loss[posidx,1])[min(length(posidx), lambda)]
  neg_lambda = sort(loss[negidx,1])[min(length(negidx), lambda)]
  
  #it is like first sorting sampled based on the metric and then select top lambda1_rank
  if(length(unique(loss[posidx]))!=1){
    selectedposidx <- posidx[which(loss[posidx,1] <= pos_lambda)]
  }else{
    selectedposidx <- sample(posidx, size = min(lambda, length(posidx)), replace = FALSE)
  }
  if(length(unique(loss[negidx]))!=1){
    selectednegidx <- negidx[which(loss[negidx,1] <= neg_lambda)]
  }else{
    selectednegidx <- sample(negidx, size = min(lambda, length(negidx)), replace = FALSE)
  }
  
  selecedidx = c(selectedposidx, selectednegidx)
  
  return(selecedidx)
  
  
}

selfpaced.single <- function(splitid, lambda, uplambda, Iter_num) {
  
  trainidx <- which(splits[,splitid]==0)
  testidx <- which(splits[,splitid]==1)
  X_train <- X[trainidx,]		
  Y_train <- Y[trainidx]
  X_test <- X[testidx,]
  Y_test <- Y[testidx]
  
  #the list storing the result for each iteration
  valpredmatrix <- list()
  evlpredmatrix <- list()
  coefmatrix <- list()
  nonzerocoefmatrix <- list()
  valmaps <- replicate(Iter_num,0)
  evlmaps <- replicate(Iter_num,0)
  iters.vidx = list()	
  
  #the starting values
  cvfit <- cv.glmnet(x = X_train,
                     y = Y_train,
                     alpha = 1,
                     family = "binomial",
                     type.measure = "class")
  
  valpred <- predict(cvfit, X_train, type = "response", s = "lambda.min")
  
  v.idx <- selfpace1.rank(dev_decval = valpred,
                          dev_labels = Y_train, 
                          lambda = lambda)
  this.training.vidx <- v.idx
  iters.vidx <- list()	
  
  for(iter in 1:Iter_num) {
    if(length(this.training.vidx) == length(Y_train)){break}
    cat("Starting the ",iter,"-th iteration.\t", sep = "")
    
    iters.vidx[[iter]] <- this.training.vidx
    # glmnet (Lasso, Elastic net & L2)
    cvfit<-cv.glmnet(x = X_train[this.training.vidx,],
                     y = Y_train[this.training.vidx],
                     alpha = 1,
                     family = "binomial",
                     type.measure = "class")  
    valprediction <- predict(cvfit, X_train, type = "response", s = "lambda.min")
    tsprediction <- predict(cvfit, X_test, type = "response", s = "lambda.min")
    coefprediction <- as.vector(coef(cvfit, s = "lambda.min")[-1])
    numbernonzerocoef <- length(which(coefprediction!=0))
    
    #self-paced learning
    selectedidx <- selfpace1.rank(dev_decval = valprediction,
                                  dev_labels = Y_train, 
                                  lambda = lambda)
    this.training.vidx <- selectedidx
    cat("Select ", length(selectedidx), " samples.\t", sep = "")
    
    #change the parameter accoding to the step size
    lambda <- lambda + uplambda
    
    #store the prediction for this iteration
    valpredmatrix[[iter]] <- as.numeric(valprediction)
    evlpredmatrix[[iter]] <- tsprediction
    coefmatrix[[iter]] <- coefprediction
    nonzerocoefmatrix[[iter]] <- numbernonzerocoef
    iters.vidx[[iter]] = selectedidx
    
    #evaluate the training and test error
    val_loss <- sum((valpredmatrix[[iter]] - Y_train)^2)
    evl_loss <- sum((evlpredmatrix[[iter]] - Y_test)^2)
    valmaps[iter] <- val_loss
    evlmaps[iter] <- evl_loss
    
    cat("Finish the ",iter,"-th iteration.\n", sep = "")
  }
  
  results <- list("valpredmatrix" = valpredmatrix, 
                  "evlpredmatrix" = evlpredmatrix,
                  "valmaps" = valmaps,
                  "evlmaps" = evlmaps,
                  "itervidx" = iters.vidx, 
                  "Coef"=coefmatrix, 
                  "Coef.num"=nonzerocoefmatrix)
  
  return(results)
}

predict.label <- function(predictlabels,truelabels,View_num){
  # v in each iter
  label_num = length(truelabels)
  yk_label1 = 1
  yk_label0 = 0
  
  pre_label_final = vector(mode = "numeric", length = label_num)
  
  for(i in 1:label_num){
    pre_label = predictlabels[i,]
    loss1 = 0
    loss0 = 0
    for(j in 1:View_num){
      loss1 = loss1 + (pre_label[j] - yk_label1)^2
      loss0 = loss0 + (pre_label[j] - yk_label0)^2
    }
    if(loss1>=loss0){
      pre_label_final[i] = 0
    }else{
      pre_label_final[i] = 1
    }
  }
  return(pre_label_final)
}


MVSPL.rank <- function(dev_decval, dev_labels, v_iter, View_id, lambda, gamma, View_num, num_add) {
  
  
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
  
  ## sort sample
  selectedposidx = list()
  selectednegidx = list()
  selecedidx = list()
  pos_lambda = matrix(0, nrow = 1, ncol = View_num)
  neg_lambda = matrix(0, nrow = 1, ncol = View_num)
  V_iter = matrix(0, nrow = length(dev_labels), ncol = View_num)
  posidx = which(dev_labels==1)	#postive id mapping
  negidx = which(dev_labels==0)	#negative id mapping
  
  #calculate pos_lambda and neg_lambda according to the rank
  for(i in 1:View_num){
    
    if(length(which(v_iter[posidx,i]==1))!=0){  # a certain class has samples
      pos_lambda[i] = sort(loss[posidx,i][which(v_iter[posidx,i]==1)])[min(length(posidx), num_add, length(which(v_iter[posidx,i]==1)))]
    }
    if(length(which(v_iter[negidx,i]==1))!=0){
      neg_lambda[i] = sort(loss[negidx,i][which(v_iter[negidx,i]==1)])[min(length(negidx), num_add, length(which(v_iter[negidx,i]==1)))]
    }
    
    if(length(unique(loss[posidx,i]))!=1){  # select samples
      selectedposidx[[i]] <- intersect(posidx[which(v_iter[posidx,i] == 1)],   ## v_iter = 1 && loss is small
                                       posidx[which(loss[posidx,i] <= pos_lambda[i])])
    }else{
      selectedposidx[[i]] <- sample(posidx, size = min(num_add, length(posidx)), replace = FALSE)
    }
    
    if(length(unique(loss[negidx,i]))!=1){
      selectednegidx[[i]] <- intersect(negidx[which(v_iter[negidx,i] == 1)],
                                       negidx[which(loss[negidx,i] <= neg_lambda[i])])
    }else{
      selectednegidx[[i]] <- sample(negidx, size = min(num_add, length(negidx)), replace = FALSE)
    }
    
    
    selecedidx[[i]] = c(selectedposidx[[i]], selectednegidx[[i]])
    V_iter[selecedidx[[i]],i] = 1
  }
  cat("The ",View_id, "-th view select ",length(selecedidx[[i]])," samples.\t", sep = "")
  #return the result
  return(V_iter)
  
  
  
  ## Original code
  # # initialize the loss function
  # loss = matrix(0, nrow = length(dev_labels), ncol = View_num)
  # 
  # #calculate the loss
  # for(m in 1: View_num){
  #   if(m != View_id){
  #     loss[,m] = (dev_decval[,m] - dev_labels)^2    #squared error
  #   }else{
  #     next;
  #   }
  # }
  # 
  # # Update v(View_num-j)
  # for(m in 1:View_num){
  #   if(m != View_id){
  #     for(i in 1:length(dev_labels)){
  #       if(loss[i,m] < lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m])){
  #         v_iter[i,m] = 1
  #       }else{
  #         v_iter[i,m] = 0
  #       }
  #     }
  #   }
  # }
  # 
  # # Update vj
  # loss[,View_id] = (dev_decval[,View_id] - dev_labels)^2
  # for(i in 1:length(dev_labels)){
  #   if(loss[i,View_id] < lambda[View_id] + gamma * (sum(v_iter[i,])-v_iter[i,View_id])){
  #     v_iter[i,View_id] = 1
  #   }else{
  #     v_iter[i,View_id] = 0
  #   }
  # }
  # #return the result
  # return(v_iter)
}
