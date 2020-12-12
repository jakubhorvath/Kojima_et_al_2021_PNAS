# A wrapper of support vector machine.
# Load libraries. ----------------------------------------------------------
library(e1071)

# Tune hyperparameters.
tune.SupportVectorMachine <- function(data, labels, params){
  model <- tune.svm(as.factor(labels)~., data=data.frame(data, labels), gamma=params$gamma, cost=params$cost,
                    tunecontrol = tune.control(sampling = "cross", cross = params$cross))
  return(model$best.model)
}

# Support Vector Machine
fit.SupportVectorMachine <- function(train.data, train.labels, test.data, test.labels, params){
  # Construct SVM classifier.
  fit <- svm(as.factor(train.labels)~., data=data.frame(train.data, train.labels), gamma=params$gamma, cost=params$cost)
  
  # Predict class of test dataset.
  if(nrow(as.matrix(t(test.data))) == 1){
    class.num <- length(unique(train.labels))
    ConfusionMatrix <- matrix(0, class.num, class.num)
    
    pred <- predict(fit, t(test.data))
    ii <- pred
    jj <- as.numeric(test.labels)
    ConfusionMatrix[ii, jj] <- ConfusionMatrix[ii, jj] + 1
    accuracy <- NULL
    mat <- NULL
  }else{
    pred <- predict(fit, test.data)
    # Make confusion matrix
    mat <- table(pred, test.labels)
    class.num <- length(unique(test.labels))
    ConfusionMatrix <- matrix(as.numeric(mat), class.num, class.num)
    
    accuracy <- mean(pred == test.labels)
  }
  return(list(ConfusionMatrix=ConfusionMatrix, accuracy=accuracy, predclass=pred, mat=mat))
}
