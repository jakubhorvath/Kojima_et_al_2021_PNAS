# A wrapper of multiclass logistic classifier.
# Load libraries. ----------------------------------------------------------
library(glmnet)

# Tune hyperparameters.
cv.multiclassLogistic <-function(data, labels, family, alpha){
  cv.fit <- cv.glmnet(data, labels, family=family, alpha=alpha)
  return(cv.fit)
}

# Multiclass logistic classifier.
fit.multiclassLogistic <- function(data, labels, family, alpha, lambda){
  fit <- glmnet(data, labels, family=family, alpha=alpha, lambda=lambda)
  xxx <- predict(fit, data)[,,1]
  predict_p <- sweep( exp(xxx), 1, apply(exp(xxx), 1, sum), FUN="/" )
  
  predict_P <- rep(0, nrow(predict_p))
  for(i in 1:length(predict_P)) predict_P[i] <- which.max(predict_p[i, ])
  
  ConfusionMatrix <- matrix(0, ncol(predict_p), ncol(predict_p))
  for(i in 1:length(predict_P))
  {
    ii <- predict_P[i]
    jj <- as.numeric(labels[i])
    ConfusionMatrix[ii, jj] <- ConfusionMatrix[ii, jj] + 1 
  }
  accuracy <- mean(predict_P == labels)
  return(list(fit=fit, ConfusionMatrix=ConfusionMatrix, accuracy=accuracy))
}
