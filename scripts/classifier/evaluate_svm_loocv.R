# Evaluate Support Vector Machine classifier by leave-one-out cross-validation.
#
# Set working directory.
#setwd("~")
# Set seed value.
set.seed(1)
# Load libraries. ----------------------------------------------------------
source("auxiliary.R")
source("algorithms/SupportVectorMachine.R")

# Load *.fasta dataset. ---------------------------------------------------
# Load training dataset. k is the numebr of k-mer.
train.dataset.org <- loadDataset(K = 3)
train.datalist <- train.dataset.org$datalist

# Adjust the training dataset. A sample size for datasets except for nrEVE is set as 950.
train.dataset.fit <- fitSampleSize(train.dataset.org, 950)
train.data <- train.dataset.fit$data
train.labels <- as.numeric(train.dataset.fit$labels)
class.num <- length(unique(train.labels))

# Standardize training dataset.
train.data.mean <- apply(train.data, 2, mean)
train.data.sd <- apply(train.data, 2, sd)
train.data <- sweep(sweep(train.data, 2, train.data.mean), 2, train.data.sd, FUN="/")

# Analysis of Support Vector Machine --------------------------------------
# Tuning hyperparameters list.
params <- list("gamma" = 10^(-3:3),
               "cost" = 10^(-2:2),
               "cross" = 2
)

# Leave-one-out cross-validation.
cat("# Leave-one-out cross-validation...\n")
fit <- NULL
for(itr in 1:length(train.labels)){
  x_train <- train.data[-itr, ]
  x_test <- train.data[itr, ]
  y_train <- train.labels[-itr]
  y_test <- train.labels[itr]
  
  tune.svm <- tune.SupportVectorMachine(x_train, y_train, params)
  best.params <- list("gamma" = tune.svm$gamma,
                      "cost" = tune.svm$cost)
  
  # Construct SVM classifier, and predict test dataset.
  fit[[itr]] <- fit.SupportVectorMachine(x_train, y_train, x_test, y_test, best.params)
}

# Construct a confusion matrix of the classifier.
ConfusionMatrix <- matrix(0, class.num, class.num)
for(itr in 1:length(train.labels)){
  ConfusionMatrix <- ConfusionMatrix + fit[[itr]]$ConfusionMatrix
}

# Calculate precision and recall.
svm_pr <- calc_precision_recall(ConfusionMatrix)

# Summary of the results.
cat("--- Summary ---\n")
cat("Analysis Method : Support Vector Machine (leave-one-out cross-validation)\n")
cat("Precision:",svm_pr$precision, "\n")
cat("Recall:",svm_pr$recall, "\n")
cat("Confusion Matrix:\n")
print(ConfusionMatrix)
