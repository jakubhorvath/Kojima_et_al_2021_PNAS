# Predict Class via Support Vector Machine --------------------------------
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

# Load test dataset., k the numebr of k of k-mer.
test.dataset.org <- loadTestDataset(K = 3)
test.data <- test.dataset.org$data
test.labels <- as.numeric(test.dataset.org$labels)
test.id <- test.dataset.org$id

# Standardize test dataset. on the mean and variance of training data.
test.data <- sweep(sweep(test.data, 2, train.data.mean), 2, train.data.sd, FUN="/")

# Analysis of Support Vector Machine --------------------------------------
# Tune hyperparameters.
params <- list("gamma" = 10^(-3:3),
               "cost" = 10^(-2:2),
               "cross" = 2)
cat("# Tune hyperparameters...\n")
tune.svm <- tune.SupportVectorMachine(train.data, train.labels, params)

# Set the optimal hyperparameters.
params <- list("gamma" = tune.svm$gamma,
               "cost" = tune.svm$cost)

# Construct SVM classifier, and predict test dataset.
cat("# Construct the SVM classifier, and predict test dataset...\n")
fit <- fit.SupportVectorMachine(train.data, train.labels, test.data, test.labels, params)

# Calculate precision and recall.
svm_pr <- calc_precision_recall(fit$mat)

# Summary of the results.
cat("# --- Summary ---\n")
cat("# Analysis Method : Support Vector Machine\n")
cat("# Precision:",svm_pr$precision, "\n")
cat("# Recall:",svm_pr$recall, "\n")
cat("# Confusion Matrix:\n")
print(fit$mat)

# Display prediction classes for test dataset.
if(FALSE){
for(i in 1:nrow(test.data)){
  cat(test.id[i])
  cat("\n")
  switch(as.character(fit$predclass[i]),
         "1"=cat("EVE"),
         "2"=cat("Coding"),
         "3"=cat("InterGenic"),
         "4"=cat("Intron"),
         "5"=cat("Noncoding"),
         "6"=cat("Pseudogene"),
         "7"=cat("TSS")
  )
  cat("\n")
}
}
