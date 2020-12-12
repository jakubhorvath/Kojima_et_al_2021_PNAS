# calculate contibution rate of each k-mer by multiclass logistic classifier
#
# Set working directory.
#setwd("~")
# Set seed value.
set.seed(1)
# Load libraries. ----------------------------------------------------------
source("auxiliary.R")
source("algorithms/MulticlassLogistic.R")

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


# Analysis of multiclass logistic classifier ------------------------------
feature.num <- ncol(train.data)
# Initialize the histogram of selected features.
hist.selected.features = array(rep(0, (class.num * ncol(train.data))), dim=c(class.num, ncol(train.data)))

# Calculate the selection probability of k-mer features by the bootstrap method.
for(b in 1:1000){
  cat(b," ")
  # Resampling.
  resampeld.itr = sample(1:nrow(train.data), nrow(train.data), replace=TRUE)
  resampled.data = train.data[resampeld.itr, ]
  resampled.label = train.labels[resampeld.itr]
  while(length(unique(resampled.label))<length(unique(train.labels))){
    resampeld.itr = sample(1:nrow(train.data), nrow(train.data), replace=TRUE)
    resampled.data = train.data[resampeld.itr, ]
    resampled.label = train.labels[resampeld.itr]
  }
  # Standardize resampled dataset.
  resampled.data.mean <- apply(resampled.data, 2, mean)
  resampled.data.sd <- apply(resampled.data, 2, sd)
  resampled.data<- sweep(sweep(resampled.data, 2, resampled.data.mean), 2, resampled.data.sd, FUN="/")
  
  # Tune hyperparameters.
  cv.fit <- cv.multiclassLogistic(resampled.data, resampled.label, family="multinomial", alpha=1)
  
  tryCatch({
    # Fitting.
    fit.mclogistic <- fit.multiclassLogistic(resampled.data, resampled.label, family="multinomial", alpha=0.9, lambda=cv.fit$lambda.min)
    
    # Store numeber of selected features for each class.
    for(class in  1:class.num){
      # selected.features is flag vector of selected features.
      selected.features <- numeric(feature.num)
      selected.features[fit.mclogistic[["fit"]][["beta"]][[class]]@i + 1] <- selected.features[fit.mclogistic[["fit"]][["beta"]][[class]]@i + 1] + 1
      hist.selected.features[class,] = hist.selected.features[class,] + selected.features
    }
  },
  error = function(e){
    message("Error! ")
    b <- b-1
  },
  warning = function(e){
    message("Warning! ") 
    b <- b-1
  }
  )
}# end of bootstrap.

### output of contribution rates
K_mer.names <- K_mer_name_list(K=3)
hist.selected.features = hist.selected.features/b
colnames(hist.selected.features) <- K_mer.names
write.table(hist.selected.features, "contribution_rates.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
