# Fetching data

data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/
GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111)
SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE)

## 1) Use the R package kernlab to fit a Gaussian process classification model for fraud on the training data. Use the default kernel
## and hyperparameters. Start using only the covariates varWave and skewWave in the model. Plot contours of the prediction probs
## over a suitable grid of values of varWave and skewWave. Overlay the training data for fraud=1 (as blue points) and fraud=0
## (as red points). You can reuse code from the file KernLabDemo.R available on course website. Compute the confusion matrix
## for the classifier and its accuracy. 

