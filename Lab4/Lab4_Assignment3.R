# Fetching data
library(AtmRay)

data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/
GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111)
SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE)
train=data[SelectTraining,]
test=data[-SelectTraining,]

## 1) Use the R package kernlab to fit a Gaussian process classification model for fraud on the training data. Use the default kernel
## and hyperparameters. Start using only the covariates varWave and skewWave in the model. Plot contours of the prediction probs
## over a suitable grid of values of varWave and skewWave. Overlay the training data for fraud=1 (as blue points) and fraud=0
## (as red points). You can reuse code from the file KernLabDemo.R available on course website. Compute the confusion matrix
## for the classifier and its accuracy. 

model = gausspr(fraud ~ varWave+skewWave, data=train)
predictedTrain = predict(model, newdata=train)
predictedTrain
confusionMatrix = table(predictedTrain, train[,5]) # confusion matrix

# class probabilities 
probPreds <- predict(model, train, type="probabilities")
x1 <- seq(min(train$varWave),max(train$varWave),length=100)
x2 <- seq(min(train$skewWave),max(train$skewWave),length=100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))

gridPoints <- data.frame(gridPoints)
names(gridPoints) <- names(subset(train, select=c("varWave", "skewWave")))
probPreds <- predict(model, gridPoints, type="probabilities")

# Plotting for Prob(Fraud)
contour(x1,x2,matrix(probPreds[,1],100,byrow = TRUE), 20, xlab = "varWave", ylab = "skewWave", main = 'Prob(Fraud) - Fraud is red')
points(train[train[,5]==1,1],train[train[,5]==1,2],col="red")
points(train[train[,5]==0,1],train[train[,5]==0,2],col="blue")

# Confusion matrix and accuracy
confusionMatrix
sum(diag(confusionMatrix))/sum(confusionMatrix)

## 2) Using the estimated model from 1), make predictions on the test set.

predictedTest = predict(model, newdata=test)
confusionMatrix_test = table(predictedTest, test[,5])
sum(diag(confusionMatrix_test))/sum(confusionMatrix_test)

## 3) Train a model using all four covariates. Make predictions on the test set and compare the accuracy to the model with only two
## covariates

model2 = gausspr(fraud ~., data=train)
predictedTest2 = predict(model2, newdata=test)
confusionMatrix_test2 = table(predictedTest2, test[,5])
sum(diag(confusionMatrix_test2))/sum(confusionMatrix_test2)
