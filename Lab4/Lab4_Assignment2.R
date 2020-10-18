## GP regression with kernlab

# Fetch data from source

data = read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/
Code/TempTullinge.csv", header=TRUE, sep=";")

# Restructuring data

data$date = as.Date(data$date, "%d/%m/%y")
time=seq(1:2190)
day=time %% 365
day[which(day == 0)] = 365
id = seq(from=1, to=2186, 5)
time=time[id]
day=day[id]

## Assignment 1: Familirarization of the kernlab library. Calculating covariance matrix from two input vectors.

# This is just to test how one evaluates a kernel function
# and how one computes the covariance matrix from a kernel function.
library(kernlab)
X <- as.vector(c(1,3,4)) # Simulating some data.
Xstar <- as.vector(c(2,3,4))
ell <- 1

# Covariance function
SquaredExpKernel <- function(x1,x2,sigmaF=1,l=3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

# This function is a nested function which returns an object of class kernel. 
NestedSquaredExpKernel <- function(sigmaF=1,l=3){
  EvaluExpKernel = function(x, xStar) {
    n1 <- length(x)
    n2 <- length(xStar)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigmaF^2*exp(-0.5*( (x-xStar[i])/l)^2 )
    }
    return(K)
  }
  class(EvaluExpKernel)='kernel'
  return(EvaluExpKernel)
}

SEkernel <- NestedSquaredExpKernel() # Note how I reparametrize the rbfdot (which is the SE kernel) in kernlab.
SEkernel(1,2) # Just a test - evaluating the kernel in the points x=1 and x'=2.
# Computing the whole covariance matrix K from the kernel. Just a test.
kernelMatrix(kernel = SEkernel, x = X, y = Xstar) # So this is K(X,Xstar).

## Assignment 2: Consider model: temp=f(time)+epsilon with epsilon~N(0, sigmaN^2). f~GP(0, k(time, time'))
## Let sigmaN be the residual variance from a simple quadratic regression fit. Estimate the above gaussian process regression
## model using the squared exponentional function from 1) with sigmaF=20 and l=0.2. Use the predict function n R to compute the
## posterior mean at every data point in the training dataset. Make a scatterplot of the data and superimpose the posterior mean
## of f as a curve (use type="l" in the plot function). Play around with different values on sigmaF and l.

temp=data$temp[id]
fit = lm(temp ~ time+time^2)
sigmaN=sd(fit$residuals)

sigmaF=20
l=0.2
SEkernel = NestedSquaredExpKernel(sigmaF, l)
model = gausspr(time, temp, kernel=SEkernel, var=sigmaN^2)
predictedMean=predict(model, newdata=time)
plot(time, temp, type="p", main="Time vs temperature")
lines(time, predictedMean, type="l", lwd=2, xlab="Time", ylab="Temp", col="red")
legend("bottomright", legend=c("Data", "Predicted mean"), pch=c(1, NA), lty=c(NA, 1), lwd=c(NA, 2), col=c("black", "red"))

## Assignment 3: Kernlab can compute posterior variance of f, but it seems to be a bug in the code. So, do your own computations of the
## posterior variance of f and plot the 95 % probability (pointwise) bands for f. Superimpose these bands on a figure with the
## posterior mean that you obtained in (2)

# Function for returning the posterior mean and variance of a given dataset (gaussian process regression)
posteriorGP = function(X, y, XStar, sigmaNoise, k, ...) {
  n = length(X)
  K=k(X, X, ...)
  kStar=k(X,XStar, ...)
  L = t(chol(K + sigmaNoise^2*diag(n)))
  alpha=solve(t(L),solve(L, y))
  predMean=t(kStar)%*%alpha
  v=solve(L, kStar)
  predVar=k(XStar, XStar, ...)-t(v)%*%v
  return(list(mean=predMean, var=predVar))
}

# Storing mean and var so that the data can be rescaled to original size at later stage
scale_mean=mean(temp)
scale_var=sqrt(var(temp))

posterior = posteriorGP(scale(time), scale(temp), scale(time), sigmaN, SquaredExpKernel, sigmaF, l)
postVar = posterior$var
postMean = posterior$mean
plot(time, temp, type="p", main="Time vs temperature")
lines(time, predictedMean, type="l", lwd=2, xlab="Time", ylab="Temp", col="red")
lines(time, postMean*scale_var+scale_mean+1.96*sqrt(diag(postVar)), lwd=2, lty=2, col="blue")
lines(time, postMean*scale_var+scale_mean-1.96*sqrt(diag(postVar)), lwd=2, lty=2, col="blue")
legend("bottomright", legend=c("Data", "Predicted mean", "95% probability bands"), pch=c(1, NA, NA), lty=c(NA, 1, 2),
       lwd=c(NA, 2, 2), col=c("black", "red", "blue"))

## Assignment 4: Consider now the following model: temp = f(day) + epsilon with epsilon ~ N(0, sigmaN^2) and f~GP(0, k(day, day'))
## Estimate the model using the squared exponential function with sigmaF=20 and l=0.2. Superimpose the posterior mean
## from this model on the posterior mean from the model in (2). Note that this plot should also have time variables on the 
## horizontal axis. Compare the results of both models. What are the pros and cons of each model?

sigmaF=20
l=0.2
SEkernel = NestedSquaredExpKernel(sigmaF, l)
model2 = gausspr(day, temp, kernel=SEkernel, var=sigmaN^2)
predictedMean2=predict(model2, newdata=day)
plot(time, temp, type="p", main="Time vs temperature")
lines(time, predictedMean, type="l", lwd=2, xlab="Time", ylab="Temp", col="red")
lines(time, predictedMean2, type="l", lwd=2, col="blue")
legend("bottomright", legend=c("Data", "Predicted mean time", "Predicted mean day"), pch=c(1, NA, NA), lty=c(NA, 1, 1),
       lwd=c(NA, 2, 2), col=c("black", "red", "blue"))

## Assignment 5: Finally, implement a generalization of the periodic kernel given in the lectures. Note that Note that we have two different
## length scales here, and `2 controls the correlation between the same day in different years. Estimate the GP model using the
## time variable with this kernel and hyperparameters sigmaF = 20, l1 = 1, l2 = 10 and d = 365/sd(time).
## The reason for the rather strange period here is that kernlab standardizes the inputs to have standard deviation of 1. 
## Compare the fit to the previous two models (with sigmaF = 20 and l = 0.2). Discuss the results.

PeriodicKernel <- function(sigmaF=1,l1=1, l2=10, d){
  val = function(x, xStar) {
    diff = abs(x-xStar)
    return(sigmaF^2*exp(-((2*sin(pi*diff))/d)/l1^2)*exp(-0.5*diff^2/l2^2))
  }
  class(val)='kernel'
  return(val)
}

sigmaF=20
l1=1
l2=10
d=365/sd(time)

PerKernel = PeriodicKernel(sigmaF, l1, l2, d)
model3 = gausspr(time, temp, kernel=PerKernel, var=sigmaN^2)
predictedMean3=predict(model3, newdata=time)
lines(time, predictedMean3, type="l", lwd=2, col="green")
legend("bottomright", legend=c("Data", "Predicted mean time", "Predicted mean day", "Periodic kernel"), 
       pch=c(1, NA, NA, NA), lty=c(NA, 1, 1, 1),
       lwd=c(NA, 2, 2, 2), col=c("black", "red", "blue", "green"))
