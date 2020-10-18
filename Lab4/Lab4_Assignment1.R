library(mvtnorm)

## Assignment 1: Implement algorithm 2.1. Write your own code for simulating from the posterior distribution of f using the squared
## exponential kernel (f is distributed according to a gaussian process regression model). 

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

posteriorGP = function(X, y, XStar, sigmaNoise, k, ...) {
  n=length(X)
  K=k(X, X, ...)
  kStar=k(X,XStar, ...)
  L = t(chol(K + sigmaNoise^2*diag(n)))
  alpha=solve(t(L),solve(L, y))
  predMean=t(kStar)%*%alpha
  v=solve(L, kStar)
  predVar=k(XStar, XStar, ...)-t(v)%*%v
  return(list(mean=predMean, var=predVar))
}

# Assignment 2: Plot one draw from posterior with observation (x,y)=(0.4,0.719).

xTest = seq(-1, 1, length=100)

# Plotting one draw
obs=data.frame(x=0.4, y=0.719)
sigmaF <- 1
sigmaN=0.1
l <- 0.3
posteriorSim=posteriorGP(obs$x, obs$y, xTest, sigmaN, SquaredExpKernel, sigmaF, l)
plot(xTest, posteriorSim$mean, type="l", 
     ylim=c(min(posteriorSim$mean - 1.96*sqrt(diag(posteriorSim$var))), max(posteriorSim$mean + 1.96*sqrt(diag(posteriorSim$var)))),
     col="red", main="Plot of posterior mean (red) with 95 % probability bands", xlab="x", ylab="Posterior mean",
     sub="One observation as prior")
lines(xTest, posteriorSim$mean - 1.96*sqrt(diag(posteriorSim$var)), col = "gray", lwd = 2, lty=21)
lines(xTest, posteriorSim$mean + 1.96*sqrt(diag(posteriorSim$var)), col = "gray", lwd = 2, lty=21)

## Assignment 3: Now update posterior with two observations. 

x=c(0.4, -0.6)
y=c(0.719, -0.044)
obs2=data.frame(x=x, y=y)
posteriorSim2=posteriorGP(obs2$x, obs2$y, xTest, sigmaN, SquaredExpKernel, sigmaF, l)
plot(xTest, posteriorSim2$mean, type="l", 
     ylim=c(min(posteriorSim2$mean-1.96*sqrt(diag(posteriorSim2$var))), max(posteriorSim2$mean+1.96*sqrt(diag(posteriorSim2$var)))),
     col="red", main="Plot of posterior mean (red) with 95 % probability bands", xlab="x", ylab="Posterior mean",
     sub="Two observations as prior")
lines(xTest, posteriorSim2$mean - 1.96*sqrt(diag(posteriorSim2$var)), col = "gray", lwd = 2, lty=21)
lines(xTest, posteriorSim2$mean + 1.96*sqrt(diag(posteriorSim2$var)), col = "gray", lwd = 2, lty=21)

## Assignment 4: Now use 5 observations and plot the posterior.

x=c(-1, -0.6, -0.2, 0.4, 0.8)
y=c(0.768, -0.044, -0.940, 0.719, -0.664)
obs2=data.frame(x=x, y=y)
posteriorSim2=posteriorGP(obs2$x, obs2$y, xTest, sigmaN, SquaredExpKernel, sigmaF, l)
plot(xTest, posteriorSim2$mean, type="l", 
     ylim=c(min(posteriorSim2$mean-1.96*sqrt(diag(posteriorSim2$var))), max(posteriorSim2$mean+1.96*sqrt(diag(posteriorSim2$var)))),
     col="red", main="Plot of posterior mean (red) with 95 % probability bands", xlab="x", ylab="Posterior mean",
     sub="Five observations as prior")
lines(xTest, posteriorSim2$mean - 1.96*sqrt(diag(posteriorSim2$var)), col = "gray", lwd = 2, lty=21)
lines(xTest, posteriorSim2$mean + 1.96*sqrt(diag(posteriorSim2$var)), col = "gray", lwd = 2, lty=21)

## Assignment 5. Repeat 4 with hyperparam sigmaF=1 and l=1. Compare the results.

sigmaF=1
l=1
posteriorSim2=posteriorGP(obs2$x, obs2$y, xTest, sigmaN, SquaredExpKernel, sigmaF, l)
plot(xTest, posteriorSim2$mean, type="l", 
     ylim=c(min(posteriorSim2$mean-1.96*sqrt(diag(posteriorSim2$var))), max(posteriorSim2$mean+1.96*sqrt(diag(posteriorSim2$var)))),
     col="red", main="Plot of posterior mean (red) with 95 % probability bands", xlab="x", ylab="Posterior mean",
     sub=expression(paste("Five observations as prior and ", sigma[F], "=1 and l=1")))
lines(xTest, posteriorSim2$mean - 1.96*sqrt(diag(posteriorSim2$var)), col = "gray", lwd = 2, lty=21)
lines(xTest, posteriorSim2$mean + 1.96*sqrt(diag(posteriorSim2$var)), col = "gray", lwd = 2, lty=21)
