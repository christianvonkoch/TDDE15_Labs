## 1. Show that multiple runs of the hill-climbing algorithm can return non-equivalent Bayesian network (BN) structures. 
## Explain why this happens. Use the Asia dataset which is included in the bnlearn package. To load the data, 
## run data("asia").

library(bnlearn)
library(gRain)
library(RBGL)
library(Rgraphviz)
data("asia")
init_dag1 = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
init_dag2 = model2network("[B][D][X][E|D:B][L|B][T|L:B:X][S|X][A|B]")
plot(init_dag1)
plot(init_dag2)
random_restarts=c(1,10,100)
scores=c("aic", "bic", "bde")
iss=c(1,10,100)

# Testing different initial structures
model1_init=hc(asia, start=init_dag1)
model2_init=hc(asia, start=init_dag2)
plot(model1_init)
plot(model2_init)
all.equal(model1_init, model2_init)

# Changed the resulting network

# Testing different random restarts
model1_random=hc(asia, restart=random_restarts[1])
model2_random=hc(asia, restart=random_restarts[2])
model3_random=hc(asia, restart=random_restarts[3])
plot(model1_random)
plot(model2_random)
plot(model3_random)
all.equal(model1_random, model2_random)
all.equal(model1_random, model3_random)
all.equal(model2_random, model3_random)

# Changed the resulting network

# Testing different score models
model1_score=hc(asia, score=scores[1])
model2_score=hc(asia, score=scores[2])
model3_score=hc(asia, score=scores[3])
plot(model1_score)
plot(model2_score)
plot(model3_score)
all.equal(model1_score, model2_score)
all.equal(model1_score, model3_score)
all.equal(model2_score, model3_score)

# Fundamentally changed the resulting network

# Testing different imaginary sample sizes
model1_iss=hc(asia, score="bde", iss=iss[1])
model2_iss=hc(asia, score="bde", iss=iss[2])
model3_iss=hc(asia, score="bde", iss=iss[3])
plot(model1_iss)
plot(model2_iss)
plot(model3_iss)
all.equal(model1_iss, model2_iss)
all.equal(model1_iss, model3_iss)
all.equal(model2_iss, model3_iss)

# Fundamentally changed the resulting network

## 2. Learn a BN from 80 % of the Asia dataset. The dataset is included in the bnlearn
## package. To load the data, run data("asia"). Learn both the structure and the
## parameters. Use any learning algorithm and settings that you consider appropriate.
## Use the BN learned to classify the remaining 20 % of the Asia dataset in two classes:
##  S = yes and S = no. In other words, compute the posterior probability distribution of S
## for each case and classify it in the most likely class. To do so, you have to use exact
## or approximate inference with the help of the bnlearn and gRain packages, i.e. you
## are not allowed to use functions such as predict. Report the confusion matrix, i.e.
## true/false positives/negatives. Compare your results with those of the true Asia BN,
## which can be obtained by running
## dag = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]").

set.seed(12345)
data("asia")
n=dim(asia)[1]
id=sample(1:n, floor(n*0.8))
train=asia[id,]
test=asia[-id,]

# Function for calculating marginal likelihoods for different target variables.
# Inputs:
# juncTree - The function tree obtained from the Laurent Spielberger algorithm
# data - The data which will be used as a basis for the prediction
# features - The features that should be the basis for the prediction
# target - The target variable for which the prediction should be made
predictNet <- function(juncTree, data, features, target){
  predArray <- matrix(nrow=nrow(data),ncol=1)
  for(i in 1:nrow(data)){
    obsStates <- NULL
    for(p in features){
      if(data[i,p]=="yes"){
        obsStates <- c(obsStates,"yes")
      } else{
        obsStates <- c(obsStates,"no")
      }
    }
    
    
    obsEvidence <- setEvidence(object = juncTree,
                               nodes = features,
                               states = obsStates)
    obsPredProb <- querygrain(object = obsEvidence,
                              nodes = target)$S
    predArray[i] <- if(obsPredProb["yes"]>=0.5) "yes" else "no"
  }
  return(predArray)
}

# Train model
BNmodel=hc(train)
real_dag = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
fit=bn.fit(BNmodel, train)
fit_real=bn.fit(real_dag, train)
fitTable=as.grain(fit)
fitTable_real=as.grain(fit_real)
junctionTree=compile(fitTable)
junctionTree_real=compile(fitTable_real)

obsVars=c("A", "T", "L", "B", "E", "X", "D")
tarVar=c("S")
predictionTest <- predictNet(junctionTree, test, obsVars, tarVar)
predictionTrue <- predictNet(junctionTree_real, test, obsVars, tarVar)
confTable=table(predictionTest, test$S)
confTable_real=table(predictionTrue, test$S)
confTable
confTable_real

## 3, In the previous exercise, you classified the variable S given observations for all the
## rest of the variables. Now, you are asked to classify S given observations only for the
## so-called Markov blanket of S, i.e. its parents plus its children plus the parents of its
## children minus S itself. Report again the confusion matrix.

markovBlanket=mb(fit_real, node="S")
print(markovBlanket)
predictMarkovBlanket=predictNet(junctionTree, test, markovBlanket, tarVar)
confTable_MB=table(predictMarkovBlanket, test$S)
confTable_MB

## 4. Repeat the exercise (2) using a naive Bayes classifier, i.e. the predictive variables are
## independent given the class variable. See p. 380 in Bishop's book or Wikipedia for
## more information on the naive Bayes classifier. Model the naive Bayes classifier as a
## BN. You have to create the BN by hand, i.e. you are not allowed to use the function
  ## naive.bayes from the bnlearn package

set.seed(12345)
naive_dag=model2network("[S][A|S][B|S][X|S][T|S][L|S][E|S][D|S]")
naive_fit=bn.fit(naive_dag, train)
fitNaive=as.grain(naive_fit)
naive_JunctionTree=compile(fitNaive)
naivePred=predictNet(naive_JunctionTree, test, obsVars, tarVar)
naiveConf=table(naivePred, test$S)
naiveConf

## 5. Explain why you obtain the same or different results in the exercises (2-4).

## Conclusion: Markov Blanket (MB) consists of the nodes of importance for the target nodes' dependencies in the network.
## If we do not change the MB, the outcome of the prediction on the target node will not change.
