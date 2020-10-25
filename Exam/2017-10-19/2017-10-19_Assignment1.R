## Learn a Bayesian network (BN) from the Asia dataset that is included in the bnlearn
## package. To load the data, run data("asia"). Learn both the structure and the
## parameters. Use any learning algorithm and settings that you consider appropriate.
## Identify a d-separation in the BN learned and show that it indeed corresponds to an
## independence in the probability distribution represented by the BN. To do so, you may
## want to use exact or approximate inference with the help of the bnlearn and gRain
## packages. (2.5 p)

library(bnlearn)
library(gRain)
library(RBGL)
library(Rgraphviz)
data("asia")

BNmodel=hc(asia)
fit=bn.fit(BNmodel, asia)
fitTable=as.grain(fit)
fitTable_real=as.grain(fit_real)
junctionTree=compile(fitTable)