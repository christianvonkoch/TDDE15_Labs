##  Show that multiple runs of the hill-climbing algorithm can return non-equivalent Bayesian network (BN) structures. 
## Explain why this happens. Use the Asia dataset which is included in the bnlearn package. To load the data, 
## run data("asia").

library(bnlearn)
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

