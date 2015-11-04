#################
# F2_Nematodes
# GLMM Power Analysis
# Quentin Schorpp
# 04.11.2015
#################


fakeData <- data.frame(ID = c(1:12,13,14,15,1:12,16,17,18),
                       time = rep(c("A","B"), each=15),
                       age_class = rep(rep(c("C","D","E","F","G"),each=3),2),
                       response = rpois(30,4)/100,
                       N = rep(100,30)
                       )

fakeModel <- glmer(response ~ time*age_class + (1|ID), family = binomial, weights=N, data=fakeData)

datfun <- function(nindiv=15, nperindiv=2, nage_class=3, nsim=3) {
                    data.frame(ID=factor(rep(1:nindiv, nperindiv)),
                               time = factor(rep(c("A", "B"), each=nindiv)),
                               age_class = rep(rep(c("C","D","E","F","G"),each=nage_class),nperindiv),
                               response = rpois(nindiv*nperindiv, 4)/100,
                               N = rep(100, nindiv*nperindiv))
}

simdat <- simulate(formula(fakeModel), newdat=datfun(),
                   newparam=list(theta=getME(fakeModel,"theta"),
                                 beta=getME(fakeModel,"beta")),
                   family=binomial,
                   weights=N)

getpower <- function(nsim,params=list(),alpha=0.05) {
  out1 <- replicate( nsim, do.call("sim1",as.list(params)))
  mean(out1<alpha)
}




