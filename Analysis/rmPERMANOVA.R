###########################
# F2 Nematodes
# PCA and RDA Analysis
# Quentin Schorpp
# 06.08.2015
###########################


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fam.rp <- fam[!env1$crop=="Maize",]

time <- as.ordered(rep(1:2,each=12))
site <- gl(12,1,24)

cbind(site, time, fam.rp)


### computing the true R2-value
### (btw, using dist() defaults to euclidean distance):
print(fit <- adonis(fam.rp ~ time, method="bray", permutations=1)) # "NO MORE dist() is needed"
str(fit$aov.tab)

### number of perms
B <- 4095
# using the whole set of enumartions (4096), what would be referred to as "exact test", would be best. 
# But, with as much as 1999 permutation you yield a good approximation of the null-distribution.

### setting up frame which will be populated by
### random r2 values:
pop <- rep(NA, B + 1)

### the first entry will be the true r2:
pop[1] <- fit$aov.tab[1, 5]


### set up a "permControl" object:
### we turn off mirroring as time should only flow in one direction
ctrl <- how(blocks = site, within = Within(type = "series", mirror = FALSE)) #"permControl" is not used at all

### Number of observations:
nobs <- nrow(fam.rp)

### check permutation (...rows represent the sample id):
### ..they are ok!
### within in each repeated sample (= sites) timepoints are shuffled,
### with keeping the sequence intact (e.g., for site 1: 1,2,3 - 2,3,1 - 3,2,1)
shuffle(nobs, control = ctrl)

### loop:
### in adonis(...) you need to put permutations = 1, otherwise
### adonis will not run
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(fam.rp ~ time[idx], method="bray", permutations = 1) # "NO MORE dist() is needed"
  pop[i] <- fit.rand$aov.tab[1, 5]
}


### get the p-value:
print(pval <- sum(pop >= pop[1]) / (B + 1))
### [1] 0.0445 # <- 4.5 % of the p-values are >= pop[1]

### the sign. p-value supports the H1 (->there is a time effect).
### ..and the fact that samples are not iid is allowed by
### the customized perms - so this p-value is trustworthy as opposed
### to tests not acknowledging dependency of data points..

## make a histogram to see random R2-values and the true one:
hist(pop, xlab = "Population R2")
abline(v = pop[1], col = 2, lty = 3)
text(0.052, 300, paste("true R2,\np = ", pval, sep = ""))

adonis(fam.rp ~ time, method="bray", strata=site, permutations=4095)


# The F2_Nematode Analysis has 2 time points and 12 sites, hence a minimal p-value of 1/4096 = 0.00024 can be obtained

# We say it doesn't matter that all samples at t1 were taken at one time and the ones of t2 at another. 
# This means we ignore the correlation structure of all samples taken at t1
# Still there would be the repeated measures and the fact that the samples were taken in a temporal sequence. 

# we will not permute across these samples taken at one site, hence the ordering will only be 1,2 or 2,1

# toroidal shift which is keeping the "neighborhood" of values intact

# in this case, where only two timepoints exist, this needs probably another idea (i.e. two sample T-Test?)



# Eample from biobucket

### species matrix with 20 species abundances (mean = 50, sd = 10)
### one time variable, with 3 timepoints, which should be tested
### and a factor denoting sites that were repeatedly sampled (site)

## Load packages
require(vegan)

### Data:
sp <- matrix(rnorm(3 * 6 * 20, 50, 10), nrow = 3 * 6, ncol = 20,
             dimnames = list(1:18, paste("Sp", 1:20, sep = "")))

time <- as.ordered(rep(1:3, 6))
site <- gl(6, 3)
cbind(site, time, sp)

### add time effect at timepoint 3,
### this will effect will be tested by adonis():
sp_1 <- sp
sp_1[time==3,] <- sp[time==3,] + rnorm(20, 10, 1)
cbind(site, time, sp_1)

### choose which species set to test:
test_sp <- sp_1

### computing the true R2-value

### (btw, using dist() defaults to euclidean distance):
print(fit <- adonis(test_sp ~ time, permutations=1)) # "NO MORE dist() is needed"

### number of perms
B <- 1999

### setting up frame which will be populated by
### random r2 values:
pop <- rep(NA, B + 1)

### the first entry will be the true r2:
pop[1] <- fit$aov.tab[1, 5]

### set up a "permControl" object:
### we turn off mirroring as time should only flow in one direction
ctrl <- how(blocks = site, within = Within(type = "series", mirror = FALSE)) #"permControl" is not used at all

### Number of observations:
nobs <- nrow(test_sp)

### check permutation (...rows represent the sample id):
### ..they are ok!
### within in each repeated sample (= sites) timepoints are shuffled,
### with keeping the sequence intact (e.g., for site 1: 1,2,3 - 2,3,1 - 3,2,1)
shuffle(nobs, control = ctrl)

### loop:
### in adonis(...) you need to put permutations = 1, otherwise
### adonis will not run
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(test_sp ~ time[idx],permutations = 1) # "NO MORE dist() is needed"
  pop[i] <- fit.rand$aov.tab[1, 5]
}

### get the p-value:
print(pval <- sum(pop >= pop[1]) / (B + 1))
### [1] 0.0105

### the sign. p-value supports the H1 (->there is a time effect).
### ..and the fact that samples are not iid is allowed by
### the customized perms - so this p-value is trustworthy as opposed
### to tests not acknowledging dependency of data points..

### test sp set without an effect:
### replace test_sp with sp set without effect:
test_sp <- sp

### now re-run the script and see the result:
### it is insign. - as expected:

### setting up frame which will be populated by
### random r2 values:
pop <- rep(NA, B + 1)

### computing the true R2-value:
print(fit <- adonis(test_sp ~ time, permutations = 1))# "NO MORE dist() is needed"

### the first entry will be the true r2:
pop[1] <- fit$aov.tab[1, 5]

### run the loop:
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(test_sp ~ time[idx], permutations = 1)# "NO MORE dist() is needed"
  pop[i] <- fit.rand$aov.tab[1, 5]
}
print(pval <- sum(pop >= pop[1]) / (B + 1))
### [1] 0.7605

## make a histogram to see random R2-values and the true one:
hist(pop, xlab = "Population R2")
abline(v = pop[1], col = 2, lty = 3)
text(0.08, 300, paste("true R2,\np = ", pval, sep = ""))



# Probably this solves my Problems:

# See if there is a n increae over time
# glmm(y ~ age + (1|field.ID), data)

# And then see if there are differences in proportions of families between age_classes

require(lme4)

nema1 <- cbind(fam.fin, env.fin)
nema2 <- nema1[!nema1$crop=="Maize",]

fit1 <- glmer(Aphelenchidae ~ age + (1|field.ID), family="poisson", nema2)
fit1 <- glmer(Mononchidae ~ age + (1|field.ID), family="poisson", nema2)
fit1 <- glmer(Plectidae ~ age + (1|field.ID), family="poisson", nema2)
fit1 <- glmer(Cephalobidae ~ age + (1|field.ID), family="poisson", nema2)

summary(fit1)
Anova(fit1)

nema.slope <- cbind(fam.slope, env.av[1:12,])

nema.slope$age_class <- as.ordered(nema.slope$age_class)

fit2 <- lm(Hoplolaimidae ~ age_class, nema.slope)
 
summary(fit2)



