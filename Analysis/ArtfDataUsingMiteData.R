
# Load Data
data(mite)
data(mite.env)

spe <- mite[1:30,1:25]
env1 <- mite.env[1:30,]

#add my original factors
env <- data.frame(matrix(NA,30,4))
colnames(env) <- c("field.ID", "age_class", "SamCam","crop")
env[,1] <- factor(c(seq(1,12,1),13,14,15,seq(1,12,1),16,17,18))
env[,2] <- factor(rep(rep(letters[1:5],each=3),2))
env[,3] <- factor(rep(c(1,2),each=15))
env[,4] <- factor(rep(c(rep("A",12),rep("B",3)),2))

env <- cbind(env,env1)

#add age_class effect
spe[env$age_class=="a",1:5] <-  spe[env$age_class=="a",1:5] + sample(c(14:31),7, replace=T)
spe[env$age_class=="d",10:15] <-  spe[env$age_class=="d",10:15] - sample(c(14:31),7, replace=T)
spe[spe<0] <- 0

# perform partial db_RDA

#A with condition(field.ID)
spe.db.repmes0 <- capscale(spe ~ 1 + Condition(field.ID), distance="bray", data=env, add=T) 
spe.db.repmes1 <- capscale(spe ~ . + Condition(field.ID), distance="bray", data=env, add=T) 
#spe.db.repmes1 <- capscale(spe ~ age_class + Condition(field.ID), distance="bray", data=env, add=T) 
stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
anova(stepping)

# no factor age_class in the plots:
x11()
par(mfrow=c(1,2))
plot(spe.db.repmes1, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(spe.db.repmes1, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")


# without Condition, while removing field.ID from the env dataset
spe.db.repmes0 <- capscale(spe ~ 1 , distance="bray", data=env[,-1], add=T) 
spe.db.repmes1 <- capscale(spe ~ . , distance="bray", data=env[,-1], add=T) 
stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
anova(stepping)



# age_class has a significant influence on community composition!!


spe <- droplevels(spe[!env$age_class=="e",])
env <- droplevels(env[!env$age_class=="e",])

spe.db.repmes0 <- capscale(spe ~ 1 + Condition(field.ID), distance="bray", data=env, add=T) 
spe.db.repmes1 <- capscale(spe ~ . + Condition(field.ID), distance="bray", data=env, add=T) 
spe.db.repmes1 <- capscale(spe ~ age_class + Condition(field.ID), distance="bray", data=env, add=T) 
stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
anova(stepping)

x11()
par(mfrow=c(1,2))
plot(spe.db.repmes1, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(spe.db.repmes1, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")





data(mite)

data(mite.env)

str(mite.env)

mite.hel = decostand(mite, "hel")

mod0 <- rda(mite.hel ~ 1 + Condition(Topo), mite.env)  # Model with intercept only

mod1 <- rda(mite.hel ~ . + Condition(Topo), mite.env)  # Model with all explanatory variables

step.res <- ordiR2step(mod0, mod1, perm.max = 200)
step.res
step.res$anova
