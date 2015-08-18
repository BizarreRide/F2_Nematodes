

# Create an artificial dataset as reproducible example
spe <- matrix(NA,30,15)
colnames(spe) <- paste("Sp",c(1:15), sep=".")
I <- c(29,15,2,2,1,22,18,75,36,2,1,11,12,2,5,5,7,3,2,3,52,1,63,76,7) # max. observed
for(i in 1:15){
  J <- sample(I,1,replace = F)
  spe[,i] <- sample(c(0:J), 30, replace = TRUE, prob = c(0.013,rep(sample(seq(0.008,0.013333,length.out=20),1,replace=F),J)))
}


env.means <- structure(c(6.337, 19.39, 0.17, 25.54, 11.52, 1.70, -3.72e-18), 
                       .Names = c("pH","mc", "n", "clay", "temp", "prec", "intensity"))

env.sd <- structure(c(0.73, 5.86, 0.07, 7.16, 3.18, 0.63, 0.65),
                    .Names = c("pH", "mc", "n", "clay", "temp", "prec", "intensity"))


env <- data.frame(matrix(NA,30,10))
colnames(env) <- c("field.ID", "age_class", "SamCam","pH", "mc", "n", "clay", "temp", "prec", "intensity")
env[,1] <- factor(c(seq(1,12,1),13,14,15,seq(1,12,1),16,17,18))
env[,2] <- factor(rep(rep(letters[1:5],each=3),2))
env[,3] <- factor(rep(c(1,2),each=15))

for(i in 4:10){
    env[,i] <- round(rnorm(30, mean=env.means[i-3], sd=env.sd[i-3]),2)
  }

#env <- data.frame(env, stringsAsFactors=T)

### add treatment effects,
spe[env$age_class %in% c("a", "b","c","d"),] <-  spe[env$age_class %in% c("a", "b","c","d"),] + sample(c(3:5),5, replace=T)
spe[env$age_class=="a",1:5] <-  spe[env$age_class=="a",] + sample(c(9:12),5, replace=T)
spe[env$age_class=="a",5:10] <-  spe[env$age_class=="a",] + sample(c(4:6),5, replace=T)
spe[env$age_class=="a",c(3,5,7)] <-  spe[env$age_class=="a",] + sample(c(6:9),5, replace=T)
spe[env$age_class=="e",4:12] <-  spe[env$age_class=="e",] - sample(c(3:5),5, replace=T)
spe[env$age_class=="e",10:15] <-  spe[env$age_class=="e",] + sample(c(3:5),5, replace=T)

spe[spe<0] <- 0
spe <- data.frame(spe, stringsAsFactors=F)
str(spe)
str(env)

env$n <- env$n + (1/spe[,5])*env$n
env$temp <- env$temp + 0.07*spe[,5]*env$temp


#env$field.ID <- as.numeric(env$field.ID)



spe.db.repmes <- capscale(spe ~ . + Condition(field.ID), distance="bray", data=env, add=T)
summary(spe.db.repmes)
stepping <- ordiR2step(capscale(spe ~ 1 + Condition(field.ID),distance="bray", data=env, add=T), scope=formula(spe.db.repmes),direction="forward",pstep=1000,trace=F)
anova(stepping)

spe.db.repmes <- capscale(spe ~ age_class + n + temp + Condition(field.ID), distance="bray", data=env, add=TRUE)
summary(spe.db.repmes, display=NULL)

anova(spe.db.repmes, step=1000, perm.max=1000)
anova(spe.db.repmes, by="axis", step=1000, perm.max=1000)
# 2 axes are significant
# However field.ID was not used as a factor!!



#Triplots:
x11()
par(mfrow=c(1,2))
plot(spe.db.repmes, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(spe.db.repmes, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")

x11()
par(mfrow=c(1,2))
plot(spe.db.repmes, scaling=1, main="Triplot RDA - scaling 1 - wa scores")
plot(spe.db.repmes, scaling=2, main="Triplot RDA - scaling 2 - wa scores")

