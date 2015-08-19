# Plot #####
x11()
par(mfrow=c(1,2))
plot(spe.db.repmes1, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(spe.db.repmes1, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")

#####

Neugier <- matrix(0,7,20)

for(k in 1:20) {
    
  source("Analysis/ArtificialData.R")
  
  
  spe1 <- spe
  env1 <- env
  num1 <- paste("spe",k,sep=".")
  assign(num1,spe1)
  num2 <- paste("env",k,sep=".")
  assign(num2,env1)
  
  #Superduper1 <- spe.15
  #Superduper2 <- env.15
  
  #write.csv(Superduper1,"superduper1.csv")
  #write.csv(Superduper2,"superduper2.csv")
  
  #A
  spe.db.repmes0 <- capscale(spe1 ~ 1 + Condition(field.ID), distance="bray", data=env1, add=T) 
  spe.db.repmes1 <- capscale(spe1 ~ . + Condition(field.ID), distance="bray", data=env1, add=T) 
  stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
  Neugier[1,k] <-  anova(stepping)$"F"[1]
  
  
  #B
  spe.db.repmes0 <- capscale(spe1 ~ 1 , distance="bray", data=env1[,-1], add=T) 
  spe.db.repmes1 <- capscale(spe1 ~ . , distance="bray", data=env1[,-1], add=T) 
  stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1, direction="forward",pstep=1000,trace=F)
  Neugier[2,k] <- anova(stepping)$"F"[1]
  
  #F
  spe.db.repmes0 <- capscale(spe1 ~ 1 + Condition(clay), distance="bray", data=env1[,-1], add=T) 
  spe.db.repmes1 <- capscale(spe1 ~ . + Condition(clay), distance="bray", data=env1[,-1], add=T) 
  stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
  Neugier[3,k] <- anova(stepping)$"F"[1]
  
  #C
  spe.db.repmes0 <- capscale(spe1 ~ 1 , distance="bray", data=env1[,-1], add=T) 
  spe.db.repmes1 <- capscale(spe1 ~ age_class , distance="bray", data=env1[,-1], add=T) 
  stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1, direction="forward",pstep=1000,trace=F)
  Neugier[4,k] <- anova(stepping)$"F"[1]
  
  spe1 <- spe[!env$crop=="B",]
  spe1 <- droplevels(spe1)
  env1 <- env[!env$crop=="B",]
  env1 <- droplevels(env1)
  
  #D
  spe.db.repmes0 <- capscale(spe1 ~ 1 + Condition(field.ID), distance="bray", data=env1, add=T)
  spe.db.repmes1 <- capscale(spe1 ~ . + Condition(field.ID), distance="bray", data=env1, add=T)
  stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
  Neugier[5,k] <- anova(stepping)$"F"[1]
  
  #E
  spe.db.repmes0 <- capscale(spe1 ~ 1 , distance="bray", data=env1[,-1], add=T) 
  spe.db.repmes1 <- capscale(spe1 ~ . , distance="bray", data=env1[,-1], add=T) 
  stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
  Neugier[6,k] <- anova(stepping)$"F"[1]
  
  #G
  spe.db.repmes0 <- capscale(spe1 ~ 1 + Condition(clay), distance="bray", data=env1[,-1], add=T) 
  spe.db.repmes1 <- capscale(spe1 ~ . + Condition(clay), distance="bray", data=env1[,-1], add=T) 
  stepping <- ordiR2step(spe.db.repmes0 , spe.db.repmes1,direction="forward",pstep=1000,trace=F)
  Neugier[7,k] <- anova(stepping)$"F"[1]

}


x11()
par(mfrow=c(1,2))
plot(spe.db.repmes1, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(spe.db.repmes1, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")






spe.db.repmes <- capscale(spe ~ age_class + n + temp + Condition(field.ID), distance="bray", data=env, add=TRUE)
summary(spe.db.repmes, display=NULL)

anova(spe.db.repmes, step=1000, perm.max=1000)
anova(spe.db.repmes, by="axis", step=1000, perm.max=1000)
# 2 axes are significant
# However field.ID was not used as a factor!!



#Triplots:
x11()
par(mfrow=c(1,2))
plot(spe.db.repmes1, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(spe.db.repmes1, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")

x11()
par(mfrow=c(1,2))
plot(spe.db.repmes, scaling=1, main="Triplot RDA - scaling 1 - wa scores")
plot(spe.db.repmes, scaling=2, main="Triplot RDA - scaling 2 - wa scores")