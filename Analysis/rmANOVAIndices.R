###########################
# F2 Nematodes
# Univaraite multiple repeated measurements Tests
# Quentin Schorpp
# 11.08.2015
###########################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env.fin$agsam <- with(env.fin, interaction(age_class,samcam))


# model formula: y ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1 + Error(field.ID)

for (i in 1:18) {
 x <-  fligner.test(indices[,i] ~ with(env.fin, interaction(age_class,samcam)))
 print(x)
}


f.stats <- matrix(NA,12,19)
p.stats <- matrix(NA,12,19)
for (i in 1:18) {
  x <-  aov(indices[,i] ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1 + Error(field.ID), env.fin)
  print(names(indices[i]))
  print(summary(x))
  nam <- paste("rmAOV",i,names(indices[i]), sep = ".")
  assign(nam, x)
  f.stats[,i+1] <- round(summary(x)$'Error: Within'[[1]]$'F value',2)
  p.stats[,i+1] <- round(summary(x)$'Error: Within'[[1]]$'Pr(>F)',3)
}

f.stats[,1] <- p.stats[,1]  <- row.names(summary(x)$'Error: Within'[[1]])
colnames(f.stats)[2:19] <- colnames(p.stats)[2:19] <- colnames(indices)
colnames(f.stats)[1] <- colnames(p.stats)[1] <- "Env"


aggregate(indices,list(env.fin$age_class),mean)

a <- aggregate(indices,list(env.fin$age_class,env.fin$field.ID),mean)
a <- aggregate(a, list(a$Group.1),mean)
a

# TukeyHSD test

for(i in 1:6) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
  nam <- paste("o", i, sep = ".")
  assign(nam, x)
}

for (i in 1:18){
  df <- get(paste(i,"rmAOV",names(indices[i]), sep = "."))
  x <- TukeyHSD(df)
  print(x)
}


indices2 <- cbind(indices,env.fin)

f.stats.lme <- matrix(NA,12,19)
p.stats.lme <- matrix(NA,12,19)

for(i in 1:18) {
  indices2$x <- indices2[,i]
  model <- lme(x ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID,data=indices2)
  print(summary(model))
  nam <- paste("lme",i,names(indices[i]), sep = ".")
  assign(nam, model)
  f.stats.lme[,i+1] <- round(anova(model)$"F-value",2)
  p.stats.lme[,i+1] <- round(anova(model)$"p-value",3)
}

anova(lme.1.N)$"F-value"
anova(lme.1.N)$"p-value"

Letters <- matrix(NA,5,18)
for(i in 1:18) {
  indices2$x <- indices2[,i]
  model <- lme(x ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID, data=indices2)
  x <- summary(glht(model, linfct=mcp(age_class="Tukey")), test = adjusted(type = "fdr"))
  nam <- paste("glht",i,names(indices[i]), sep = ".")
  assign(nam, x)
  x <- cld(x)$mcletters
  xx <- x$Letters
  Letters[,i] <- xx
  row.names(Letters) <- levels(env.fin$age_class)
}

summary(glht.1.N)

