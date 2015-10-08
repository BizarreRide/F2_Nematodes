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


# Data processing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# create interaction factor for interaction between age_class and samcam, 
# to see if age_class effects stay constant
env.fin$agsam <- with(env.fin, interaction(age_class,samcam))

indices2 <- droplevels(cbind(indices,env.fin))
str(indices2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# A Closer Look on NCR
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
means = tapply(indices2$NCR,
               list(indices2$samcam, indices2$age_class),mean)
means

x11()
x.axis = unique(indices2$samcam)
par(fin=c(6.0,6.0),pch=18,mkh=.1,mex=1.5,
    cex=1.2,lwd=3)
matplot(c(1,4), c(0,1), type="n", 
        xlab="Sampling Campaign (seasons)", ylab="NCR",
        main= "Observed NCR Means")   
matlines(x.axis,means,type='l',lty=c(1,2,3,4,5))
matpoints(x.axis,means)    
legend(0.5,1,legend=c("A_Cm","SP_Y","SP_I1","SP_I2","SP_O"),lty=c(1,2,3),col=1:5,bty='n')



# Model formulas to be used
# model formula: y ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1 + Error(field.ID)
# model formula2: y ~ age + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1 + Error(field.ID)


#NCR is a precentage and should be analysed with binomial glm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test for Heteroscedasticity:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b <- c(1:18)
for (i in 1:18) {
 a <-  fligner.test(indices[,i] ~ env.fin$agsam)
 b[i] <- a$p.value
}
any(b<0.05)
# There is no violation of heterosecedasticity for any of the indices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Repeated measurments ANOVA 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# create Matrix to extract F-Values
f.stats <- matrix(NA,11,19)
# create Matrix to extract p-Values
p.stats <- matrix(NA,11,19)

# rmANOVAs
for (i in 1:18) {
  x <-  aov(indices[,i] ~ age_class + Error(field.ID/samcam), env.fin) #*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1 
  print(names(indices[i]))
  print(summary(x))
  nam <- paste("rmAOV",i,names(indices[i]), sep = ".")
  assign(nam, x)
  f.stats[,i+1] <- round(summary(x)$'Error: field.ID'[[1]]$'F value',2)
  p.stats[,i+1] <- round(summary(x)$'Error: field.ID'[[1]]$'Pr(>F)',3)
}

f.stats[,1] <- p.stats[,1]  <- row.names(summary(x)$'Error: field.ID'[[1]])
colnames(f.stats)[2:19] <- colnames(p.stats)[2:19] <- colnames(indices)
colnames(f.stats)[1] <- colnames(p.stats)[1] <- "Env"



# TukeyHSD test

for (i in 1:18) {
  df <- get(paste("rmAOV",i,names(indices[i]), sep = "."))
  x <- TukeyHSD(df)
  print(cld(x))
  }

# No Tukey Tests for aov objects!!!

y <- indices2[,1]
ev <- env.fin
y <- indices[!env.fin$crop=="Maize",1]
ev <- env.fin[!env.fin$crop=="Maize",]
x <- lme(y ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random= ~1 | field.ID, correlation=corCompSymm(form=~1|field.ID),method="ML", ev)
anova(x)
x <- lme(y ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random=list(field.ID=pdCompSymm(~samcam-1)), method="REML", ev)
anova(x)


indicesx 
indices2 <- indices2[!env.fin$crop=="Maize",]
indices2 <- indicesx

#**********************************************************************************************************************************************

# Repeat with glm from nlme package####
f.stats.lme <- matrix(NA,11,20)
p.stats.lme <- matrix(NA,11,20)

for(i in 1:18) {
  indices2$y <- indices2[,i]
  model <- lme(y ~ age_class*samcam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1, random=list(field.ID=pdCompSymm(~samcam-1)), method="REML", indices2)
  print(summary(model))
  nam <- paste("lme",i,names(indices[i]), sep = ".")
  assign(nam, model)
  f.stats.lme[,i+2] <- round(Anova(model, type = "II")$"Chisq",2)
  p.stats.lme[,i+2] <- round(Anova(model, type = "II")$"Pr(>Chisq)",3)
}
f.stats.lme[,1] <- p.stats.lme[,1]  <- row.names(Anova(model, type="II"))
f.stats.lme[,2] <- p.stats.lme[,2]  <- Anova(model, type = "II")$"Df"
colnames(f.stats.lme)[3:20] <- colnames(p.stats.lme)[3:20] <- colnames(indices)
colnames(f.stats.lme)[1] <- colnames(p.stats.lme)[1] <- "Env"
colnames(f.stats.lme)[2] <- colnames(p.stats.lme)[2] <- "Df"


#write.csv(p.stats.lme, file="p-Wertelme.csv")
#write.csv(f.stats.lme, file="Chi2-Wertelme.csv")


# Repeat with glm from lme4 package####
# Anova Type II
require(lme4)
f.stats.lmer2 <- matrix(NA,11,20)
p.stats.lmer2 <- matrix(NA,11,20)

for(i in 1:18) {
  indices2$y <- indices2[,i]
  model <- lmer(y ~ age_class*samcam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1 + (1|field.ID), indices2)
  print(summary(model))
  nam <- paste("lme",i,names(indices[i]), sep = ".")
  assign(nam, model)
  f.stats.lmer2[,i+2] <- round(Anova(model, type = "II")$"Chisq",2)
  p.stats.lmer2[,i+2] <- round(Anova(model, type = "II")$"Pr(>Chisq)",3)
}
f.stats.lmer2[,1] <- p.stats.lmer2[,1]  <- row.names(Anova(model, type="II"))
f.stats.lmer2[,2] <- p.stats.lmer2[,2]  <- Anova(model, type = "II")$"Df"
colnames(f.stats.lmer2)[3:20] <- colnames(p.stats.lmer2)[3:20] <- colnames(indices)
colnames(f.stats.lmer2)[1] <- colnames(p.stats.lmer2)[1] <- "Env"
colnames(f.stats.lmer2)[2] <- colnames(p.stats.lmer2)[2] <- "Df"


#write.csv(p.stats.lmer2, file="p-WertelmerII.csv")
#write.csv(f.stats.lmer2, file="Chi2-WertelmerII.csv")

# Anova Type III
f.stats.lmer3 <- matrix(NA,9,20)
p.stats.lmer3 <- matrix(NA,9,20)

for(i in 1:18) {
  indices2$y <- indices2[,i]
  #model <- gls(y ~ age_class*samcam + pH + c + mc + ata1, indices2) # reduced model, since fety delivers better reults with glm() and no random term
  model <- lmer(y ~ age_class*samcam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1 + (1|field.ID), indices2)
  print(summary(model))
  nam <- paste("lme",i,names(indices[i]), sep = ".")
  assign(nam, model)
  f.stats.lmer3[,i+2] <- round(Anova(model, type = "III")$"Chisq",2)
  p.stats.lmer3[,i+2] <- round(Anova(model, type = "III")$"Pr(>Chisq)",3)
}
f.stats.lmer3[,1] <- p.stats.lmer3[,1]  <- row.names(Anova(model, type="III"))
f.stats.lmer3[,2] <- p.stats.lmer3[,2]  <- Anova(model, type = "III")$"Df"
colnames(f.stats.lmer3)[3:20] <- colnames(p.stats.lmer3)[3:20] <- colnames(indices)
colnames(f.stats.lmer3)[1] <- colnames(p.stats.lmer3)[1] <- "Env"
colnames(f.stats.lmer3)[2] <- colnames(p.stats.lmer3)[2] <- "Df"

#write.csv(p.stats.lmer3, file="p-WertelmerIII.csv")
#write.csv(f.stats.lmer3, file="Chi2-WertelmerIII.csv")

# Which Type of Anova?




# Post Hoc Tukey Test with false discovery rate adjustment
Letters <- matrix(NA,10,18)
for(i in 1:18) {
  indices2$x <- indices2[,i]
  #model <- lm(x ~ agsam + pH + c + mc + ata1, data=indices2)
  #model <- lme(x ~ agsam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID, data=indices2)
  model <- lmer(x ~ agsam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1 + (1|field.ID), indices2) # No Difference
  x <- summary(glht(model, linfct=mcp(agsam="Tukey")), test = adjusted(type = "bonferroni"))
  nam <- paste("glht",i,names(indices[i]), sep = ".")
  assign(nam, x)
  x <- cld(x)$mcletters
  xx <- x$Letters
  Letters[,i] <- xx
  row.names(Letters) <- levels(env.fin$agsam)
}
colnames(Letters) <- colnames(indices)

par(mar=c(5,5,7,2))
plot(cld(glht.4.PPI))

# NCR ####

fety2$bacfu <- fety2$bacterivore + fety2$fungivore
fety2$age_class <- as.ordered(fety2$age_class)

model <- glmer(fungivore/bacfu ~ age_class*samcam + pH + c + mc + ata1 + (1|field.ID),weights=bacfu,family="binomial",fety2, control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
#model <- glmer(y ~ age_class*samcam + pH + c + mc + ata1 + (1|field.ID),weights=N,family="binomial",fety2,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
print(summary(model))
Anova(model, type = "III")

model <- glmer(fungivore/bacfu ~ agsam + pH + c + mc + ata1 + (1|field.ID),weights=bacfu,family="binomial",fety2, control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
x <- summary(glht(model, linfct=mcp(agsam="Tukey")), test = adjusted(type = "bonferroni"))
par(mar=c(5,5,8,2))
plot(cld(x))



#write.csv(Letters, file="LettersFDR.csv")
#write.csv(Letters, file="LettersBonferroni.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




##### Plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

se <- function(x) {
  return(c(ymin=(mean(x) - (sd(x)/sqrt(length(x)))), ymax=(mean(x) +
                                                             (sd(x)/sqrt(length(x))))))
}

require(reshape2)
indices.melt <- melt(indices2, id.vars=19:35)
g1 <- ggplot(indices.melt, aes(x = age_class, y = value)) +
  stat_summary(fun.data = 'se', geom = 'errorbar', width = 0.2, size = 1) + 
  #stat_summary(fun.y = mean, geom = 'point', size = 3, color = 'red') +
  stat_summary(fun.y = mean, geom = 'bar', size = 1, color = 'gray20', fill="lightpink4", alpha=1/3) +
  facet_grid( ~ samcam) +
  #ylab(levels(variable)) +
  #facet_wrap(~ variable, scales="free") +
  theme_bw()

pl1 = plyr::dlply(indices.melt, "variable", `%+%`, e1 = g1)
pl1$PPI

do.call(gridExtra::grid.arrange, pl)
do.call(gridExtra::grid.arrange, pl[c(4,7,15,16,17,18)])
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




# Feeding Types ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env.sc <- env.fin[,sapply(env.fin,is.numeric)]
env.sc <- data.frame(scale(env.sc, center=TRUE, scale=TRUE))
env.sc <- cbind(env.fin[,!sapply(env.fin,is.numeric)],env.sc)
str(env.sc)


fety.av <- fety/indices$N
fety2 <- droplevels(cbind(fety.av,env.sc,N=indices$N))

fety2 <- droplevels(cbind(fety,env.sc, N=indices$N))
str(fety2)

f.stats.glmer3 <- matrix(NA,7,7)
p.stats.glmer3 <- matrix(NA,7,7)

for(i in 1:5) {
  fety2$y <- fety2[,i]
  model <- glm(y ~ age_class*samcam + pH + c + mc + ata1,weights=N,family="binomial",fety2)
  #model <- glmer(y ~ age_class*samcam + pH + c + mc + ata1 + (1|field.ID),weights=N,family="binomial",fety2,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
  print(summary(model))
  nam <- paste("lme",i,names(fety2[i]), sep = ".")
  assign(nam, model)
  f.stats.glmer3[,i+2] <- round(Anova(model, type = "III")$"LR Chisq",2)
  p.stats.glmer3[,i+2] <- round(Anova(model, type = "III")$"Pr(>Chisq)",3)
}
f.stats.glmer3[,1] <- p.stats.glmer3[,1]  <- row.names(Anova(model, type="III"))
f.stats.glmer3[,2] <- p.stats.glmer3[,2]  <- Anova(model, type = "III")$"Df"
colnames(f.stats.glmer3)[3:7] <- colnames(p.stats.glmer3)[3:7] <- colnames(fety)
colnames(f.stats.glmer3)[1] <- colnames(p.stats.glmer3)[1] <- "Env"
colnames(f.stats.glmer3)[2] <- colnames(p.stats.glmer3)[2] <- "Df"

#write.csv(p.stats.glmer3, file="p-WerteglmerIII.csv")
#write.csv(f.stats.glmer3, file="Chi2-WerteglmerIII.csv")

# Post Hoc Tukey Test with false discovery rate adjustment
LettersFety <- matrix(NA,10,5)
for(i in 1:5) {
  fety2$x <- fety2[,i]
  #model <- glmer(x ~ agsam + (1|field.ID),weights=N,family="binomial",fety2,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000))) 
  model <- glm(x ~ agsam + pH + c + mc + ata1,weights=N,family="binomial",fety2) 
  x <- summary(glht(model, linfct=mcp(agsam="Tukey")), test = adjusted(type = "none"))
  nam <- paste("glht",i,names(fety[i]), sep = ".")
  assign(nam, x)
  x <- cld(x)$mcletters
  xx <- x$Letters
  LettersFety[,i] <- xx
  row.names(LettersFety) <- levels(env.fin$agsam)
}
colnames(LettersFety) <- colnames(fety)

glht.3.fungivore
#write.csv(LettersFety, file="LettersFety.csv")

par(mar=c(5,4,10,2))
plot(cld(glht.3.fungivore))
plot(cld(glht.4.herbivore))


require(reshape2)
fety.melt <- melt(fety2, id.vars=6:21)
g2 <- ggplot(fety.melt, aes(x = age_class, y = value)) +
  stat_summary(fun.data = 'se', geom = 'errorbar', width = 0.2, size = 1) + 
  #stat_summary(fun.y = mean, geom = 'point', size = 3, color = 'red') +
  stat_summary(fun.y = mean, geom = 'bar', size = 1, color = 'gray20', fill="lightpink4", alpha=1/3) +
  facet_grid( ~ samcam) +
  #ylab() +
  #facet_wrap(~ variable, scales="free") +
  theme_bw()
g2
pl2 = plyr::dlply(fety.melt, "variable", `%+%`, e1 = g2)

do.call(gridExtra::grid.arrange, pl2)

with(fety, hist(log1p(fety$carnivore)))
with(fety, hist(fety$bacterivore))
with(fety, hist(fety$herbivore))
with(fety, hist(fety$omnivore))
with(fety, hist(fety$fungivore))




q <- aov(fungivore ~ age_class*samcam, fety2)
summary(q)

# NCR ####

indices2$bacfu <- fety2$bacterivore + fety2$fungivore

model <- glmer(fungivore/bacfu ~ age_class*samcam + pH + c + mc + ata1 + (1|field.ID),weights=bacfu,family="binomial",fety2, control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
#model <- glmer(y ~ age_class*samcam + pH + c + mc + ata1 + (1|field.ID),weights=N,family="binomial",fety2,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
print(summary(model))
Anova(model, type = "III")

model <- glmer(fungivore/bacfu ~ agsam + pH + c + mc + ata1 + (1|field.ID),weights=bacfu,family="binomial",fety2, control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
x <- summary(glht(model, linfct=mcp(agsam="Tukey")), test = adjusted(type = "none"))
par(mar=c(5,5,8,2))
plot(cld(x))





