###########################
# F2 Nematodes
# Univaraite multiple repeated measurements Tests
# Quentin Schorpp
# 11.08.2015
###########################


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 

asinTransform <- function(p) { asin(sqrt(p)) }
library(influence.ME)
library(lmerTest)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# as needed
# source("Data/DataProcessing/MaturityIndices.R") 
# source("Data/DataProcessing/FaunalProfileIndices.R") 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create interaction factor for interaction between age_class and samcam, 
# to see if age_class effects stay constant

env.fin$agsam <- with(env.fin, interaction(age_class,samcam))
#env.fin$c <-  env1$c # include carbon

#fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)

indices <- droplevels(cbind(env.fin, location=env1$location, N=rowSums(fam.org)))
indices$ID <- 1:nrow(indices)
indices$nsamcam <- as.numeric(factor(indices$samcam))

indices.backup <- indices
indices <- indices.backup
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis of FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


fety.backup <- fety/indices$N
fety2 <- fety.backup[,-c(1:4)]

p <- ncol(fety2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(2:10)
for (i in 2:10) {
 fk <-  fligner.test(fety[,i] ~ env.fin$agsam)
 pvector[i] <- fk$p.value
 boxplot(fety[,i]~env.fin$agsam)
}
any(pvector<0.05)
pvector

# There is no violation of heterosecedasticity for any of the indices
# Null Hypothesis: All variances are equal can't be abgelehnt at the 5% level
# The probability to judge real diffeences in varainces although the NULL is true is too high
# in more than 5 of a hundred cases we would judge wrong if we say there are differences
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data summary ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- asinTransform(fety2[,i])
  plot(indices$y)
  boxplot(indices$y)
  hist(indices$y, main="")
  plot(indices$y^2)
  title(names(fety2)[i],outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Detecting Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- asinTransform(fety2[,i])
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$samcam)
  car::Boxplot(indices$y ~ indices$agsam)
  title(names(fety2)[i],outer=TRUE)
}

# outliers:
# fungivores: 2
# omnivores: 1


# age_class:
# bacterivores decrease
# herbivores increase
# fungivores have an polynomial relationship

# samcam:
# carnivores increase in the second year, only in Silphie!
# no big differences for the others
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#indices <- indices.backup[-c(9,15,6,23,20),]
#fety2 <- fety2[-c(9,15,6,23,20),]

indices$age_class2 <- indices$age_class
indices$age_class <- as.ordered(indices$age_class)
indices$age_class <- indices$age_class2

indices$samcam2 <- indices$samcam
indices$samcam <- indices$nsamcam
indices$samcam <- indices$samcam2

# Für AIC und die Residuean macht es keinen Unterschied ob die Faktoren geordnet sind oder nicht, 
# oder ob samcam numerisch ist oder nicht


# Bi☺nomial GLMM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.fety.biglmer <- matrix(NA,7,2+p)
colnames(p.fety.biglmer) <- c("Env", "DF", colnames(fety2)[1:p])

f.fety.biglmer <- p.fety.biglmer

fety.biglmer <- list()

for(i in 1:p) {
  indices$y <- fety2[,i]
  model <- glmer(y ~ age_class*samcam + pH + mc + n + ata1 + (1|field.ID), family=binomial(link=logit), weights=N, indices, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("fety",i,names(fety2)[i], sep = ".")
  assign(name, model)
  fety.biglmer[[i]] <- assign(name, model)
  f.fety.biglmer[,i+2] <- round(car::Anova(model)$"Chisq",2)
  p.fety.biglmer[,i+2] <- round(car::Anova(model)$"Pr(>Chisq)",3)
}
f.fety.biglmer[,1] <- p.fety.biglmer[,1]  <- row.names(Anova(model))
f.fety.biglmer[,2] <- p.fety.biglmer[,2]  <- Anova(model)$"Df"

mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(fety2)[i], sep = "."))}
names(fety.biglmer)[1:p] <- mod.names

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



### binomial GLMM - Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:p){ 
  # print(list(summary(fety.biglmer[[k]]),Anova(fety.biglmer[[k]], type="III")))
  #corvif(fety.biglmer[[k]])
  
  E1 <- resid(fety.biglmer[[k]], type="pearson")
  E2 <- resid(fety.biglmer[[k]], type="response")
  F1 <- fitted(fety.biglmer[[k]], type="response")
  P1 <- predict(fety.biglmer[[k]], type="response")
  
  par(mfrow=c(2,2),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  )
  # Plot fitted vs. residuals
  scatter.smooth(F1, E1, cex.lab = 1.5, xlab="Fitted values", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot predicted vs. residuals
  scatter.smooth(P1, E1, cex.lab = 1.5, xlab="Predicted values", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot fitted vs. predicted
  #scatter.smooth(F1, P1, cex.lab = 1.5, xlab="Fitted values", ylab="Predicted")
  abline(h = 0, v=0, lty=2)
  
  # Histogram of Residuals
  hist(E1, prob=TRUE, main = "", breaks = 20, cex.lab = 1.5, xlab = "Response Residuals", col="PapayaWhip")
  lines(density(E1), col="light blue", lwd=3)
  lines(density(E1, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
  
  title(names(fety.biglmer)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  # plot age_class vs. residuals
  boxplot(E1 ~ indices$age_class, cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$samcam,cex.lab = 1.5, xlab="samcam", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$agsam,cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  title(names(fety.biglmer)[k], outer=TRUE)
  
  indices$y <- fety2[,k]
  par(mfrow=c(1,1))
  scatter.smooth(F1,indices$y, cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
  #text(F1,indices$y,label=interaction(indices$field.ID,indices$samcam),col='red')
  abline(h = 0, v=0, lty=2)
  
  title(names(fety.biglmer)[k], outer=TRUE)
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Residuals against variables not in the model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(k in 1:p){ 
  #corvif(fety.biglmer[[k]])
  
  E1 <- resid(fety.biglmer[[k]], type="pearson")
  E2 <- resid(fety.biglmer[[k]], type="response")
  F1 <- fitted(fety.biglmer[[k]], type="response")
  P1 <- predict(fety.biglmer[[k]], type="response")
  
  par(mfrow=c(2,2),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  )
  
  # Plot fitted vs. residuals
  scatter.smooth(indices$n, E1, cex.lab = 1.5, xlab="Nitrogen", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot predicted vs. residuals
  scatter.smooth(indices$pH, E1, cex.lab = 1.5, xlab="pH", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot fitted vs. predicted
  scatter.smooth(indices$mc, E1, cex.lab = 1.5, xlab="mc", ylab="Predicted")
  abline(h = 0, v=0, lty=2)
  
  # plot fitted vs. predicted
  scatter.smooth(indices$ata1, E1, cex.lab = 1.5, xlab="Temperature", ylab="Predicted")
  abline(h = 0, v=0, lty=2)
  
  title(names(fety.biglmer)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Influence measures ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infl <- list()
for (k in 1:p) {
  infl[[k]] <- influence(fety.biglmer[[k]], obs=TRUE)
}

for(k in 1:p) {
  p <- plot(infl[[k]], which = "cook", main=mod.names[k])
  print(p)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Multimodel Averaging ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# First calculate new global model

library(MuMIn)
fety.dredge <- list()
options(na.action = "na.fail")
for(i in 1:p) {
  indices$y <- fety2[,i]
  model <- glmer(y ~ age_class*samcam + pH + mc + n + ata1 + (1|field.ID), family=binomial(link="logit"), weights=N, REML=FALSE, indices)
  M.dredge <- dredge(model, fixed=c(c("age_class", "samcam")))
  name <- paste("fetydredge.fx.on",i,names(fety2[i]), sep = ".")
  assign(name, M.dredge)
  fety.dredge[[i]] <- assign(name, M.dredge)
  print(head(M.dredge,10))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#####

fety.dredge.fixed.fafac <- fety.dredge

head(fety.dredge.facfac[[1]])
head(fety.dredge.ordnum[[1]])
head(fety.dredge.fixed.fafac[[1]])
head(fety.dredge.fixed.ordnum[[1]])

# with glmm from nlme package####
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


# Repeat with glmm from lme4 package####
# Anova Type II
require(lme4)
f.stats.lmer2 <- matrix(NA,11,20)
p.stats.lmer2 <- matrix(NA,11,20)

for(i in 1:18) {
  indices2$y <- indices2[,i]
  model <- lmer(y ~ age_class*samcam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1 + (1|field.ID) + (1|location), indices2)
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



