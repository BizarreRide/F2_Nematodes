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

env.sc <- env.fin[,sapply(env.fin,is.numeric)]
env.sc <- data.frame(scale(env.sc, center=TRUE, scale=TRUE))
env.sc <- cbind(env.fin[,!sapply(env.fin,is.numeric)],env.sc)

indices2 <- droplevels(cbind(indices,env.sc))
indices2$fungi <- fety$fungivore
indices2$bacfu <- fety$bacterivore + fety$fungivore

fety.av <- fety/indices$N
fety2 <- droplevels(cbind(fety.av,env.sc,N=indices$N))
fety2$bacfu <- fety2$bacterivore + fety2$fungivore
#fety2$age_class <- as.ordered(fety2$age_class)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Repeated measurments Analysis 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Anova Type III
f.stats.lmer3 <- matrix(NA,12,20)
p.stats.lmer3 <- matrix(NA,12,20)

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

# Special Case NCR ####
lme.12.NCR <- glmer(fungi/bacfu ~ age_class*samcam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1 + (1|field.ID),weights=bacfu,family="binomial",indices2, control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
#model <- glmer(y ~ age_class*samcam + pH + c + mc + ata1 + (1|field.ID),weights=N,family="binomial",fety2,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
f.stats.lmer3[,14] <- round(Anova(lme.12.NCR, type = "III")$"Chisq",2)
p.stats.lmer3[,14] <- round(Anova(lme.12.NCR, type = "III")$"Pr(>Chisq)",3)

glht.12.NCR <- glmer(fungi/bacfu ~ agsam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1 + (1|field.ID),weights=bacfu,family="binomial",indices2, control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
#glht.12.NCR <- glmer(fungi/bacfu ~ agsam + clay + (1|field.ID),weights=bacfu,family="binomial",indices2, control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)))
x <- summary(glht(glht.12.NCR, linfct=mcp(agsam="Tukey")), test = adjusted(type = "bonferroni"))
x <- cld(x)$mcletters
xx <- x$Letters
Letters[,12] <- xx

par(mar=c(5,5,8,2))
plot(cld(x))

require(MuMIn)
options(na.action = "na.fail")
dredge.2.MI <- dredge(lme.2.MI, trace=FALSE, rank="AICc")
head(dredge.2.MI)
dredge.12.NCR <- dredge(lme.12.NCR, trace=FALSE, rank="AICc")
head(dredge.12.NCR)
save(dredge.12.NCR,file="dredgeNCR.Rdata")

#write.csv(Letters, file="LettersFDR.csv")
#write.csv(Letters, file="LettersBonferroni.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Feeding Types ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


lm1 <- lm(fungivore ~ field.ID, fety2)
lm2 <- lm(fungivore ~ samcam, fety2)

plot(resid(lm1), resid(lm2))
cov(resid(lm1), resid(lm2))

