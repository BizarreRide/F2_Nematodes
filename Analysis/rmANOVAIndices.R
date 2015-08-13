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

# ANalysis of Nematode Indices

# Data processing
indices2 <- cbind(indices,env.fin)
indices2$samcam <- as.factor(indices2$samcam)
indices2$field.ID <- as.ordered(indices2$field.ID)
str(indices2)


# create interaction factor for interaction between age_class and samcam, 
# to see if age_class effects stay constant
env.fin$agsam <- with(env.fin, interaction(age_class,samcam))



# model formula: y ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1 + Error(field.ID)
# model formula2: y ~ age + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1 + Error(field.ID)


# Test for Heteroscedasticity:
for (i in 1:18) {
 x <-  fligner.test(indices[,i] ~ env.fin$agsam)
 print(x)
}


# Repeated measurments ANOVA 
# create Matrix to extract F-Values
f.stats <- matrix(NA,12,19)
# create Matrix to extract F-Values
p.stats <- matrix(NA,12,19)

# rmANOVAs
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



# TukeyHSD test

for (i in 1:18) {
  df <- get(paste("rmAOV",i,names(indices[i]), sep = "."))
  x <- TukeyHSD(df)
  print(cld(x))
  }

# No Tukey Tests for aov objects



# Repeat with glm
f.stats.lme <- matrix(NA,12,20)
p.stats.lme <- matrix(NA,12,20)

for(i in 1:18) {
  indices2$x <- indices2[,i]
  model <- lme(x ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID,data=indices2)
  print(summary(model))
  nam <- paste("lme",i,names(indices[i]), sep = ".")
  assign(nam, model)
  f.stats.lme[,i+2] <- round(Anova(model, type = "III")$"Chisq",2)
  p.stats.lme[,i+2] <- round(Anova(model, type = "III")$"Pr(>Chisq)",3)
}
f.stats.lme[,1] <- p.stats.lme[,1]  <- row.names(Anova(model, type="III"))
f.stats.lme[,2] <- p.stats.lme[,2]  <- Anova(model, type = "III")$"Df"
colnames(f.stats.lme)[3:20] <- colnames(p.stats.lme)[3:20] <- colnames(indices)
colnames(f.stats.lme)[1] <- colnames(p.stats.lme)[1] <- "Env"
colnames(f.stats.lme)[2] <- colnames(p.stats.lme)[2] <- "Df"

write.csv(p.stats.lme, file="p-Werte.csv")
write.csv(f.stats.lme, file="Chi2-Werte.csv")

Letters <- matrix(NA,10,18)
for(i in 1:18) {
  indices2$x <- indices2[,i]
  model <- lme(x ~ agsam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID, data=indices2)
  x <- summary(glht(model, linfct=mcp(agsam="Tukey")), test = adjusted(type = "fdr"))
  nam <- paste("glht",i,names(indices[i]), sep = ".")
  assign(nam, x)
  x <- cld(x)$mcletters
  xx <- x$Letters
  Letters[,i] <- xx
  row.names(Letters) <- levels(env.fin$agsam)
}
colnames(Letters) <- colnames(indices)
write.csv(Letters, file="Letters.csv")


model <- lme(NCR ~ age_class*samcam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID,data=indices2)
Anova(model)

model <- lme(NCR ~ agsam + pH + n + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID,data=indices2)
Anova(model)

indices2[,12]
indices2$NCR
# Plots

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


fety2 <- cbind(fety,env.fin)
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

pl2 = plyr::dlply(fety.melt, "variable", `%+%`, e1 = g2)

do.call(gridExtra::grid.arrange, pl2)



