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


ggplot(indices2, aes(x=age_class, y=N)) + 
  facet_grid(. ~ samcam) +
  stat_summary(fun.data = 'mean_sdl', geom = 'errorbar', width = 0.2, size = 1) +
  stat_summary(fun.y = mean, geom = 'bar', size = 1, color = 'black') +
  stat_summary(fun.y = mean, geom = 'line', size = 1, color = 'red') +
  coord_cartesian( y=c(0,150))
 

indices3 <- aggregate(indices, list(age_class=indices2$age_class, samcam=indices2$samcam), mean)
indices3.sd <- aggregate(indices, list(age_class=indices2$age_class, samcam=indices2$samcam), sd)

ggplot(indices3, aes(x=age_class, y=N)) + 
  facet_grid(. ~ samcam) +
  geom_bar(stat="identity") +
  geom_errorbar(data = indices3.sd, aes(y = indices3$N, ymin = indices3$N - N, ymax = indices3$N + N),
                width = 0.3, size = 1) +
  coord_cartesian( y=c(0,150))


se <- function(x) {
  return(c(ymin=(mean(x) - (sd(x)/sqrt(length(x)))), ymax=(mean(x) +
                                                             (sd(x)/sqrt(length(x))))))
}

gg_list <- list()
gg_list[[i]]  <-
  
  
for (i in 1:18) {
  nam <- paste("Plot",i, sep = ".")
  gg_list[[i]]  <- g <-  print(ggplot(indices2, aes(x = age_class, y = indices2[,i])) +
    stat_summary(fun.data = 'se', geom = 'errorbar', width = 0.2, size = 1) + 
    #stat_summary(fun.y = mean, geom = 'point', size = 3, color = 'red') +
    stat_summary(fun.y = mean, geom = 'bar', size = 1, color = 'gray20', fill="lightpink4", alpha=1/3) +
    facet_wrap(~ samcam) +
    ggtitle(names(indices2[i])) +
    theme_bw())
    assign(nam, g)
}

do.call(gridExtra::grid.arrange, gg_list)

Plot.1
Plot.2

Plot.4

gg_list[[2]]


i=6


g

