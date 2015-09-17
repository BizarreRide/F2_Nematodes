#%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2 Nematodes
# EDA
# Quentin Schorpp
# 06.08.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/RMakeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#***********************************************Environmental Data**************************************************************************####

# Cross Tables to reveal the nested structure and/or balanced designs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tapply(indices$N,list(env1$age_class, env1$samcam),length)
tapply(indices$N,list(env1$field.ID,env1$age_class),length) # field.ID is nested in age_class
tapply(indices$N,list(env1$field.ID, env1$samcam),length)   # field.ID is partly nested in samcam (Maize field!)

tapply(indices$N,list(env1$crop, env1$samcam),length)
tapply(indices$N,list(env1$field.ID, env1$crop),length)
tapply(indices$N,list(env1$age_class, env1$crop),length)   # crop and age_class are correlated


tapply(indices$N,list(env1$field.ID, env1$age_class),length)
tapply(indices$N,list(env1$field.ID, env1$crop),length)
tapply(indices$N,list(env1$field.ID, env1$n),length)
tapply(indices$N,list(env1$field.ID, env1$c),length)
tapply(indices$N,list(env1$field.ID, env1$clay),length)

tapply(indices$N,list(env1$field.ID, env1$pH),length)
tapply(indices$N,list(env1$field.ID, env1$intensity),length)
tapply(indices$N,list(env1$field.ID, env1$fertilisation),length)
tapply(indices$N,list(env1$field.ID, env1$ata1),length)
tapply(indices$N,list(env1$field.ID, env1$ata2),length)
tapply(indices$N,list(env1$field.ID, env1$hum2),length)
tapply(indices$N,list(env1$field.ID, env1$prec1),length)

# also possible with 
# table(env1$field.ID, env1$prec1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Pair Plots to check variable correlations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Tes set
env.num <- env.fin[sapply(env.fin,is.numeric)]

# Pearson r linear correlation among environmentla variables
env.pearson <- cor(env.num)
round(env.pearson,2)

# Reorder the variables prior to plotting
env.o <- gclus::order.single(env.pearson)

# Pearson Correlation matrix
op <- par(mfrow=c(1,1), pty="s")
pairs(env.fin, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist, main="Pearson Correlation Matrix")
par(op)

# Kendall Correlation matrix
env.ken <- cor(env.num, method="kendall")

# Reorder the variables prior to plotting
env.o <- gclus::order.single(env.ken)

# Kendall Correlation matrix
op <- par(mfrow=c(1,1), pty="s")
pairs(env.fin, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist, main="Kendall Correlation Matrix")
par(op)

# C and N together is not valid
# Age_class and intensity, is doubtful
# Clay and age_class is ok, if the treshhold is above 0.5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





# Creating a grouped Data Object ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data <- cbind(indices, env.fin)

data <- groupedData(PPI ~ age_class | field.ID, data)
str(data)
#Purpose??
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#***********************************************Families/response Data**********************************************************************####


# Here I Also try to figure out, if averaging over repeated measurements confounds my results.

# Question: In which direction change the nematode species from 2012 to 2013?
# Are there general trends?

# Change data ####
# Change (slope) between two samples from year 1 and year 2. X1=1, X2=2  <- nominator is always 1, hence the slope is the change in y ( i.e. the difference between y2 - y1)
fam.slope <- fam[!env1$crop=="Maize",]

fam.slope <- fam.slope[13:24,] - fam.slope[1:12,]

par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.slope[,i] ~ env1$age_class[1:12], main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode rel. abundance",outer=T)
}

par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.slope[,i], main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode rel. abundance",outer=T)
}



# Aphelenchidae and Rhabditidae seem to periodically increase and decrease

# Aphelenchidae have a similar but less pronounced periodicity to APhelenchidae

# Qudsianematidae seem to have the antagonistic periodicity to Rhabdititdae

# Plectidae and Mononchidae increase more an more towards SP_O Stages

# Mononchidae seem to increae in general

# Hoplolaimidae decrease in SP_old stages

# Telotylenchidae increase in SP_Y Stages

# Thornenematidae increase in SP_O Stages

# Tylenchidae increase in SP_O stages

# Telotylenchidae and Hoplolaimidae seem antagonistic


# Averaged Data ####

par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.av[,i] ~ env.av$age, main=names(fam.av)[i])
  abline(h=0,col="red")
  title(main="Mean Nematode rel. abundance",outer=T)
}

par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.av[,i], main=names(fam.av)[i])
  abline(h=0,col="red")
  title(main="Mean Nematode rel. abundance",outer=T)
}


# Qudsianematiden zeigen zwischen den jahren bis auf SP_I2 Flächen, wo eine Zunahme zu erkennen ist, immer eine abnahme. In SP_I2 Flächen sind ihre Anteile außerdem am Höchsten

# There is a rise and Fall of Aphelenchidae, Aphelenchoididae, Panagrolaimidae, Cephalobidae, Rhabditidae and Mononchidae with peaking rel. abundances in SP_I1
# And a bit in Alaimidae, Discolaimidae and Qudsianematidae

# SP_Y seem to have the highest proportion of carnivorous nematodes

# SP_O stages semm characterized by high Hoplolaimidae, high thornenematidae, low Telotylenchidae and low Tylenchidae

# Nordiidae were only found in maize fields,

# Diphterophoridae make a difference between Maize/SP_Y and older SP stages

# Pratylenchiden kommen nur in Mais und in SP_I2 stages vor. 


# Wie müssen die Changes in den mean abundances berücksichtigt werden?,
# Denkbar wäre eine fehlinterpretation, wenn man einen kontinuielichen Anstieg zu alten Stufen in den Means sieht, aber die changes zwischen den Jahren ausschließlich negativ sind.
# das wäre ein Wiederspruch.

# Umgekehrt wäre jdeoch ein kontinuierlicher Trend durch durchweg positive changes bestätigt


par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam[,i] ~ env1$age, main=names(fam.av)[i])
  abline(h=0,col="red")
  title(main="True Nematode rel. Abundance",outer=T)
}

# Sowohl für Aphelenchiden, Mononchiden als auch für Cephalobiden ist der Trend zwischen den beiden Jahren immer eine Zunahme. Dennoch haben beide NICHT ihre höchsten Anteile in den SP_O stages.
# Es ist jedoch nicht auszuschlißen, dass sie beide kontinuierlich zunehmen, ihre Anteile jedoch auf Kosten anderer Anteile sinken.
# Hier spielt es eine entscheidende Rolle, dass immer nur 100 Individuen aus der Population bestimmt wurden.

# Hier ist eine Fehlinterpretation der mean daten zu vermuten, da sie ein Rise and fall muster zeigen, wo eigentlich ein kontinuierliche Zuwachsrate ist.

par(mfrow=c(1,1))
  boxplot(counts$counts.av ~ env2$samcam+env2$age_class, main="total abundance")
  abline(h=0,col="red")

# Von Frühjahr zu herbst nehmen die abundanzen ab

# guckt man sich die Frühjahre an, würde man sagen, man hat eine stetige Abnahme bis zu den ältesten Stufen

# guckt man sich die herbste an würde man sagen, man hat eine leichte Zunahme, mit Möglichen periodischen Schwankungen, da sich Zu und Abnahmen zwischen
# den Jahren abwechseln

# Also ein linearer trend für age_class und ein (CO)Sinuskurviger für samcam

par(mfrow=c(1,1))
boxplot(counts[!env2$samcam==1,]$counts.av ~ env2[!env2$samcam==1,]$samcam + env2[!env2$samcam==1,]$age_class, main="total abundance")
abline(h=0,col="red")

par(mfrow=c(1,1))
boxplot(counts[!env2$samcam==1,]$counts.av ~ env2[!env2$samcam==1,]$age, main="total abundance")
abline(h=0,col="red")


# Upscaled data ####
par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.usc[,i] ~ env1$age_class, main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode rel. abundance",outer=T)
}

par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.usc[,i]~ env1$age, main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode rel. abundance",outer=T)
}


# Growth rates
fam.usc.slope <- fam.usc[!env1$crop=="Maize",]

fam.usc.slope <- fam.usc.slope[13:24,] - fam.usc.slope[1:12,]

fam.usc.gr <- fam.usc.slope/fam.usc[1:12,]
fam.usc.gr[is.na(fam.usc.gr)] <- 0
fam.usc.gr[mapply(is.infinite, fam.usc.gr)]  <- 0

par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.usc.gr[,i] ~ env1$age_class[1:12], main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode rel. abundance",outer=T)
}

par(mfrow=c(5,5))
for (i in 1:25) {
  boxplot(fam.usc.gr[,i], main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode rel. abundance",outer=T)
}


