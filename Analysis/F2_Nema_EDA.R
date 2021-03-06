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

tapply(biodiv$N,list(env1$age_class, env1$samcam),length)
tapply(biodiv$N,list(env1$field.ID,env1$age_class),length) # field.ID is nested in age_class
tapply(biodiv$N,list(env1$field.ID, env1$samcam),length)   # field.ID is partly nested in samcam (Maize field!)

tapply(biodiv$N,list(env1$crop, env1$samcam),length)
tapply(biodiv$N,list(env1$field.ID, env1$crop),length)
tapply(biodiv$N,list(env1$age_class, env1$crop),length)   # crop and age_class are correlated


tapply(biodiv$N,list(env1$field.ID, env1$n),length)       # C and n was only measured once for all silphie fields
tapply(biodiv$N,list(env1$field.ID, env1$c),length)       # Only Maize fields have individual measurements
tapply(biodiv$N,list(env1$field.ID, env1$clay),length)    # similar for clay
                                                          # c, n and grain size is probably more useful in Analysis of averaged data

        
tapply(biodiv$N,list(env1$field.ID, env1$intensity),length)
tapply(biodiv$N,list(env1$field.ID, env1$fertilisation),length) # fertilisation should probably be treated as factor

tapply(biodiv$N,list(env1$field.ID, env1$pH),length)
tapply(biodiv$N,list(env1$field.ID, env1$mc),length)
tapply(biodiv$N,list(env1$field.ID, env1$ata1),length)    
tapply(biodiv$N,list(env1$field.ID, env1$ata2),length)
tapply(biodiv$N,list(env1$field.ID, env1$hum2),length)
tapply(biodiv$N,list(env1$field.ID, env1$prec1),length)     # climate data and soil chemistry is sample specific

# also possible with 
# table(env1$field.ID, env1$prec1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Pair Plots to check variable correlations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Test set
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
fam.slope <- fam.org[!env1$crop=="Maize",]

fam.slope <- fam.slope[13:24,] - fam.slope[1:12,]

par(mfrow=c(5,5),
    mar=c(2,2,2,2),
    oma=c(0.5,0.5,2,0.5))

for (i in 1:25) {
  boxplot(fam.slope[,i] ~ env1$age_class[1:12], main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode dominance",outer=T)
}

for (i in 1:25) {
  boxplot(fam.slope[,i], main=names(fam.slope)[i])
  abline(h=0,col="red")
  title(main="Change in Nematode dominance",outer=T)
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


#***********************************************Indices**********************************************************************####

x11()
par(mfrow=c(4,5))

for(i in 1:18) {
  boxplot(indices[,i] ~ interaction(env1$age_class,env1$samcam), indices)
  title(colnames(indices)[i])
}

x11()
par(mfrow=c(4,5))

for(i in 1:18) {
  scatter.smooth(env1$age, indices[,i])
  title(colnames(indices)[i])
}

x11()
par(mfrow=c(3,2))

for(i in 1:5) {
  scatter.smooth(env1$age, fety[,i])
  title(colnames(fety)[i])
}





indices.av <- aggregate(indices,list(env1$field.ID),mean)
fety.av <- aggregate(fety,list(env1$field.ID),mean)

x11()
par(mfrow=c(4,5))

for(i in 2:19) {
  boxplot(indices.av[,i] ~ env.av$age_class, indices.av)
  title(colnames(indices.av)[i])
}


x11()
par(mfrow=c(2,3))

for(i in 2:6) {
  boxplot(fety.av[,i] ~ env.av$age_class, fety.av)
  title(colnames(fety.av)[i])
}

#***********************************************Feeding types**********************************************************************####

boxplot(fety$fungivore ~ interaction(env1$age_class,env1$samcam))

boxplot(fety$fungivore~env1$field.ID)
boxplot(fam.rel$Hoplolaimidae~env.fin$agsam)
boxplot(fam$Telotylenchidae~env.fin$agsam)
boxplot(fam$Tylenchidae~env.fin$agsam)
boxplot(fam$env.fin$agsam)
boxplot(indices$NCR~env1$field.ID)
boxplot(indices$NCR~env.fin$agsam)

boxplot(indices$PPI~env.fin$agsam)
boxplot(indices$MI~env.fin$agsam)
boxplot(indices$MI2.5~env.fin$agsam)
boxplot(indices$summedMI~env.fin$agsam)
boxplot(indices$MI~env.fin$age_class)

plot(indices2$SI, indices2$EI)

hefubac <- subset(fety, select = c("herbivore","fungivore","bacterivore")) 
ade4::triangle.plot(hefubac)


#***********************************************Circle Diagrams**********************************************************************####
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") # Load original datasets nema, env, counts and define factors, etc.
rm(counts.env,counts.org,counts3SC,env.org,nema) # remove unnecessary datasets

source("Data/DataProcessing/EnvDataProcessing.R") # mainly slice and categorize, extract orthogonal variables in subset env.fin
source("Data/DataProcessing/FamDatProcessing.R") # upscaled (counts * relative abundance in subset of 100 Ind.), presene absence data and bioiv indices on upscales data 
biodiv <- biodiv.fun(round(fam.usc,0))

source("Data/DataProcessing/AverageData.R") 
fam.av.usc2 <- (fam.av.org/rowSums(fam.av.org))*counts.av$counts
#????????????????


# A Families without samcam,...
# nee upscaled abundances....

fam.av.usc1 <- round(fam.av.usc[,c("Psilenchidae", "Tylenchidae", "Criconematidae", "Telotylenchidae", "Hoplolaimidae", "Pratylenchidae")])
                                #, "Hoplolaimidae", "Cephalobidae", "Plectidae", "Telotylenchidae", "Rhabditidae", "Aporcelaimidae", "Aphelenchoididae", "Panagrolaimidae")],0))

df.fam <- cbind(age_class=env.av$age_class,
                N=rowSums(round(fam.av.usc1,0)),
                fam.av.usc1)

nema3 <- aggregate(. ~ age_class, df.fam, mean)
nema3 <- reshape2::melt(nema3, id.vars=c(1,2))

require(ggplot2)
Ne3 <- ggplot(nema3, aes(x=N/2, y = value, fill = variable, width = N)) +
  geom_bar(position="fill", stat="identity", col="black") + 
  facet_grid(.~ age_class) + 
  coord_polar("y") + mytheme + theme(legend.position="bottom") 
Ne3


# B Feeding types

data=fam.av.usc
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae

df.fety <- round(fety[,-c(1:4,5,9,11)],0) # fungivores2
df.fety <- round(fety[,-c(1:4,9,10,12)],0) # herbivores2
df.fety <- round(fety[,-c(1:4,11,12)],0) # Tylenchidae
df.fety <- cbind(age_class=env.av$age_class,
                 N = rowSums(round(df.fety,0)),
                 df.fety)

nema3 <- aggregate(. ~ age_class, df.fety, mean)
nema3 <- reshape2::melt(nema3, id.vars=c(1,2))

require(ggplot2)
Ne3 <- ggplot(nema3, aes(x=N/2, y = value, fill = variable, width = N)) +
  geom_bar(position="fill", stat="identity") + 
  facet_grid(.~ age_class) + 
  coord_polar("y") + mytheme + theme(legend.position="bottom") 
Ne3

