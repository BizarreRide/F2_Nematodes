
#$$$$$$$$$$$$$$$$$$$$$$$$$
# F2 Nematodes
# Calculate Indices in R
# Quentin Schorpp
# 07.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# take the values of a data frame and multiply them by a coefficient in column B of the mastertable, if their colnames matches the name of column A in the mastertable

master <- read.csv("Data/F2_Nema_MastertableFam.csv", sep=";")

master$bwt <- 0
master$bwt[master$FeedingType %in% c(2,3) & master$c.p==2] <- 0.8

master$ewt <- 0
master$ewt[master$FeedingType==2 & master$c.p==2] <- 0.8
master$ewt[master$FeedingType==3 & master$c.p==1] <- 3.2

master$swt <- 0
master$swt[master$FeedingType == 5 & master$c.p==4] <- 3.2
master$swt[master$FeedingType == 5 & master$c.p==5] <- 5
master$swt[master$FeedingType == 8 & master$c.p==4] <- 3.2
master$swt[master$FeedingType == 8 & master$c.p==5] <- 5
master$swt[master$FeedingType == 2 & master$c.p==3] <- 1.8
master$swt[master$FeedingType == 2 & master$c.p==4] <- 3.2
master$swt[master$FeedingType == 5 & master$c.p==4] <- 3.2

str(master)





data(varespec, package = "vegan")
attributes <- data.frame(
  Species = c(colnames(varespec), "spec1", "spec2"),
  Attribute = c(rep(c("MI", "PI"), c(14, 30)), "MI", "PI"),
  index1=c(rep(c(1,2,3,4),each=11),3,3)
)

attributes$new[attributes$Attribute=="MI" & attributes$index1==c(2,3)] <- 0.8



MI.matrix <- matrix(NA,30,25)

for (i in 1:25){
  for (j in 1:30){
    if (colnames(fam)[i] == master[j,2]) {
      MI.matrix[,i] <- fam[,i]*master[j,4]/
    }
  }
}

for (i in 1:25){
  for (j in 1:30){
    if (colnames(fam)[i] == master[j,2]) {
      if(master[j,5] == "MI") {
        vector1 <- apply(fam,1,sum)
      }
    }
  }}

vector1 - apply(fam,1,sum)
sum(fam)


tfam <- data.frame(t(fam))

for (i in 1:25){
  for (j in 1:30){
44444444444444444444447if (rownames(tfam)[i]==rownames(master)[j]){
      tfam$mat<- master[,5]
}}}




# Species data from vegan package:
data(varespec)

# create attributes table
attributes <- matrix(NA, length(varespec), 2)
attributes[,1] <- colnames(varespec)
attributes[,2] <- c(rep("MI",14),rep("PI",30))

# add species to the attribute table
x <- c("spec1","MI")
y <- c("spec2","PI")

attributes <- rbind(attributes, x, y)
row.names(attributes) <- c(1:46)

# calculate rowsums only for species contained in the attributes table 
# and having the entry "MI" in the attributes table
for (i in 1:44){
  for (j in 1:46){
    if ((colnames(fam)[i] == master[j,1]) & (master[j,2] == "MI")) {
      apply(varspec,1,sum)
    }
  }}


library(dplyr)
library(tidyr)

MI.sum <- fam %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  filter(MI.PPI == "MI") %>% 
  group_by(ID.x) %>% 
  summarise(Total = sum(Value))

MI.sum$ID.x <- as.numeric(MI.sum$ID.x)
MI.sum <- MI.sum[order(MI.sum$ID.x, decreasing = FALSE),]

library(dplyr)
library(tidyr)
data(varespec, package = "vegan")
attributes <- data.frame(
  Species = c(colnames(varespec), "spec1", "spec2"),
  Attribute = c(rep(c("MI", "PI"), c(14, 30)), "MI", "PI")
)
varespec %>% 
  add_rownames("ID") %>% 
  gather(Species, Value, -ID) %>% #convert to long format
  inner_join(attributes, by = "Species") %>% 
  filter(Attribute == "MI") %>% 
  group_by(ID) %>% 
  summarise(Total = sum(Value))


gather(fam, Species,orks,-ID)

fam <- add_rownames(fam,"ID")

fam <- gather(fam, Species,orks,-ID)

filter(fam, Species=="Pratylenchidae")
