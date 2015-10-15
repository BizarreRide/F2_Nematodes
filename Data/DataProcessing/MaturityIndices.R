#$$$$$$$$$$$$$$$$$$$$$$$$$
# F2 Nematodes
# Calculate Indices in R
# Quentin Schorpp
# 07.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# take the values of a data frame and multiply them by a coefficient in column B of the mastertable, if their colnames matches the name of column A in the mastertable

master <- read.csv("Data/F2_Nema_MastertableFam.csv", sep=";")
master$c.p <- as.factor(master$c.p)




# Maturity Index ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N <- rowSums(data)

data.sum <- data %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family")

MI.sum <- data.sum %>% 
  filter(MI.PPI == "MI") %>% 
  group_by(ID.x) %>% 
  summarise(MI = sum(Value))

MI25.sum <- data.sum %>% 
  filter(MI.PPI == "MI", c.p %in% c(2:5)) %>% 
  group_by(ID.x) %>% 
  summarise(MI25 = sum(Value))

sigmaMI25.sum <- data.sum %>% 
  filter(c.p %in% c(2:5)) %>% 
  group_by(ID.x) %>% 
  summarise(sigmaMI25 = sum(Value))

PPI.sum <- data.sum %>% 
  filter(MI.PPI == "PPI") %>% 
  group_by(ID.x) %>% 
  summarise(PPI = sum(Value))

data.sum <- cbind(MI.sum, MI25.sum[,2], sigmaMI25.sum[,2], PPI.sum[,2],N)

data.sum$ID.x <- as.numeric(data.sum$ID.x)
data.sum <- data.sum[order(data.sum$ID.x),]
data.sum$ID.x <- rownames(data.sum)




# Maturity Index 
master$c.p <- as.numeric(master$c.p)

MI.matrix <- matrix(NA,nrow(data),ncol(data))
colnames(MI.matrix) <- colnames(data)

# c-p weighting of abundances
for (i in 1:ncol(data)){
  for (j in 1:nrow(data)){
    if (colnames(data)[i] == master[j,2]) {
      MI.matrix[,i] <- (data[,i]*master[j,4])/data.sum$MI  # multiply abundances by cp values and divide by sum of group members
    }
  }
}

MI <- MI.matrix %>% 
  as.data.frame() %>%
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  filter(MI.PPI == "MI") %>% 
  group_by(ID.x) %>% 
  summarise(MI = sum(Value))
MI$ID.x <- as.numeric(MI$ID.x)
MI <- MI[order(MI$ID.x),]

# Maturity Index cp 2-5
MI25.matrix <- matrix(NA,nrow(data),ncol(data))
colnames(MI25.matrix) <- colnames(data)

# c-p weighting of  abundances
for (i in 1:ncol(data)){
  for (j in 1:nrow(data)){
    if (colnames(data)[i] == master[j,2]) {
      MI25.matrix[,i] <- (data[,i]*master[j,4])/data.sum$MI25
    }
  }
}

MI25 <- MI25.matrix %>% 
  as.data.frame() %>%
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  filter(MI.PPI == "MI", c.p %in% c(2:5)) %>% 
  group_by(ID.x) %>% 
  summarise(MI25 = sum(Value))
MI25$ID.x <- as.numeric(MI25$ID.x)
MI25 <- MI25[order(MI25$ID.x),]


# Summed Maturity Index 
sigmaMI.matrix <- matrix(NA,nrow(data),ncol(data))
colnames(sigmaMI.matrix) <- colnames(data)

# c-p weighting of  abundances
for (i in 1:ncol(data)){
  for (j in 1:nrow(data)){
    if (colnames(data)[i] == master[j,2]) {
      sigmaMI.matrix[,i] <- (data[,i]*master[j,4])/data.sum$N
    }
  }
}

sigmaMI <- sigmaMI.matrix %>% 
  as.data.frame() %>%
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  group_by(ID.x) %>% 
  summarise(sigmaMI = sum(Value))
sigmaMI$ID.x <- as.numeric(sigmaMI$ID.x)
sigmaMI <- sigmaMI[order(sigmaMI$ID.x),]

# Summed Maturity Index for c-p 2-5
sigmaMI25.matrix <- matrix(NA,nrow(data),ncol(data))
colnames(sigmaMI25.matrix) <- colnames(data)

# c-p weighting of  abundances
for (i in 1:ncol(data)){
  for (j in 1:nrow(data)){
    if (colnames(data)[i] == master[j,2]) {
      sigmaMI25.matrix[,i] <- (data[,i]*master[j,4])/data.sum$sigmaMI25
    }
  }
}

sigmaMI25 <- sigmaMI25.matrix %>% 
  as.data.frame() %>%
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  filter(c.p %in% c(2:5)) %>%
  group_by(ID.x) %>% 
  summarise(sigmaMI25 = sum(Value))
sigmaMI25$ID.x <- as.numeric(sigmaMI25$ID.x)
sigmaMI25 <- sigmaMI25[order(sigmaMI25$ID.x),]



# Plant Parasite Index
PPI.matrix <- matrix(NA,nrow(data),ncol(data))
colnames(PPI.matrix) <- colnames(data)

# c-p weighting of abundances
for (i in 1:length(data)){
  for (j in 1:30){
    if (colnames(data)[i] == master[j,2]) {
      PPI.matrix[,i] <- (data[,i]*master[j,4])/PPI.sum$PPI
    }
  }
}

PPI <- PPI.matrix %>% 
  as.data.frame() %>%
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  filter(MI.PPI == "PPI") %>% 
  group_by(ID.x) %>% 
  summarise(PPI = sum(Value))
PPI$ID.x <- as.numeric(PPI$ID.x)
PPI <- PPI[order(PPI$ID.x),]

MaturityIndices <- data.frame(
  MI = MI[,2],
  PPI = PPI[,2],
  MI25 = MI25[,2],
  sigmaMI = sigmaMI[,2],
  sigmaMI25 = sigmaMI25[,2],
  PPIMI = PPI[,2]/MI[,2]
)

rm(MI.matrix, PPI.matrix, MI25.matrix, sigmaMI25.matrix, sigmaMI.matrix, MI.sum, PPI.sum, MI25.sum, sigmaMI25.sum)
rm(MI, MI25, sigmaMI, sigmaMI25, PPI)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Example from Thierry ####
# calculate summed abundances only for species that count for the MI
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# library(dplyr)
# library(tidyr)
# data(varespec, package = "vegan")
# attributes <- data.frame(
#   Species = c(colnames(varespec), "spec1", "spec2"),
#   Attribute = c(rep(c("MI", "PI"), c(14, 30)), "MI", "PI")
# )
# varespec %>% 
#   add_rownames("ID") %>% 
#   gather(Species, Value, -ID) %>% #convert to long format
#   inner_join(attributes, by = "Species") %>% 
#   filter(Attribute == "MI") %>% 
#   group_by(ID) %>% 
#   summarise(Total = sum(Value))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
