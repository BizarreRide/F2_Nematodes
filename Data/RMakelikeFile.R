###########################
# F2 Nematodes
# Load Data - Makelike File
# Quentin Schorpp
# 05.08.2015
###########################


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Community Data (Family level)
spe <- read.delim("Data/spe.txt")

indices <- spe[,which(names(spe)=="N"):47] # Nematode Indices - without gneral diversity indices (besides total abundance N)
fam <- spe[, which(names(spe)=="Aphelenchidae"):which(names(spe)=="Pratylenchidae")]
fety <- spe[, which(names(spe)=="carnivore"):which(names(spe)=="omnivore")]

env <- read.delim("Data/env.txt")
env1 <- env[16:45,]

env2 <- env[rep(seq_len(nrow(env)), each=3),-1]

counts <- read.delim("Data/counts.txt")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Calculate biodiversity Indices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# families Richness
indices$SR <- rowSums(fam >0)

# rarefaction
indices$rarefy <- vegan::rarefy(fam,90,se=F, MARGIN=1)


# Shannon entropy
indices$H <- vegan::diversity(fam, index="shannon")

# simpson dominance
indices$D <- vegan::diversity(fam, index="simpson")

# Pielou Evenness
indices$J <- indices$H/log(indices$SR)

# Hill's N1
indices$H1 <- exp(indices$H)

# camargo's diversity

# McIntosh dominance
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
