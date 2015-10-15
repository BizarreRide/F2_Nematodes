#§§§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Data Processing
# Calculate Biodiversity Indices
# Quentin Schorpp
# 06.08.2015
#§§§§§§§§§§§§§§§§§§§§§§§§§


# Indices ####
# Calculate and add biodiversity Indices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biodiv.fun <- function(data) data.frame(
  SR = rowSums(data >0),                              # families Richness
  rarefy = vegan::rarefy(data,90,se=F, MARGIN=1),     # rarefaction
  H = vegan::diversity(data, index="shannon"),        # Shannon entropy
  D = vegan::diversity(data, index="simpson"),        # simpson dominance
  J = vegan::diversity(data, index="shannon")/log(rowSums(data >0)),     # Pielou Evenness
  H1 = exp(vegan::diversity(data, index="shannon")))                     # Hill's N1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%