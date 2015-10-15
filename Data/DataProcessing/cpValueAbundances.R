#§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Calculate cp value abundances
# Quentin Schorpp
# 13.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$

master <- read.csv("Data/F2_Nema_MastertableFam.csv", sep=";")
master$c.p <- as.factor(master$c.p)

cpvalues <- data %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>%
  group_by(ID.x, c.p) %>% 
  summarise(Value=sum(Value)) %>%
  spread(c.p,Value) 
cpvalues <- cpvalues[order(as.numeric(cpvalues$ID.x)),]