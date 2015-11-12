#§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Calculate Feeding Types
# Quentin Schorpp
# 07.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$

master <- read.csv("Data/F2_Nema_MastertableFam.csv", sep=";")
master$c.p <- as.factor(master$c.p)
master$FeedingType <- as.character(master$FeedingType)
master[29,"FeedingType"] <- "T"
master$FeedingType <- factor(master$FeedingType)

fety <- data %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>%
  group_by(ID.x, FeedingType) %>% 
  summarise(Value=sum(Value)) %>%
  spread(FeedingType,Value) 
fety <- fety[order(as.numeric(fety$ID.x)),]

