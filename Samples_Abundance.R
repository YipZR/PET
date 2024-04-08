library(ggplot2)
library(tidyverse)
rm(list = ls())
setwd('D:/PET/Abundance')

dt <- read.csv('Merged.csv')
dt <- na.omit(dt)
dt_1 <- count(dt,Type)
colnames(dt_1) <- c('Type','Count')
rm(p)
p = ggplot(dt,aes(x=Type , y=Abundance , color = Type)) +
  geom_boxplot()+
  coord_flip()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)
p

dt1 <- subset(dt,dt$Abundance>0)
dt1_1 <- count(dt1,Type)
colnames(dt1_1) <- c('Type','Count_0')

dt2 <- subset(dt,dt$Abundance>10)
dt2_1 <- count(dt2,Type)
colnames(dt2_1) <- c('Type','Count_10')

dt_add <- data.frame(Type = c('Bog','Hotspring','freshwater','marine'),
                     Count_10 = c(0,0,0,0))
dt2_1 <- rbind(dt2_1,dt_add)

dt_count <- merge(dt_1,dt1_1,by='Type')
dt_count <- merge(dt_count,dt2_1,by='Type')

dt_count$percent_0 <- dt_count$Count_0/dt_count$Count
dt_count$percent_10 <- dt_count$Count_10/dt_count$Count

write.csv(dt_count,'SamplesCounts.csv')
p2 <- ggplot(dt_count,aes(x = Type , y = percent_0 , fill = Type)) +
  geom_bar(stat="summary",position='dodge')+
  labs(x=NULL,y='Percent_GT0')+
  coord_flip()
#  ggtitle("") +
#  theme(plot.title = element_text(hjust = 0.5))
p2

p3 <- ggplot(dt_count,aes(x = Type , y = percent_10 , fill = Type)) +
  geom_bar(stat="summary",position='dodge')+
  labs(x=NULL,y='Percent_GT10')+
  coord_flip()
p3
