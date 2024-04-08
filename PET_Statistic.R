#install.packages("tidyverse")
library(tidyverse)

dt_organism <- read.csv('D:/PET_Fasta/Ortho_result/Organism_PET.csv',header = F)
colnames(dt_organism) <- c('organism','Gene_ID')
write.csv(dt_organism,'D:/PET_Fasta/Ortho_result/Organism_PET.csv',row.names = F)
dt_Gene <- dt_organism$Gene_ID
dt_Organism <- as.data.frame(dt_organism$organism)
dt_statis <- map_dfr(dt_Organism,table)
dt_statis <- as.data.frame(t(dt_statis))
dt_statis <- cbind(rownames(dt_statis),dt_statis)
colnames(dt_statis) <- c('organism','count')
dt_statis <- dt_statis[order(dt_statis$count,decreasing = TRUE),]
SUM <- sum(dt_statis$count)
percentage <- dt_statis$count/SUM
dt_statis$Percentage = percentage
dim <- dim(dt_statis)
dim1 <- dim[1]
rownames(dt_statis) <- seq(from=1,to=dim1)

write.csv(dt_statis,'D:/PET_Fasta/Ortho_result/PET_Ortho_Statistic.csv',row.names = F)
