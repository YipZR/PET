library(dplyr)
library(reshape2)
library(ggplot2)
ALL_Run <- read.table('D:/PET/OTU/All_OTU/All_OTU_Phylum_Run.txt')

colnames(ALL_Run) <- ALL_Run[1,]
ALL_Run <- ALL_Run[-1,]
colnames(ALL_Run) <- c('Run','OTU_Run')
setwd('D:/PET/Abundance/')

PET <- read.csv('PET_merged.csv')

Run_PET <- merge(ALL_Run,PET,by.x = 'Run',by.y = 'Run',all.x = T)
Run_PET <- na.omit(Run_PET)
colnames(Run_PET) <- c('Run','OTU_Run','Abundance','Type')
table(Run_PET$Type)
bog_Run <- Run_PET[Run_PET$Type == 'Bog',]
compost_Run <- Run_PET[Run_PET$Type == 'compost',]
freshwater_Run <- Run_PET[Run_PET$Type == 'freshwater',]
hotspring_Run <- Run_PET[Run_PET$Type == 'Hotspring',]
landfill_Run <- Run_PET[Run_PET$Type == 'landfill',]
marine_Run <- Run_PET[Run_PET$Type == 'marine',]
marinesediment_Run <- Run_PET[Run_PET$Type == 'marinesediment',]
petroleum_Run <- Run_PET[Run_PET$Type == 'petroleum',]
phyllosphere_Run <- Run_PET[Run_PET$Type == 'phyllosphere',]
rhizosphere_Run <- Run_PET[Run_PET$Type == 'rhizosphere',]
volcano_Run <- Run_PET[Run_PET$Type == 'volcano',]
wetland_Run <- Run_PET[Run_PET$Type == 'wetland',]

phylum_dt <- read.csv('D:/PET/OTU/Run_OTU/OTU_Abundance_Phylum.csv')
clade_name <- phylum_dt$clade_name
Run_ID <- colnames(phylum_dt)
Run_ID <- Run_ID[-1]
phylum_dt2 <- t(phylum_dt)
phylum_dt2 <- as.data.frame(phylum_dt2[-1,])
phylum_dt2 <- na.omit(phylum_dt2)
colnames(phylum_dt2) <- clade_name
phylum_dt2$Run <- rownames(phylum_dt2)
#write.csv(phylum_dt2,'D:/PET/OTU/Run_OTU/phylum_dt2.csv',row.names = F)
phylum_dt2 <- read.csv('D:/PET/OTU/Run_OTU/phylum_dt2.csv')

phylum_ana <- merge(phylum_dt2,Run_PET,by='Run')

phylum_ana <- phylum_ana[!duplicated(phylum_ana$Run),]
phylum_ana_info <- phylum_ana[,c(1,49,50,51)]
#median(phylum_ana_info$Abundance)
#phylum_ana_info <- phylum_ana_info[phylum_ana_info$Abundance>=2.775761,]
phylum_ana <- merge(phylum_ana_info,phylum_ana,by='Run',all.x = T)

phylum_ana_dt <- phylum_ana[,-c(1:3,52:54)]
type_ana_dt <- phylum_ana_dt$Type.x
phylum_ana_dt <- phylum_ana_dt[,-1]
phy <- colnames(phylum_ana_dt)

xx <- NULL
for (i in 1:47) {
  a <- unlist(strsplit(phy[i],split = '\\.'))
  b <- a[length(a)]
  xx <- append(xx,b)
}
xx[13] <- 'p__Absconditabacteria'
xx[14] <- 'p__Bipolaricaulota'
xx[15] <- 'p__Cloacimonadota'
xx[16] <- 'p__Fervidibacteria'
xx[17] <- 'p__Omnitrophota'
xx[18] <- 'p__Saccharibacteria'

phylum_ana_dt <- as.data.frame(apply(phylum_ana_dt,2,as.numeric))
colnames(phylum_ana_dt) <- xx
phylum_ana_dt$type <- type_ana_dt
table(phylum_ana_dt$type)
xx
data_aggregated <- phylum_ana_dt %>%
  group_by(type) %>%
  summarise(across(p__Acidobacteriota:p__Verrucomicrobiota, sum, na.rm = TRUE))
data_aggregated
data_aggregated <- as.data.frame(t(data_aggregated))
colnames(data_aggregated) <- data_aggregated[1,]
data_aggregated <- data_aggregated[-1,]
data_aggregated <- as.data.frame(apply(data_aggregated,2,as.numeric))
rownames(data_aggregated) <- xx
data_aggregated$sum <- rowSums(data_aggregated)
data_aggregated_order <- data_aggregated[order(data_aggregated$sum,decreasing=T),]
data_aggregated_order <- data_aggregated_order[,-13]
data_aggregated_colsum <- colSums(data_aggregated_order)
data_aggregated_colsum

data_top10 <- data_aggregated_order[1:10,]
d <- NULL
for (i in 1:ncol(data_top10)) {
  x <- data_top10[,i]/data_aggregated_colsum[i]
  d <- cbind(d,x)  
}
colnames(d) <- colnames(data_top10)
rownames(d) <- rownames(data_top10)
col_sum_d <- colSums(d)
col_sum_d
data_others <- 1-col_sum_d
data_merge <- rbind(d,data_others)

kk <- rownames(data_merge)
kk[11] <- 'Others'
kk
rownames(data_merge) <- kk
bb <- colnames(data_merge)
data_merge <- as.data.frame(data_merge)

data_merge$phy <- kk
data_melt <- melt(data_merge,id='phy')

p <- ggplot(data_melt,aes(x=variable,y=100*value,fill=phy))+
  geom_col(position = 'stack', width = 0.6) +
  scale_y_continuous(expand=c(0, 0))+
  scale_fill_manual(values =  rev(c('#FFD700', 
                                    '#FF8C00', '#FF4500', '#40E0D0', '#00CED1',
                                    '#20B2AA', '#7FFF00', '#87CEEB', '#6495ED', 
                                    '#800080', '#CD5C5C'))) + 
  labs(x = NULL, y = 'Relative Abundance(%)') + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11),axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
  p
