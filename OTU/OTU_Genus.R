library(DESeq2)
library(ggplot2)
library(dplyr)
library(edgeR)
library(patchwork)
library(vegan)
library(ape)
library(phangorn)
library(ggsci)
library(tidytree)
setwd('D:/PET/OTU/All_OTU')
# <- read.csv('All_OTU_Genus_Run.CSV')
#Genus_Run <- colnames(Genus_Run)
#Genus_Run <- as.data.frame(Genus_Run)
#Genus_Genus_Run <- Genus_Run
#write.csv(Genus_Run,'Genus_Run.csv')

Genus_Run <- read.csv('Genus_Run.csv')

PET <- read.csv('D:/PET/Abundance/PET_merged.csv')
PET_Genus <- merge(PET,Genus_Run,by='Run',all.y=T)
PET_Genus <- na.omit(PET_Genus)
PET_Genus <- PET_Genus[!duplicated(PET_Genus$Run),]

Genus <- read.csv('All_OTU_Genus.csv')
Genus_clade <- Genus$clade_name
xx <- NULL
for (i in 1:length(Genus_clade)) {
  a <- unlist(strsplit(Genus_clade[i],split = '\\|'))
  b <- a[length(a)]
  xx <- append(xx,b)
}

Genus_clade <- xx
Genus_Genus_Run <- colnames(Genus)
Genus_Genus_Run <- Genus_Genus_Run[-1]


Genus <- as.data.frame(apply(Genus,1,as.numeric))
Genus <- Genus[-1,]
colnames(Genus) <- Genus_clade
rownames(Genus) <- Genus_Genus_Run
Genus$Genus_Run <- Genus_Genus_Run
Genus <- na.omit(Genus)

Genus_df <- merge(Genus,PET_Genus,by='Genus_Run',all.y = T)
Genus_df <- na.omit(Genus_df)
table(Genus_df$Type)
Genus_dt <- Genus_df[,-c(1,1946:1948)]
Genus_info <- Genus_df[,c(1,1946:1948)]
Genus_dt_type <- Genus_df[,-c(1,1946,1947)]

result <- Genus_dt_type %>%
  group_by(Type) %>%
  summarise(across('g__Candidatus Babela':g__Methylacidimicrobium, sum, na.rm = TRUE))
result_type <- result$Type
result_genus <- colnames(result)
result_genus <- result_genus[-1]
result_t <- as.data.frame(t(result))
result_t <- result_t[-1,]
colnames(result_t) <- result_type
result_t <- as.data.frame(apply(result_t,2,as.numeric))
result_t$Sum <- rowSums(result_t)
rownames(result_t) <- result_genus
type_genus_sum <- colSums(result_t)

PET_Genus_landfill <- PET_Genus[PET_Genus$Type == 'landfill',]
median_landfill <- median(PET_Genus_landfill$Abundance)
PET_Genus_landfill <- PET_Genus_landfill %>%
  mutate(Group = ifelse(Abundance > median_landfill,"high","low")) # landfill = 171

PET_Genus_compost <- PET_Genus[PET_Genus$Type == 'compost',]
median_compost <- median(PET_Genus_compost$Abundance)
PET_Genus_compost <- PET_Genus_compost %>%
  mutate(Group = ifelse(Abundance > median_compost,"high","low")) # compost = 269

PET_Genus_marinesediment <- PET_Genus[PET_Genus$Type == 'marinesediment',]
median_sediment <- median(PET_Genus_marinesediment$Abundance)
PET_Genus_marinesediment <- PET_Genus_marinesediment %>%
  mutate(Group = ifelse(Abundance > median_sediment,"high","low")) # marinesediment = 197

PET_Genus_phyllosphere <- PET_Genus[PET_Genus$Type == 'phyllosphere',]
median_phyllosphere <- median(PET_Genus_phyllosphere$Abundance)
PET_Genus_phyllosphere <- PET_Genus_phyllosphere %>%
  mutate(Group = ifelse(Abundance > median_phyllosphere,"high","low")) # phyllosphere = 126

PET_Genus_rhizosphere <- PET_Genus[PET_Genus$Type == 'rhizosphere',]
median_rhizosphere <- median(PET_Genus_rhizosphere$Abundance)
PET_Genus_rhizosphere <- PET_Genus_rhizosphere %>%
  mutate(Group = ifelse(Abundance > median_rhizosphere,"high","low")) # rhizosphere = 190

table(PET_Genus$Type)
landfill <- merge(Genus,PET_Genus_landfill,by = 'Genus_Run',all.y = T)
compost <- merge(Genus,PET_Genus_compost,by='Genus_Run',all.y = T)
sediment <- merge(Genus,PET_Genus_marinesediment,by='Genus_Run',all.y = T)
phyllosphere <- merge(Genus,PET_Genus_phyllosphere,by='Genus_Run',all.y = T)
rhizosphere <- merge(Genus,PET_Genus_rhizosphere,by='Genus_Run',all.y = T)

#landfill
dt <- landfill[,-c(1,1946:1949)]
rownames(dt) <- landfill$Run
group <- landfill$Group
dt <- as.data.frame(t(dt))
dt <- dt*10000
dt <- round(dt)
group_file <- data.frame(Run = landfill$Run,Group = landfill$Group)
group_file$Group <- as.factor(group_file$Group)
nrow(group_file)
ncol(dt)


dds <- DESeq2::DESeqDataSetFromMatrix(countData = dt,
                                      colData = group_file,
                                      design = ~ Group)
dds_res <- DESeq2::DESeq(dds,sfType = 'poscounts')
resultsNames(dds_res)
res <- results(dds_res,tidy=T,format='DataFrame',contrast = c('Group',"high","low"))

DEG <- res
logFC_cutoff <- 2
DEG$change <- as.factor(ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
                               ifelse(DEG$log2FoldChange>logFC_cutoff,"UP","DOWN"),
                               'NOT'))
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up OTUs is ',nrow(DEG[DEG$change == 'UP',]),
                     '\nThe number of down OTUs is ',nrow(DEG[DEG$change == 'DOWN',]))
DEG <- na.omit(DEG)
landfill_DEG <- DEG
ggplot(data=DEG,aes(x=log2FoldChange,
                    y=-log10(pvalue),
                    color=change))+
  geom_point(alpha=0.8,size=3)+
  labs(x='log2 fold change')+ylab('-log10 pvalue')+
  ggtitle(this_title) + theme_bw(base_size = 20)+
  theme(plot.title = element_text(size=15,hjust = 0.5),)+
  scale_color_manual(values = c('#a121f0','#bebebe','#ffad21')) -> p1
p1

# compost
dt <- compost[,-c(1,1946:1949)]
rownames(dt) <- compost$Run
group <- compost$Group
dt <- as.data.frame(t(dt))
dt <- dt*10000
dt <- round(dt)
group_file <- data.frame(Run = compost$Run,Group = compost$Group)
group_file$Group <- as.factor(group_file$Group)
nrow(group_file)
ncol(dt)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = dt,
                                      colData = group_file,
                                      design = ~ Group)
dds_res <- DESeq2::DESeq(dds,sfType = 'poscounts')
resultsNames(dds_res)
res <- results(dds_res,tidy=T,format='DataFrame',contrast = c('Group',"high","low"))

DEG <- res
logFC_cutoff <- 2
DEG$change <- as.factor(ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
                               ifelse(DEG$log2FoldChange>logFC_cutoff,"UP","DOWN"),
                               'NOT'))
compost_DEG <- DEG
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up OTUs is ',nrow(DEG[DEG$change == 'UP',]),
                     '\nThe number of down OTUs is ',nrow(DEG[DEG$change == 'DOWN',]))
DEG <- na.omit(DEG)
ggplot(data=DEG,aes(x=log2FoldChange,
                    y=-log10(pvalue),
                    color=change))+
  geom_point(alpha=0.8,size=3)+
  labs(x='log2 fold change')+ylab('-log10 pvalue')+
  ggtitle(this_title) + theme_bw(base_size = 20)+
  theme(plot.title = element_text(size=15,hjust = 0.5),)+
  scale_color_manual(values = c('#a121f0','#bebebe','#ffad21')) -> p2
p2
pa <- p1+p2

# sediment
dt <- sediment[,-c(1,1946:1949)]
rownames(dt) <- sediment$Run
group <- sediment$Group
dt <- as.data.frame(t(dt))
dt <- dt*10000
dt <- round(dt)
group_file <- data.frame(Run = sediment$Run,Group = sediment$Group)
group_file$Group <- as.factor(group_file$Group)
nrow(group_file)
ncol(dt)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = dt,
                                      colData = group_file,
                                      design = ~ Group)
dds_res <- DESeq2::DESeq(dds,sfType = 'poscounts')
resultsNames(dds_res)
res <- results(dds_res,tidy=T,format='DataFrame',contrast = c('Group',"high","low"))

DEG <- res
logFC_cutoff <- 2
DEG$change <- as.factor(ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
                               ifelse(DEG$log2FoldChange>logFC_cutoff,"UP","DOWN"),
                               'NOT'))
sediment_DEG <- DEG
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up OTUs is ',nrow(DEG[DEG$change == 'UP',]),
                     '\nThe number of down OTUs is ',nrow(DEG[DEG$change == 'DOWN',]))
DEG <- na.omit(DEG)
ggplot(data=DEG,aes(x=log2FoldChange,
                    y=-log10(pvalue),
                    color=change))+
  geom_point(alpha=0.8,size=3)+
  labs(x='log2 fold change')+ylab('-log10 pvalue')+
  ggtitle(this_title) + theme_bw(base_size = 20)+
  theme(plot.title = element_text(size=15,hjust = 0.5),)+
  scale_color_manual(values = c('#a121f0','#bebebe','#ffad21')) -> p3
p3

#phyllosphere
dt <- phyllosphere[,-c(1,1946:1949)]
rownames(dt) <- phyllosphere$Run
group <- phyllosphere$Group
dt <- as.data.frame(t(dt))
dt <- dt*10000
dt <- round(dt)
group_file <- data.frame(Run = phyllosphere$Run,Group = phyllosphere$Group)
group_file$Group <- as.factor(group_file$Group)
nrow(group_file)
ncol(dt)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = dt,
                                      colData = group_file,
                                      design = ~ Group)
dds_res <- DESeq2::DESeq(dds,sfType = 'poscounts')
resultsNames(dds_res)
res <- results(dds_res,tidy=T,format='DataFrame',contrast = c('Group',"high","low"))

DEG <- res
logFC_cutoff <- 2
DEG$change <- as.factor(ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
                               ifelse(DEG$log2FoldChange>logFC_cutoff,"UP","DOWN"),
                               'NOT'))
phyllospere_DEG <- DEG
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up OTUs is ',nrow(DEG[DEG$change == 'UP',]),
                     '\nThe number of down OTUs is ',nrow(DEG[DEG$change == 'DOWN',]))
DEG <- na.omit(DEG)
ggplot(data=DEG,aes(x=log2FoldChange,
                    y=-log10(pvalue),
                    color=change))+
  geom_point(alpha=0.8,size=3)+
  labs(x='log2 fold change')+ylab('-log10 pvalue')+
  ggtitle(this_title) + theme_bw(base_size = 20)+
  theme(plot.title = element_text(size=15,hjust = 0.5),)+
  scale_color_manual(values = c('#a121f0','#bebebe','#ffad21')) -> p4
p4

p <- p1+p2+p3+p4
p

landfill_up <- landfill_DEG[landfill_DEG$change == 'UP',]
landfill_down <- landfill_DEG[landfill_DEG$change == 'DOWN',]

compost_up <- compost_DEG[compost_DEG$change =='UP',]
compost_down <- compost_DEG[compost_DEG$change == 'DOWN',]

sediment_up <- sediment_DEG[sediment_DEG$change =='UP',]
sediment_down <- sediment_DEG[sediment_DEG$change =='DOWN',]

phyllospere_up <- phyllospere_DEG[phyllospere_DEG$change =='UP',]
phyllospere_down <- phyllospere_DEG[phyllospere_DEG$change == 'DOWN',]

library(ggvenn)
up_list <- list(landfill=landfill_up$row,
                compost = compost_up$row,
                phyllosphere = phyllospere_up$row)
ggvenn(data=up_list)
down_list <- list(landfill=landfill_down$row,
                  compost = compost_down$row,
                  phyllosphere = phyllospere_down$row)
ggvenn(data=down_list)

median_Genus <- median(Genus_df$Abundance)
PET_Genus <- Genus_df %>%
  mutate(Group = ifelse(Abundance > median_Genus,"high","low"))
group <- PET_Genus$Group
group_file <- data.frame(Run = PET_Genus$Run,Group = PET_Genus$Group)
group_file$Group <- as.factor(group_file$Group)
trans_dt <- Genus_dt
rownames(trans_dt) <- PET_Genus$Run
distance <- vegdist(trans_dt,method = 'bray')
pcoa <- cmdscale(distance,k=2,eig = TRUE)
plot_data <- data.frame({pcoa$points})
eig <- pcoa$eig
data <- data.frame(group_file$Group,plot_data)  
names(data) <- c('Group','PCoA1','PCoA2')  
p11 <- ggplot(data=data,aes(x=PCoA1,y=PCoA2,color=Group))
p11 <- p11+geom_point(alpha = 1,size=2)+
  stat_ellipse(aes(fill=Group),type = 'norm',geom='polygon',alpha=0.2,color=NA)+
  labs(x=paste('PCoA1(',format(100*eig[1]/sum(eig),digits=4),"%)",sep=''),y=paste('PCoA2(',format(100*eig[1]/sum(eig),digits=4),"%)",sep = ''))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'))
p11

#lefse
#install.packages('microeco')
library(microeco)
feature_dt <- as.data.frame(t(Genus_dt))
colnames(feature_dt) <- group_file$Run

group_file <-  data.frame(Run = PET_Genus$Run,Group = PET_Genus$Group)
group_file$Group <- as.factor(group_file$Group)

sample_dt <- group_file
rownames(sample_dt) <- sample_dt$Run
tax_dt <- read.csv('Abundance_OTU_Genus_tax.csv')
tax_dt$Kingdom <- 'k__Bacteria'
tax_dt$Species <- 's__'
tax_dt <- tax_dt %>% select(Kingdom,phylum,class,order,family,genus,Species)
rownames(tax_dt) <- tax_dt$genus

tax_dt_phy <- as.data.frame(tax_dt$phylum)
colnames(tax_dt_phy) <- 'Phylum'

dataset <- microtable$new(sample_table = sample_dt,
                          otu_table = feature_dt,
                          tax_table = tax_dt)
dataset

lefse <- trans_diff$new(dataset = dataset,
                        method = 'lefse',
                        group = "Group",
                        alpha = 0.01,
                        lefse_subgroup = NULL)

lefse$plot_diff_bar(use_number = 1:20,
                    width = 0.8,
                    group_order = c('high','low'))+
  ggsci::scale_color_npg()+
  ggsci::scale_fill_npg()

lefse$plot_diff_cladogram(use_taxa_num = 200,
                          use_feature_num = 10,
                          clade_label_level = 10,
                          group_order = c('low','high'))

landfill_sig <- landfill_DEG[landfill_DEG$change != 'NOT',]
landfill_sig$Type <- 'landfill'
compost_sig <- compost_DEG[compost_DEG$change != 'NOT',]
compost_sig$Type <- 'compost'
sediment_sig <- sediment_DEG[sediment_DEG$change != 'NOT',]
sediment_sig$Type <- 'sediment'
phyllosphere_sig <- phyllospere_DEG[phyllospere_DEG$change != 'NOT',]
phyllosphere_sig$Type <- 'phyllosphere'

sig_dt <- rbind(landfill_sig,compost_sig,sediment_sig,phyllosphere_sig)
colnames(sig_dt) <- c('Genus','baseMean','log2FoldChange',
                      'lfc5E','stat','pvalue',
                      'padj','change','Type')
write.csv(sig_dt,'DESeq2_4Type_sig_Genus.csv',row.names = F)
sig_up <- sig_dt[sig_dt$change == 'UP',]
sig_down <- sig_dt[sig_dt$change == 'DOWN',]
sig_plot <- sig_dt[sig_dt$log2FoldChange > 4.8 | sig_dt$log2FoldChange < -5.1,]

ggplot(sig_plot,aes(Genus,log2FoldChange,fill = Type))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c("#69b3d2", "#7cb48f", "#f2a15c", "#a77ec0"))+
  theme_classic()+
  ylim(-8,7)+
  coord_flip() +
  geom_segment(aes(y=-8,yend= -8,x = 0, xend = 0))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line.y=element_line(linetype=1,color="black",linewidth=0.5))+
  theme(axis.ticks.y=element_line(color="black",size=0.1,lineend = 10))+
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.text.x = element_text(size = 11))+
  theme(text = element_text(size = 14))

