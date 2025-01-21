# load library ####
library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(pheatmap)
library(ViSEAGO)
library(stringr)
library(ggforce)
library(viridis)

#load data ####
counts <- read.delim("./counts.txt", row.names=1) 
coldata <- read.delim("./coldata.txt",row.names = 1)
coldata$Group<-factor(coldata$Group)

# set color ####

col<-c('grey80','black','#FFAAAA','red')
pie(rep(1,4),col=col)

# Deseq####

dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = coldata,
                                    design = ~Group)
dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)


## PCA ####

#plotPCA(vsd,intgroup='Group')

pcadat<-plotPCA(vsd,intgroup="Group",returnData=T)

all(rownames(pcadat)==rownames(coldata)) #sanity check
pcadat$Group<-coldata$Group

percentVar.vsd<-round(100*attr(pcadat,"percentVar"))

pcadat$name<-str_split_i(pcadat$name,'_',1)
pcadat$name<-str_split_i(pcadat$name,'\\.',2)


pca<-ggplot(pcadat, aes(PC1, PC2, fill=Group)) +
  geom_point(size=5,pch=21,stroke=0.5)+ 
  geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd[2],"% variance")) + 
  #stat_ellipse()+
  #geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  #scale_fill_manual(values = color4liver)+
  #scale_color_manual(values = color4liver)+
  coord_fixed()

dev.new()
pdf("pca.pdf")
pca
dev.off()

## PCA-revised ####
vsd_sub<-as.data.frame(assay(vsd))
vsd_sub<-dplyr::select(vsd_sub,include)

coldata_sub<-coldata %>%
 subset(rownames(.) %in% include)

pca2<-prcomp(t(vsd_sub))
pca2_df<-as.data.frame(pca2$x)
pca2_df$group<-coldata_sub$Group
pca2.proportionvariances <- round((pca2$sdev^2) / (sum(pca2$sdev^2))*100,2)
pca2.proportionvariances <- paste(colnames(pca2_df),"(",paste(as.character(pca2.proportionvariances),"%",")", sep=""))

pca2_plot<-ggplot(pca2_df, aes(PC1, PC2, fill=group)) +
  geom_point(size=5,pch=21,stroke=0.5)+ 
  #geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd[2],"% variance")) + 
  #stat_ellipse()+
  #geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)+
  coord_fixed()
  
dev.new()
pdf("pca2_revised.pdf")
pca2_plot
dev.off()

  
# DEGs ####
## PAO1 vs Ctrl ####
res.PAO1vctrl <- as.data.frame(results(dds, contrast=c('Group','PAO1','Ctrl'),alpha=0.05)) 
res.PAO1vctrl<-res.PAO1vctrl[order(res.PAO1vctrl$padj),]
sig.PAO1vctrl<-res.PAO1vctrl[which(res.PAO1vctrl$padj<0.05 & 
                                     (abs(res.PAO1vctrl$log2FoldChange)>1)),]

write.csv(sig.PAO1vctrl,"sig.PAO1vctrl.csv")
write.csv(res.PAO1vctrl,"res.PAO1vctrl.csv")

## Mockaso vs Ctrl ####
res.Mockasovctrl <- as.data.frame(results(dds, contrast=c('Group','Mockaso','Ctrl'),alpha=0.05)) 
res.Mockasovctrl<-res.Mockasovctrl[order(res.Mockasovctrl$padj),]
sig.Mockasovctrl<-res.Mockasovctrl[which(res.Mockasovctrl$padj<0.05 & 
                                     (abs(res.Mockasovctrl$log2FoldChange)>1)),]

write.csv(sig.Mockasovctrl,"sig.Mockasovctrl.csv")
write.csv(res.Mockasovctrl,"res.Mockasovctrl.csv")


## PAO1aso vs PAO1 ####
res.PAO1asovPAO1 <- as.data.frame(results(dds, contrast=c('Group','PAO1aso','PAO1'),alpha=0.05)) 
res.PAO1asovPAO1<-res.PAO1asovPAO1[order(res.PAO1asovPAO1$padj),]
sig.PAO1asovPAO1<-res.PAO1asovPAO1[which(res.PAO1asovPAO1$padj<0.05 & 
                                     (abs(res.PAO1asovPAO1$log2FoldChange)>1)),]

write.csv(sig.PAO1asovPAO1,"sig.PAO1asovPAO1.csv")
write.csv(res.PAO1asovPAO1,"res.PAO1asovPAO1.csv")


## Mockaso vs Ctrl ####
res.PAO1asovMockaso <- as.data.frame(results(dds, contrast=c('Group','PAO1aso','Mockaso'),alpha=0.05)) 
res.PAO1asovMockaso<-res.PAO1asovMockaso[order(res.PAO1asovMockaso$padj),]
sig.PAO1asovMockaso<-res.PAO1asovMockaso[which(res.PAO1asovMockaso$padj<0.05 & 
                                           (abs(res.PAO1asovMockaso$log2FoldChange)>1)),]

write.csv(sig.PAO1asovMockaso,"sig.PAO1asovMockaso.csv")
write.csv(res.PAO1asovMockaso,"res.PAO1asovMockaso.csv")

# Gene ID ####
Mmu.dataset<-useDataset('mmusculus_gene_ensembl',mart=useMart("ensembl",host = "https://may2024.archive.ensembl.org"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name'), 
               filters='external_gene_name',
               values=rownames(counts),mart=Mmu.dataset)
genesymbols <- tapply(Genemap$ensembl_gene_id,
                      Genemap$external_gene_name, 
                      paste, collapse="; ")


# Heatmap ####
DEG.all<-c(rownames(sig.PAO1vctrl),
           rownames(sig.Mockasovctrl),
           rownames(sig.PAO1asovMockaso),
           rownames(sig.PAO1asovPAO1))
DEG.all<-DEG.all[!duplicated(DEG.all)]

hm_mat<-assay(vsd)
hm_mat<-subset(hm_mat,rownames(hm_mat) %in% DEG.all)

coldat_hm<-coldata[order(coldata$Group),,drop=F]

hm_mat<-hm_mat[,rownames(coldat_hm)]
hm_mat<-as.data.frame(hm_mat)

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
break_hm = seq(-2, 2,length.out=100)

my_color_annotation<-list(Group= c('Ctrl'='grey80','Mockaso'='black',
                                         'PAO1'='#FFAAAA',
                                         'PAO1aso'='red'))

include<-colnames(hm_mat[,-c(3,4,8,9,11,12,18)]) # only include these
hm_mat<-dplyr::select(hm_mat,include)

hm<-pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows =T,cluster_cols =F,
         annotation_col = coldat_hm,
         #annotation_row = rowdat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45,cutree_rows = 12) 


hc<-hclust(as.dist(1-cor(t(hm_mat),method = "pearson")),method = "complete")
all(hc$order==hm$tree_row$order) #all TRUE
hm_cluster <- as.data.frame(cutree(tree = hc, k = 12)) #k must be the same as the number of cutree rows
table(hm_cluster)


order.row<-hm$tree_row$order
hm_mat <- hm_mat[order.row,]

hm_cluster<-hm_cluster[rownames(hm_mat),,drop=F]
colnames(hm_cluster)<-c('cluster')
all(rownames(hm_mat)==rownames(hm_cluster))
hm_cluster$cluster<-factor(hm_cluster$cluster)
rowdat_hm<-hm_cluster

#check heatmap
my_color_annotation<-list(Group= c('Ctrl'='grey80','Mockaso'='black',
                                   'PAO1'='#FFAAAA',
                                   'PAO1aso'='red'),
                          cluster=c('1'='#E41A1C','2'='#FDB462','3'='#4DAF4A',
                                    '4'='#FCCDE5','5'='#377EB8'))

pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm,
         annotation_row = rowdat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45) 

rowdat_hm<-hm_cluster %>% 
  mutate(cluster=case_when(cluster=='1'~'2',
                           cluster=='2'~'2',
                           cluster=='3'~'1',
                           cluster=='4'~'2',
                           cluster=='5'~'2',
                           cluster=='6'~'2',
                           cluster=='7'~'5',
                           cluster=='8'~'3',
                           cluster=='9'~'3',
                           cluster=='10'~'4',
                           cluster=='11'~'5',
                           cluster=='12'~'2'
  ))

rowdat_hm<-rowdat_hm %>%
  arrange(cluster)

hm_mat<-hm_mat[rownames(rowdat_hm),]


hm_reclus<-pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows =F,cluster_cols =F,
         annotation_col = coldat_hm,
         annotation_row = rowdat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45) 


dev.new()
pdf("hm.plot.pdf")
hm_reclus
dev.off()

#functional analysis ####

clus1<-subset(rowdat_hm,rowdat_hm$cluster=='1')
clus1$ensembl<-genesymbols[rownames(clus1)]
clus1$ensembl<-str_split_i(clus1$ensembl,pattern = ';',1)

clus2<-subset(rowdat_hm,rowdat_hm$cluster=='2')
clus2$ensembl<-genesymbols[rownames(clus2)]
clus2$ensembl<-str_split_i(clus2$ensembl,pattern = ';',1)

clus3<-subset(rowdat_hm,rowdat_hm$cluster=='3')
clus3$ensembl<-genesymbols[rownames(clus3)]
clus3$ensembl<-str_split_i(clus3$ensembl,pattern = ';',1)

clus4<-subset(rowdat_hm,rowdat_hm$cluster=='4')
clus4$ensembl<-genesymbols[rownames(clus4)]
clus4$ensembl<-str_split_i(clus4$ensembl,pattern = ';',1)

clus5<-subset(rowdat_hm,rowdat_hm$cluster=='5')
clus5$ensembl<-genesymbols[rownames(clus5)]
clus5$ensembl<-str_split_i(clus5$ensembl,pattern = ';',1)



expressed_genes<-as.data.frame(counts)
expressed_genes[expressed_genes=="0"]<-NA
expressed_genes<-expressed_genes[complete.cases(expressed_genes),]
expressed_genes$ensembl<-genesymbols[rownames(expressed_genes)]
expressed_genes$ensembl<-str_split_i(expressed_genes$ensembl,pattern = ';',1)
expressed_genes<-expressed_genes$ensembl
write.table(expressed_genes,'expressed_genes.txt')

background<-scan("expressed_genes.txt",
                 quiet=TRUE,
                 what="")

Ensembl<-ViSEAGO::Ensembl2GO(version='112')
ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "mmusculus_gene_ensembl",
  Ensembl)

#Clus1 
write.table(clus1$ensembl,'clus1.txt')

clus1<-scan(
  "clus1.txt",
  quiet=TRUE,
  what="")

GO_clus1<-ViSEAGO::create_topGOdata(
  geneSel=clus1,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus1<-topGO::runTest(
  GO_clus1,
  algorithm ="classic",
  statistic = "fisher")

clus1_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("GO_clus1","classic_clus1")))

clus1_res<-as.data.frame(clus1_sResults@data) %>%
  arrange(desc(`condition.-log10_pvalue` ))
write.csv(clus1_res,"clus1_res.csv")


#clus2 
write.table(clus2$ensembl,'clus2.txt')

clus2<-scan(
  "clus2.txt",
  quiet=TRUE,
  what="")

GO_clus2<-ViSEAGO::create_topGOdata(
  geneSel=clus2,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus2<-topGO::runTest(
  GO_clus2,
  algorithm ="classic",
  statistic = "fisher")

clus2_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("GO_clus2","classic_clus2")))

clus2_res<-as.data.frame(clus2_sResults@data) %>%
  arrange(desc(`condition.-log10_pvalue` ))
write.csv(clus2_res,"clus2_res.csv")



#clus3 
write.table(clus3$ensembl,'clus3.txt')

clus3<-scan(
  "clus3.txt",
  quiet=TRUE,
  what="")

GO_clus3<-ViSEAGO::create_topGOdata(
  geneSel=clus3,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus3<-topGO::runTest(
  GO_clus3,
  algorithm ="classic",
  statistic = "fisher")

clus3_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("GO_clus3","classic_clus3")))

clus3_res<-as.data.frame(clus3_sResults@data) %>%
  arrange(desc(`condition.-log10_pvalue` ))
write.csv(clus3_res,"clus3_res.csv")


#clus4 
write.table(clus4$ensembl,'clus4.txt')

clus4<-scan(
  "clus4.txt",
  quiet=TRUE,
  what="")

GO_clus4<-ViSEAGO::create_topGOdata(
  geneSel=clus4,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus4<-topGO::runTest(
  GO_clus4,
  algorithm ="classic",
  statistic = "fisher")

clus4_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("GO_clus4","classic_clus4")))

clus4_res<-as.data.frame(clus4_sResults@data) %>%
  arrange(desc(`condition.-log10_pvalue` ))
write.csv(clus4_res,"clus4_res.csv")


#clus5 
write.table(clus5$ensembl,'clus5.txt')

clus5<-scan(
  "clus5.txt",
  quiet=TRUE,
  what="")

GO_clus5<-ViSEAGO::create_topGOdata(
  geneSel=clus5,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus5<-topGO::runTest(
  GO_clus5,
  algorithm ="classic",
  statistic = "fisher")

clus5_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("GO_clus5","classic_clus5")))

clus5_res<-as.data.frame(clus5_sResults@data) %>%
  arrange(desc(`condition.-log10_pvalue` ))
write.csv(clus5_res,"clus5_res.csv")


#selected genes ####

selgene <- read.delim("./selgene.txt")
rownames(selgene)<-selgene$Gene
selgene$Gene<-NULL

selgene.other<-selgene %>%
  subset(Group %in% c('Fibrogenic','Oxidative')) %>%
  arrange(Pathway) %>%
  dplyr::select(Pathway)
  
selgene.interferon<-selgene %>%
  subset(Group %in% c('Interferon')) %>%
  arrange(Pathway) %>%
  dplyr::select(Pathway)

selgene.chemokine<-selgene %>%
  subset(Group %in% c('Chemokines')) %>%
  arrange(Pathway) %>%
  dplyr::select(Pathway)

selgene.cytokine<-selgene %>%
  subset(Group %in% c('Cytokines')) %>%
  arrange(Pathway) %>%
  dplyr::select(Pathway)

selgene.gap<-selgene %>%
  subset(Group %in% c('Gap')) %>%
  arrange(Pathway) %>%
  dplyr::select(Pathway)


df<-res.PAO1asovPAO1 %>%
  dplyr::select(log2FoldChange)

hm_color.2<- colorRampPalette(rev(brewer.pal(11, "PRGn")))(100)
break_hm = seq(-1, 1,length.out=100)

##Gap ####
df.gap<-df[rownames(selgene.gap),,drop=F]

hm_gap<-pheatmap(df.gap,scale='none',border_color = NA,color = hm_color.2,
                       show_rownames = T,show_colnames = F,
                       cluster_rows =F,cluster_cols =F,
                       #annotation_col = coldat_hm,
                       annotation_row = selgene.gap,
                       breaks = break_hm,
                       #annotation_colors = my_color_annotation,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       #gaps_row = 5,
                       cellwidth =  15,cellheight =  15) 

dev.new()
pdf('hm_Gap.pdf',width = 15,height = 15)
hm_gap
dev.off()





##cytokines ####
df.cytokines<-df[rownames(selgene.cytokine),,drop=F]

hm_cytokines<-pheatmap(df.cytokines,scale='none',border_color = NA,color = hm_color.2,
                        show_rownames = T,show_colnames = F,
                        cluster_rows =F,cluster_cols =F,
                        #annotation_col = coldat_hm,
                        annotation_row = selgene.cytokine,
                        breaks = break_hm,
                        #annotation_colors = my_color_annotation,
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        #gaps_row = 5,
                        cellwidth =  15,cellheight =  15) 

dev.new()
pdf('hm_cytokines.pdf',width = 15,height = 15)
hm_cytokines
dev.off()


##chemokines ####
df.chemokines<-df[rownames(selgene.chemokine),,drop=F]

hm_chemokines<-pheatmap(df.chemokines,scale='none',border_color = NA,color = hm_color.2,
                       show_rownames = T,show_colnames = F,
                       cluster_rows =F,cluster_cols =F,
                       #annotation_col = coldat_hm,
                       annotation_row = selgene.chemokine,
                       breaks = break_hm,
                       #annotation_colors = my_color_annotation,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       #gaps_row = 5,
                       cellwidth =  15,cellheight =  15) 

dev.new()
pdf('hm_chemokines.pdf',width = 15,height = 15)
hm_chemokines
dev.off()




##Interferon ####
df.interferon<-df[rownames(selgene.interferon),,drop=F]

hm_interferon<-pheatmap(df.interferon,scale='none',border_color = NA,color = hm_color.2,
                   show_rownames = T,show_colnames = F,
                   cluster_rows =F,cluster_cols =F,
                   #annotation_col = coldat_hm,
                   annotation_row = selgene.interferon,
                   breaks = break_hm,
                   #annotation_colors = my_color_annotation,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   gaps_row = 5,
                   cellwidth =  15,cellheight =  15) 

dev.new()
pdf('hm_interferon.pdf')
hm_interferon
dev.off()


##other functions ####
df.other<-df[rownames(selgene.other),,drop=F]

hm_other<-pheatmap(df.other,scale='none',border_color = NA,color = hm_color.2,
         show_rownames = T,show_colnames = F,
         cluster_rows =F,cluster_cols =F,
         #annotation_col = coldat_hm,
         annotation_row = selgene.other,
         breaks = break_hm,
         #annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         gaps_row = 13,
         cellwidth =  15,cellheight =  15) 

dev.new()
pdf('hm_other.pdf')
hm_other
dev.off()

# volcano plot ####

vol.PAO1vctrl<-res.PAO1vctrl %>%
  mutate(lab=case_when(log2FoldChange> 1 & padj<0.05 ~ 'Up',
                       log2FoldChange< -1 & padj<0.05 ~ 'Down',
                       T~NA))

volplot.PAO1vctrl<-ggplot(vol.PAO1vctrl,aes(x=log2FoldChange,y=-log10(padj),color=lab))+
  geom_point(shape=16,alpha=0.5,size=5)+
  scale_color_manual(values=c('blue','red','grey50'))+
  #scale_fill_manual(values=c('blue','red','grey50'))+
  geom_hline(yintercept = -log10(0.05),size=2,color='black',linetype='dashed')+
  geom_vline(xintercept =1,size=2,color='black',linetype='dashed' )+
  geom_vline(xintercept =-1,size=2,color='black',linetype='dashed' )+
  theme_light()

dev.new()
pdf('volplot.PAO1vctrl.pdf')
volplot.PAO1vctrl
dev.off()


vol.PAO1asovPAO1<-res.PAO1asovPAO1 %>%
  mutate(lab=case_when(log2FoldChange> 1 & padj<0.05 ~ 'Up',
                       log2FoldChange< -1 & padj<0.05 ~ 'Down',
                       T~NA))

volplot.PAO1asovPAO1<-ggplot(vol.PAO1asovPAO1,aes(x=log2FoldChange,y=-log10(padj),color=lab))+
  geom_point(shape=16,alpha=0.5,size=5)+
  scale_color_manual(values=c('blue','red','grey50'))+
  #scale_fill_manual(values=c('blue','red','grey50'))+
  geom_hline(yintercept = -log10(0.05),size=2,color='black',linetype='dashed')+
  geom_vline(xintercept =1,size=2,color='black',linetype='dashed' )+
  geom_vline(xintercept =-1,size=2,color='black',linetype='dashed' )+
  theme_light()

dev.new()
pdf('volplot.PAO1asovPAO1.pdf')
volplot.PAO1asovPAO1
dev.off()



vol.Mockasovctrl<-res.Mockasovctrl %>%
  mutate(lab=case_when(log2FoldChange> 1 & padj<0.05 ~ 'Up',
                       log2FoldChange< -1 & padj<0.05 ~ 'Down',
                       T~NA))

volplot.Mockasovctrl<-ggplot(vol.Mockasovctrl,aes(x=log2FoldChange,y=-log10(padj),color=lab))+
  geom_point(shape=16,alpha=0.5,size=5)+
  scale_color_manual(values=c('blue','red','grey50'))+
  #scale_fill_manual(values=c('blue','red','grey50'))+
  geom_hline(yintercept = -log10(0.05),size=2,color='black',linetype='dashed')+
  geom_vline(xintercept =1,size=2,color='black',linetype='dashed' )+
  geom_vline(xintercept =-1,size=2,color='black',linetype='dashed' )+
  theme_light()

dev.new()
pdf('volplot.Mockasovctrl.pdf')
volplot.Mockasovctrl
dev.off()




vol.PAO1asovMockaso<-res.PAO1asovMockaso %>%
  mutate(lab=case_when(log2FoldChange> 1 & padj<0.05 ~ 'Up',
                       log2FoldChange< -1 & padj<0.05 ~ 'Down',
                       T~NA))

volplot.PAO1asovMockaso<-ggplot(vol.PAO1asovMockaso,aes(x=log2FoldChange,y=-log10(padj),color=lab))+
  geom_point(shape=16,alpha=0.5,size=5)+
  scale_color_manual(values=c('blue','red','grey50'))+
  #scale_fill_manual(values=c('blue','red','grey50'))+
  geom_hline(yintercept = -log10(0.05),size=2,color='black',linetype='dashed')+
  geom_vline(xintercept =1,size=2,color='black',linetype='dashed' )+
  geom_vline(xintercept =-1,size=2,color='black',linetype='dashed' )+
  theme_light()

dev.new()
pdf('volplot.PAO1asovMockaso.pdf')
volplot.PAO1asovMockaso
dev.off()


