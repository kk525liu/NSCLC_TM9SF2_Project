rm(list = ls())
gc()

if(!dir.exists("./result")) dir.create("./result")

library(Seurat)
library(qs)
library(SCP)
library(ggsci)

############################################################
## All cell


sce <- readRDS( "./rds/sce_harmony.rds"  ) 

# Plot 1: Ptprc FeaturePlot
p1 <- FeatureDimPlot( sce , features = c( 'Ptprc' ))
ggsave(filename = "./result/FeatureDimPlot_Ptprc.pdf", plot = p1, width = 6, height = 5)

# Plot 2: CellDimPlot
p2 <- CellDimPlot(    
  srt = sce, group.by = c("celltype"),label = T, label_insitu = T , label.size = 3 , 
  reduction = "umap", theme_use = "theme_blank"  )
ggsave(filename = "./result/CellDimPlot_AllCell.pdf", plot = p2, width = 7, height = 6)

# Plot 3: CellStatPlot
p3 <- CellStatPlot(srt = sce, stat.by = "celltype",  group.by = "group", label = F  , label.size = 5 )
ggsave(filename = "./result/CellStatPlot_AllCell.pdf", plot = p3, width = 6, height = 5)

#############---------------------------------------------------------------------------------###################
#### Neutrophil

Neutrophil <- readRDS( './rds/Neutrophil_harmony.rds'  )
Neutrophil$celltype2 <- Neutrophil$celltype 

new.cluster.ids <- c( "Neutrophil_C0",
                      "Neutrophil_C1", "Neutrophil_C2","Neutrophil_C3","Neutrophil_C4","Neutrophil_C5",
                      "Neutrophil_C6")

Neutrophil@active.ident<- as.factor( Neutrophil$celltype2  )
names(new.cluster.ids) <- levels(Neutrophil)
Neutrophil<- RenameIdents(Neutrophil, new.cluster.ids)
Neutrophil$celltype<-Neutrophil@active.ident

Neutrophil$celltype <- factor( Neutrophil$celltype , levels = c( 
  "Neutrophil_C0",
  "Neutrophil_C1", "Neutrophil_C2","Neutrophil_C3","Neutrophil_C4","Neutrophil_C5",
  "Neutrophil_C6"  )  ) 

# Plot 4: Neutrophil CellDimPlot
p4 <- CellDimPlot(  
  srt = Neutrophil,  group.by = c("celltype"),   label = T, label_insitu = T ,    label.size = 5 ,
  reduction = "umap", theme_use = "theme_blank"   )
p4 <- p4 + scale_color_lancet()
ggsave(filename = "./result/CellDimPlot_Neutrophil.pdf", plot = p4, width = 7, height = 6)


# Plot 5: Neutrophil CellStatPlot
p5 <- CellStatPlot(srt = Neutrophil, stat.by = "celltype",  group.by = "group", label = F  , label.size = 5 )
p5 <- p5 + scale_fill_lancet()
ggsave(filename = "./result/CellStatPlot_Neutrophil.pdf", plot = p5, width = 6, height = 5)

### Il15ra umap
# Plot 6: Il15ra FeaturePlot
p6 <- FeatureDimPlot( Neutrophil , features = c( 'Il15ra' ))
ggsave(filename = "./result/FeatureDimPlot_Il15ra.pdf", plot = p6, width = 6, height = 5)


###############################################################################################
#### All cell Marker    

sce.markers <- readRDS(   './rds/sce_harmony_marker.rds'   )
library(dplyr)
top5 <- sce.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 

library(jjAnno)
library(scRNAtoolVis)

# Plot 7: All cell Marker DotPlot
p7 <- jjDotPlot(object = sce,   id = 'celltype',
                gene = unique(top5$gene),  dot.col = c('#C7E9C0'  , '#F7F7F7', '#8C2D04'   ) , 
                point.geom = F,  tile.geom = T ,   ytree =  F ,
                cluster.order =c( 
                  'NK', 'B cell','Proliferate T', 'cDC2', 'cDC1','Cd4T',
                  'Monocyte' , 'Cd8T', 'Macro', 'Neutrophil' ) , 
                legend.position =  'top' , textSize = 12, 
                rescale = T  ,   rescale.min = -2, rescale.max = 2,  midpoint = 0 )
p7 <- p7 + theme(text = element_text(size = 12 ))+ theme(axis.text = element_text(size = 12 ))  
ggsave(filename = "./result/DotPlot_AllCell_Markers.pdf", plot = p7, width = 10, height = 6)


### Neutrophil Cluster Marker 
Neutrophil.markers <- readRDS( './rds/Neutrophil_harmony_marker.rds'   )

library(dplyr)
top5 <- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  
top50 <- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)  

library(jjAnno)
library(scRNAtoolVis)

# Plot 8: Neutrophil Marker DotPlot
p8 <- jjDotPlot(object = Neutrophil,   id = 'celltype',
                gene = unique(top5$gene),  dot.col = c( '#9ECAE1' , '#F0F0F0' ,  '#CB181D' ) , 
                point.geom = F,  tile.geom = T ,   ytree =  F ,
                cluster.order =c( 'Neutrophil_C6' , 'Neutrophil_C5', 'Neutrophil_C4',
                                  'Neutrophil_C3' ,'Neutrophil_C2', 'Neutrophil_C1', 'Neutrophil_C0'   
                ) ,  legend.position =  'top' , textSize = 12, 
                rescale = T  ,   rescale.min = -2, rescale.max = 2,  midpoint = 0 )
p8 <- p8 + theme(text = element_text(size = 12 ))+ theme(axis.text = element_text(size = 12 ))
ggsave(filename = "./result/DotPlot_Neutrophil_Markers.pdf", plot = p8, width = 8, height = 6)


########################################################################################################################################################
### Differ gene  
library(clusterProfiler)
top70 <- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 70, wt = avg_log2FC)  

mygene <- top70[top70$cluster == 'Neutrophil_C4',]$gene

eego <- enrichGO(gene          = mygene,
                 keyType       =  'SYMBOL' ,
                 OrgDb         = org.Mm.eg.db,
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.1,
                 qvalueCutoff  = 0.1,
                 readable      = TRUE)

df <- data.frame(eego) %>%
  slice_head(n = 10) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,  function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio
df$Description <- factor(  df$Description,  levels = df$Description  )
df$ONTOLOGY  <- 'GO:CC'

## Plot 9: GO Enrichment
p9 <- ggplot(df) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = after_stat(index)),
                     n = 500,
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 6,shape = 21) +
  geom_line(aes(x = ratio*100,  y = Description, group = 1 ),
            orientation = "y",  linewidth = 1, color = "#FFCC00") +
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  scale_color_manual(values = c( "#729F4A")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) 
ggsave(filename = "./result/Enrichment_GO_CC.pdf", plot = p9, width = 8, height = 6)


############  KEGG  ###############################
gene.df <- bitr(mygene , fromType="SYMBOL",  toType=c("ENTREZID","ENSEMBL"),  ##基因名转换
                OrgDb = org.Mm.eg.db)

enrichKEGG <- enrichKEGG( unique( gene.df$ENTREZID  ) , organism='mmu', 
                          pvalueCutoff=0.1, pAdjustMethod='BH',qvalueCutoff=0.5,
                          minGSSize=1, maxGSSize=500, use_internal_data=F)

df <- data.frame(enrichKEGG) %>%
  slice_head(n = 15) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio
df$Description <- factor(df$Description,levels = df$Description)
df$ONTOLOGY <- 'KEGG'

# Plot 10: KEGG Enrichment
p10 <- ggplot(df) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = after_stat(index)),
                     n = 500,
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 6,shape = 21) +
  geom_line(aes(x = ratio*100,  y = Description, group = 1 ),
            orientation = "y",  linewidth = 1, color = "#FFCC00"  ) +
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  scale_color_manual(values = c( "#FDBF6F" )) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) 
ggsave(filename = "./result/Enrichment_KEGG.pdf", plot = p10, width = 8, height = 6)


############  ReactomePA  ###############################
library(ReactomePA)
library(clusterProfiler)
library(stringr)
library(openxlsx)
library(dplyr)


light_green <- "#729F4A"
light_red <- "#FDBF6F"
light_blue <- "#A6CEE3"
reactome_enrich_up = enrichPathway(gene = unique(gene.df$ENTREZID), #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
                                   organism='mouse', #物种选择 "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
                                   minGSSize = 10, #默认 10
                                   maxGSSize = 500, #默认 500
                                   pAdjustMethod = "BH", # 指定多重假设检验矫正的方法，选项包含 "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                   pvalueCutoff = 0.1,  #指定 p 值阈值（指定 1 以输出全部），默认为0.05
                                   qvalueCutoff = 0.2, #指定 q 值阈值（指定 1 以输出全部），默认0.2
                                   readable = T) # 是否将gene ID转换到 gene symbol


df <- data.frame( reactome_enrich_up  ) %>%
  slice_head(n = 15) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio
df$Description <- factor(df$Description,levels = df$Description)
df$ONTOLOGY <- 'Reactome'

# Plot 11: Reactome Enrichment
p11 <- ggplot(df) +
  ggforce::geom_link(aes(x = 0,  y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = after_stat(index)),
                     n = 500,
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 6,shape = 21) +
  geom_line(aes(x = ratio*100,  y = Description, group = 1 ),
            orientation = "y",  linewidth = 1, color = "#FFCC00") +
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  scale_color_manual(values = c( "#A6CEE3" )) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) 
ggsave(filename = "./result/Enrichment_Reactome.pdf", plot = p11, width = 8, height = 6)


############################################################################################################

### Cellchat 

cellchat_NC <- readRDS( './rds/cellchat_NC.rds'  )
cellchat_NC <- netAnalysis_computeCentrality(cellchat_NC, slot.name = "netP")  

cellchat_sg <- readRDS( './rds/cellchat_sg.rds'  )
cellchat_sg <- netAnalysis_computeCentrality(cellchat_sg, slot.name = "netP") 

levels(cellchat_NC@idents)  

object.list <- list(   NC = cellchat_NC , sgTm9sf2 = cellchat_sg )
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#####--------------------------Plot------------------------------------------------------------------------------------------------------------
# Plot 12: CellChat Circle Plot (需要用 pdf() 手动保存，因为是基础绘图函数生成)
pdf("./result/CellChat_Interaction_Circle.pdf", width = 12, height = 6)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list) ) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge=T,  
                   targets.use = c("CD4T_C1_Treg" , "CD4T_C5_Treg" , "CD8T_C0_Memory", "CD8T_C2_Effect",  "CD8T_C3_Naive" ,   
                                   "CD8T_C4_Exhausted", "CD8T_C6_Exhausted" ,"CD8T_C7_ISG" ),  sources.use =  c(   "Neutrophil_C4" ) , 
                   edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

########---------------------------------paired Signal----------------------------------##################################################################################################################################
# Plot 13: CellChat Bubble Plot
p13 <- netVisual_bubble(cellchat, sources.use = c(9 ), 
                        targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)
ggsave(filename = "./result/CellChat_Interaction_Bubble.pdf", plot = p13, width = 10, height = 8)