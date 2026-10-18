


setwd('/home/kk/Lung_cancer_Data/scRNA-seq/')  

rm(list = ls())
gc()

library(Seurat)
library(qs)
library(SCP)
library(ggsci)

############################################################
## A , B
sce <- qread( "./RDS/Li lab//All/sce_harmony.qs" )

FeatureDimPlot( sce , features = c( 'Ptprc' ))

CellDimPlot(    
  srt = sce, group.by = c("celltype"),label = T, label_insitu = T , label.size = 3 , 
  reduction = "umap", theme_use = "theme_blank"  )

CellStatPlot(srt = sce, stat.by = "celltype",  group.by = "group", label = F  , label.size = 5 )

#############---------------------------------------------------------------------------------###################
#### C, D
load( "./RDS/Li lab/All/Neutrophil_harmony.rds")    
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



p <- CellDimPlot(  
  srt = Neutrophil,  group.by = c("celltype"),   label = T, label_insitu = T ,    label.size = 5 ,
  reduction = "umap", theme_use = "theme_blank"   )
p + scale_color_lancet()
# 获取 Lancet 调色板的颜色值
pal_lancet("lanonc")(7)  # 默认返回 9 种颜色



p <-CellStatPlot(srt = Neutrophil, stat.by = "celltype",  group.by = "group", label = F  , label.size = 5 )
p + scale_fill_lancet()
save( Neutrophil , file = './RDS/Li lab/All/Neutrophil_harmony.rds' )     
load(  './RDS/Li lab/All/Neutrophil_harmony.rds' )

FeatureDimPlot( Neutrophil , features = c( 'Vegfa', 'Thbs1' ))
FeatureDimPlot( Neutrophil , features = c( 'Bmpr1a' ))




# NEJM 调色板
colors_nejm <- pal_nejm("default")(7)
# NPG 调色板 (Nature)
colors_npg <- pal_npg("nrc")(7)
# AAAS 调色板 (Science)
colors_aaas <- pal_aaas("default")(7)
# Lancet 调色板
colors_lancet <- pal_lancet("lanonc")(7)  # 或 "lanonc"-> Lancet Oncology 风格
# JAMA 调色板
colors_jama <- pal_jama("default")(7)
# JCO 调色板 (Journal of Clinical Oncology)
colors_jco <- pal_jco("default")(7)
# D3 调色板 (ggsci 自带)
colors_d3 <- pal_d3("category20")(7)  # 可选 "category10", "category20b" 等
# LocusZoom 调色板
colors_locuszoom <- pal_locuszoom("default")(7)
# IGV 调色板 (Integrative Genomics Viewer)
colors_igv <- pal_igv("default")(7)  # 可选 "alternating"
# Chicago 调色板
colors_uchicago <- pal_uchicago("default")(7)  # 可选 "light", "dark"
# Tron Legacy 调色板
colors_tron <- pal_tron("legacy")(7)  # 可选 "basic"
# Rick and Morty 调色板
colors_rick <- pal_rickandmorty("schwifty")(7)
# Simpsons 调色板
colors_simpsons <- pal_simpsons("springfield")(7)


###############################################################################################
library(ggsci)
p + scale_color_nejm()
p + scale_color_npg()
p + scale_color_aaas()
p + scale_color_lancet()
p + scale_color_jama()
p + scale_color_jco()
p + scale_color_d3()
p + scale_color_locuszoom()
p + scale_color_igv()
p + scale_color_uchicago()
p + scale_color_tron()
p + scale_color_rickandmorty()
p + scale_color_simpsons()

##
p + scale_fill_nejm()
p + scale_fill_npg()
p + scale_fill_aaas()
p + scale_fill_lancet()
p + scale_fill_jama()
p + scale_fill_jco()
p + scale_fill_d3()
p + scale_fill_locuszoom()
p + scale_fill_igv()
p + scale_fill_uchicago()
p + scale_fill_tron()
p + scale_fill_rickandmorty()
p + scale_fill_simpsons()
p + scale_fill_material()
p + scale_fill_gsea()



###############################################################################################
### E 中性粒亚群Marker
load(file = "RDS/Li lab//Neutrophil_harmony_marker.rds")
##write.csv(Neutrophil.markers ,file = './result/Neutrophil_harmony_marker.csv')
library(dplyr)
top5 <- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  
top20 <- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)  
top100 <- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(top100 ,file = './result/Neutrophil_harmony_marker_top100.csv')
library(jjAnno)
library(scRNAtoolVis)
unique( Neutrophil$celltype )
p <-jjDotPlot(object = Neutrophil,   id = 'celltype',
              gene = unique(top5$gene),  dot.col = c( '#C7E9C0'  , '#F7F7F7', '#8C2D04'  ) , 
              point.geom = F,  tile.geom = T ,   ytree =  F ,
              cluster.order =c( 'Neutrophil_C6' , 'Neutrophil_C5', 'Neutrophil_C4',
                                'Neutrophil_C3' ,'Neutrophil_C2', 'Neutrophil_C1', 'Neutrophil_C0'   
              ) ,  legend.position =  'top' , textSize = 12, 
              rescale = T  ,   rescale.min = -2, rescale.max = 2,  midpoint = 0 )
p + theme(text = element_text(size = 12 ))+ theme(axis.text = element_text(size = 12 ))

#### 所有细胞亚群的Marker
sce <- qread( "./RDS/Li lab//All/sce_harmony.qs")
load( "RDS/Li lab/All/sce_harmony_marker.rds")
library(dplyr)
top5 <- sce.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 
top2 <- sce.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 

library(jjAnno)
library(scRNAtoolVis)
unique( sce$celltype )
p <-jjDotPlot(object = sce,   id = 'celltype',
              gene = unique(top5$gene),  dot.col = c('#C7E9C0'  , '#F7F7F7', '#8C2D04'   ) , 
              point.geom = F,  tile.geom = T ,   ytree =  F ,
              cluster.order =c( 
                                 'NK', 'B cell','Proliferate T', 'cDC2', 'cDC1','Cd4T',
                                'Monocyte' , 'Cd8T', 'Macro', 'Neutrophil' ) , 
              legend.position =  'top' , textSize = 12, 
              rescale = T  ,   rescale.min = -2, rescale.max = 2,  midpoint = 0 )
p + theme(text = element_text(size = 12 ))+ theme(axis.text = element_text(size = 12 ))  

gene <- c( "S100a9"  ,   "S100a8"  ,   "Ccl24" ,     "Arg1" ,       "Ccl5" ,      "Gzma"   ,    "Lyz2"   ,    "Plac8",     
 "Inpp4b"   ,  "St6galnac3", "Ccr7"   ,    "Tmcc3"    ,  "Immp2l"  ,   "Cst3"     ,  "Stmn1" ,     "Pclaf"    ,  
 "Ebf1"    ,   "Igkc"  ,     "Clnk" ,     "Prf1"  )  

FeatureDimPlot( sce ,features =  gene  , ncol = 5 )
### '#9ECAE1' , '#F0F0F0' ,  '#CB181D'  淡色系蓝白红
##  "#ffffe0", "#00ced1", "#006400"       黄蓝青  
##   '#DEEBF7' , '#9ECAE1' , '#6A1B9A'    蓝白紫 
##  '#D9D9D9' ,  '#FFFFFF' ,  '#D94701'  浅灰 → 白 → 深橙红
###  "#DCDCDC", "#FFFFFF", "red"          灰白红  
##  '#C7E9C0'  , '#F7F7F7', '#8C2D04'   淡青 → 米白 → 深棕红
##  '#A6BDDB' ,  '#FFFFCC' ,  '#54278F'  浅蓝 → 浅黄 → 深紫

###########################################################################################################################################################################

####F
load( "./RDS/Li lab/All/Neutrophil_harmony.rds")
FeatureDimPlot( Neutrophil , features = 'Il15ra')
VlnPlot(  Neutrophil , features = 'Vegfa' )

### G 
###Slingshot





########################################################################################################################################################
### IJK  
load(file = "RDS/Li lab//Neutrophil_harmony_marker.rds")
g

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
## 横坐标为  ， 纵坐标是富集到的信号通路
## 
ggplot(df) +
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

######
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


ggplot(df) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = after_stat(index)),
                     n = 500,
                     # color = "#FF0033",
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


light_green <- "#729F4A"
light_red <- "#FDBF6F"
light_blue <- "#A6CEE3"


library(ReactomePA)
library(clusterProfiler)
library(stringr)
library(openxlsx)
library(dplyr)

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

## 横坐标为  ， 纵坐标是富集到的信号通路   
## 
ggplot(df) +
  ggforce::geom_link(aes(x = 0,  y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = after_stat(index)),
                     n = 500,
                     # color = "#FF0033",
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y                                                                                                    = Description),
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




############################################################################################################
library(qs)
cellchat <- qread('./RDS/Li lab//All/cellchat_Neu_T.qs')
library(patchwork)
groupSize <- as.numeric(table(cellchat@idents))


levels(cellchat@idents)
par(mfrow = c(1,2), xpd=TRUE)  
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, vertex.label.cex = 0.8, edge.label.cex =0.8,
                 label.edge= T, title.name = "Number of interactions", 
                 sources.use = c('Neutrophil_C4_IFIT1'),targets.use =c("CD4T_C1_Treg"  ,   "CD4T_C5_Treg" , "CD8T_C0_Memory" , "CD8T_C2_Effect"  ,
                                                                       "CD8T_C3_Naive" , "CD8T_C4_Exhausted"  , "CD8T_C6_Exhausted" ,  "CD8T_C7_ISG"    ))

weight <-  cellchat@net$weight


netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, vertex.label.cex = 0.8,
                 label.edge= T, edge.label.cex =0.8,  title.name = "Number of interactions", 
                 sources.use = c("Neutrophil_C1_Rps8" , "Neutrophil_C2", "Neutrophil_C3_Il1b" ,        
                                 "Neutrophil_C4_IFIT1" , "Neutrophil_C5", "Neutrophil_C6" ),
                 targets.use  = 'CD8T_C2_Effect')
















