

file.edit("~/.Rprofile")

rm(list = ls ())
gc()
set.seed(525)
library(SCP)
library(Seurat) 
#install.packages('ggsci')
library(ggsci) 
library(qs)
library(DelayedMatrixStats)


load( "/home/kk/Lung_cancer_Data/scRNA-seq/RDS/Li lab/All/Neutrophil_harmony.rds")
#Naive marker
FeatureDimPlot(Neutrophil , 
               features = c( 'Camp' , 'Ltf' , 'Lcn2' , 'Retnlg' , 'Mmp9' ))

###Cytotrace
matrix <- as.matrix(Neutrophil@assays$RNA@counts)
phenot <- Neutrophil$celltype
phenot <- as.character(phenot)
names(phenot) <- rownames(Neutrophil@meta.data)
emb <- Neutrophil@reductions[["umap"]]@cell.embeddings
library(CytoTRACE) #加载包

results <- CytoTRACE(matrix, ncores = 8)
plotCytoGenes(results, numOfGenes = 10)

###可视化
plotCytoTRACE(results, phenotype = phenot)
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
pancreas_sub2 <- pancreas_sub

p <- CellDimPlot(  
  srt = Neutrophil,  group.by = c("celltype"),   label = T, label_insitu = T ,    label.size = 5 ,
  reduction = "umap", theme_use = "theme_blank"   )
p + scale_color_lancet()   

###Slingshot

pancreas_sub <- Neutrophil
unique( pancreas_sub$celltype )
pancreas_sub <- RunSlingshot(srt = pancreas_sub, group.by = "celltype",show_plot = F , 
                             start = 'Neutrophil_C6' , end = 'Neutrophil_C4'   , 
                             reduction = "UMAP")

FeatureDimPlot(pancreas_sub, features = paste0("Lineage", 1:3),
               reduction = "UMAP", theme_use = "theme_blank") 

CellDimPlot(pancreas_sub, group.by = "celltype", reduction = "UMAP", lineages_linewidth =  2 , 
            palcolor = c("#00468BFF", "#ED0000FF", "#42B540FF" ,"#0099B4FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF"),
            lineages = paste0("Lineage", 1:3), lineages_span = 0.1)

CellDimPlot(pancreas_sub, group.by = "celltype", reduction = "UMAP",  split.by = 'group' , 
            lineages_linewidth =  2 , 
            palcolor = c("#00468BFF", "#ED0000FF", "#42B540FF" ,"#0099B4FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF"),
            lineages = paste0("Lineage", 1:3), lineages_span = 0.1)


rlang::last_trace()


##Dynamic features
pancreas_sub <- RunDynamicFeatures(srt = pancreas_sub, lineages = c("Lineage3"), n_candidates = 200)
ht <- DynamicHeatmap(
  srt = pancreas_sub, lineages = c("Lineage3" ),
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage3",
  species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "viridis", cell_annotation = "celltype",
  separate_annotation = list("celltype", c("Ifit1")), separate_annotation_palette = c("Paired", "Set1"),
 feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(ht$plot)

DynamicPlot(
  srt = pancreas_sub, lineages = c( "Lineage2"), group.by = "celltype",
  features = c("Cxcl10"),point_palcolor =  c("#00468BFF", "#ED0000FF", "#42B540FF" ,"#0099B4FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF"),
  compare_lineages = TRUE, compare_features = FALSE
)

