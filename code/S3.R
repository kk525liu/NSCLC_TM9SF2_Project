

##########-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Neutrophil <- readRDS( './rds/Neutrophil_harmony.rds'  )
pancreas_sub <- Neutrophil
pancreas_sub$celltype

pancreas_sub <- RunDEtest(srt = Neutrophil, group_by = "celltype", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = pancreas_sub, group_by = "celltype")


DEGs <- pancreas_sub@tools$DEtest_celltype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 0.5 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "CSPA"))
ht <- FeatureHeatmap(
  srt = pancreas_sub, group.by = "celltype", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Mus_musculus", db = c("GO_BP", "KEGG", "WikiPathway"), anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4
)
print(ht$plot)





##----------------------------------------------------------------------------------------------------------

pancreas_sub <- RunGSEA(
  srt = pancreas_sub, group_by = "celltype", db = c("GO_MF" , 'GO_BP', 'GO_CC' ,  'KEGG'  ,'Reactome' )
  , species = "Mus_musculus",  
  DE_threshold = "p_val_adj < 0.05"
)

##        
result <- GSEAPlot(
  srt = pancreas_sub, group_by = "celltype", group_use = c('Neutrophil_C4'),db =c( 'KEGG'  ),
  plot_type = "bar" , padjustCutoff = 0.05,  topTerm = 100 , 
)
data <- result$data
GSEAPlot(srt = pancreas_sub, group_by = "celltype", db = 'KEGG' ,
         group_use = "Neutrophil_C4", id_use = "mmu04612")






