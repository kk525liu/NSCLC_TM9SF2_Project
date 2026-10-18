


## install.packages("qs", repos = "https://packagemanager.posit.co/cran/__linux__/focal/2025-12-01")
library(qs)


sce <- qread(  '/home/kk/Lung_cancer_Data/scRNA-seq/data/Li_lab/LYX/sce1.2_harmony_renamed_newLevel.qs'      )

colnames(sce@meta.data)  

FeatureDimPlot( sce , features = c( 'Col6a1' )    )


CellDimPlot(    
  srt = sce, group.by = c("celltype"),label = T, label_insitu = T , label.size = 3 , 
  reduction = "umap", theme_use = "theme_blank"    )


