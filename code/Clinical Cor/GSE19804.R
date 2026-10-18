



#### bulk矩阵TPM 行为基因，列为样本  
load('./RDS/Li lab/NK_Correlation/GSE19804_expr.rds'  )
Sample <- 'GSE19804'

exp_merge <- as.matrix( result2)
head(exp_merge , c( 5 , 5 ))
expr <- exp_merge[-1,]

dim(expr)






### 相关性
##导入TM9差异基因集
DEG_limma_voom <- qread('./RDS/Li lab/bulk_RNA_Seq/LYX_cellline.qs')
colnames(DEG_limma_voom)
DEG_limma_voom <- DEG_limma_voom[which(DEG_limma_voom$P.Value < 0.05),]

order_DiffGene <- DEG_limma_voom[order(DEG_limma_voom$logFC ,  decreasing = F),]
DownGene <- rownames( order_DiffGene[1:20,] )

library(homologene)
genelist <- DownGene
gene_homo <- homologene(genelist,   inTax = 10090,   outTax = 9606  )
DownGene_homo <- gene_homo$`9606`

#eg2 <- bitr(DownGene_homo, fromType = "SYMBOL", toType =  "ENSEMBL", OrgDb = "org.Hs.eg.db")
#gene<-eg2$ENSEMBL
gene.set<-list() 
gene.set[[1]]<-DownGene_homo   

##hallmark转移基因集GSVAscore
# gsva 分析，read counts 使用泊松分布，通路至少包含 10 个基因
library(TCGAbiolinks)
library(GSVA)
library(GSEABase)
gs.exp1 <- GSVA::gsva(expr,  gene.set, kcdf = "Gaussian", method = "ssgsea",min.sz = 5)
row.names(gs.exp1)<-'TM9SF2_GENESET'



#####抗原呈递基因集
AP <- read.gmt( './data/Li_lab/WL/Srsf3/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION.v2024.1.Hs.gmt')
AP_gene <- AP$gene
AP_gene
Neu <- c('CSF3R', 'FPR1' , 'CXCR2', 'G0S2',  'NAMPT','MNDA' ,'IL1R2',  'FPR2')
AP_gene_Neu <- c( AP_gene , Neu  )

gene.set<-list()
gene.set[[1]]<-AP_gene_Neu
# gsva 分析，read counts 使用泊松分布，通路至少包含 10 个基因
gs.exp2 <- GSVA::gsva(expr, gene.set, kcdf = "Gaussian", method = "ssgsea",min.sz = 5  )
row.names(gs.exp2)<-'neutrophil_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION'  

### MHC1                     
MHC1<- read.gmt( './data/Li_lab/WL/Srsf3/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I.v2024.1.Hs.gmt')
gene<-MHC1$gene
gene.set<-list()
gene.set[[1]]<-gene
gs.exp3 <- GSVA::gsva(expr, gene.set, kcdf = "Gaussian", method = "ssgsea",min.sz = 5  )
row.names(gs.exp3)<-'neutrophil_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I'  

### MHC2                     
MHC2<- read.gmt( './data/Li_lab/WL/Srsf3/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II.v2024.1.Hs.gmt')
gene<-MHC2$gene
gene.set<-list()
gene.set[[1]]<-gene
gs.exp4 <- GSVA::gsva(expr, gene.set, kcdf = "Gaussian", method = "ssgsea",min.sz = 5  )
row.names(gs.exp4)<-'neutrophil_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II'  



table( colnames( gs.exp1  )  == colnames( gs.exp2  )   )
table( colnames( gs.exp1  )  == colnames( gs.exp3  )   )
table( colnames( gs.exp1  )  == colnames( gs.exp4  )   )

gs.exp <- as.data.frame(  t( rbind( gs.exp1, gs.exp2 , gs.exp3 ,  gs.exp4  )  ))

###########----------------------------------------------------------------------------------------------------------------------

library("ggsci")
library("ggpubr")

colnames(gs.exp)
ggscatter(gs.exp, x = "TM9SF2_GENESET", y = "neutrophil_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II", title = 'GSE19804',
          add = "reg.line", conf.int = TRUE, add.params = list(color = "blue"),
          cor.coef = TRUE, cor.method = "pearson",
          color = "grey",
)

