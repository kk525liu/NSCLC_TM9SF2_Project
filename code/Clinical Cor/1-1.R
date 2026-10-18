
##########  导入表达矩阵  

tcga_LUAD <- qread('./RDS/bulk_RNA_NSCLC/TCGA_LUAD_expr_allmerge.qs')
tcga_LUSC <- qread('./RDS/bulk_RNA_NSCLC/TCGA_LUSC_expr_allmerge.qs')

Same_gene <- intersect( colnames( tcga_LUAD) , colnames(tcga_LUSC))
tcga_LUAD <- tcga_LUAD[,Same_gene]
tcga_LUSC <- tcga_LUSC[,Same_gene]
table( colnames( tcga_LUAD) == colnames( tcga_LUSC)  )


tcga_data <- rbind(tcga_LUAD ,tcga_LUSC  )
tcga_data  <-  as.data.frame( t(tcga_data)  )
rownames(tcga_data)
expr <- as.matrix(tcga_data )
head( expr ,c(4,4) )

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
gs.exp <- GSVA::gsva(expr, gene.set, kcdf = "Gaussian", method = "ssgsea",min.sz = 5)
row.names(gs.exp)<-'Tm9sf2_geneset'
gs.exp <- as.data.frame(t( gs.exp ))



#######################--------------------------------------------------------------------------------------------------------

sample<-str_sub(  rownames(gs.exp)  )
sample<-str_replace_all(sample,"\\.","-")
rownames(gs.exp)<-sample


#####设置表型信息
samplename <- as.data.frame( rownames(gs.exp))
samplename$sample <- str_split(samplename$`rownames(gs.exp)` , '-' , simplify = T)[,4] 
table(samplename$sample)

gs.exp$Sample <-samplename$sample 

for (i in c('01A' , '01B' , '01C' , '02A' )) {
  
  gs.exp[which(gs.exp$Sample==i),]$Sample= "Tumor"
}
for (i in c('11A' ,'11B')) {
  
  gs.exp[which(gs.exp$Sample==i),]$Sample= "Normal"
}


table(gs.exp$Sample)

gs.exp$Sample <-  factor(gs.exp$Sample , levels = c('Normal'  , 'Tumor'  ))


library(ggpubr)
library(ggsci)
## 指定自己想要的比较
my_comparisons <- list(c( 'Normal'  , 'Tumor' ) )
p1<-ggboxplot(gs.exp, x = "Sample", y = "Tm9sf2_geneset", width = 0.4 , title = 'TCGA NSCLC',
              color = "Sample",add = "jitter", palette = "jama")+ 
  stat_compare_means(comparisons = my_comparisons)+scale_color_lancet()

p1

#############-------------------------------------------------------------------------------------------------------

meta_LUAD <- read.csv('./data/Bulk-data/TCGA-LUAD.survival (1).tsv' , sep = '') 
meta_LUSC <- read.csv('./data/Bulk-data/TCGA-LUSC.survival.tsv' , sep = '') 
meta_LUSC <- dplyr::select( meta_LUSC , colnames(meta_LUAD ))
meta <- na.omit( rbind(meta_LUAD , meta_LUSC )   )



meta2 <- meta[match(  rownames(gs.exp)  , meta$sample ),  ]
table( meta2$sample == rownames( gs.exp ) )

sur_merge <- cbind(meta2 , gs.exp )

library(survival)
library(survminer)
exprSet=sur_merge  #套流程
sur_merge$gs.exp_level = ifelse(sur_merge[,'Tm9sf2_geneset']>median(sur_merge[,'Tm9sf2_geneset']),'high','low')
sfit1=survfit(Surv(OS.time, OS)~gs.exp_level, data=sur_merge)
ggsurvplot(sfit1,pval =TRUE, data = sur_merge,palette = c('red' ,'black') 
           ,
           risk.table = TRUE )


# 将表达水平分为四个等级（四分位法）
quantiles <- quantile(sur_merge$Tm9sf2_geneset, probs = c(0, 0.25, 0.5, 0.75, 1))
sur_merge$gs.exp_level <- cut(
  sur_merge$Tm9sf2_geneset,
  breaks = quantiles,
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)

# 确保分组变量是正确的因子类型
sur_merge$gs.exp_level <- factor(sur_merge$gs.exp_level, levels = c("Q1", "Q2", "Q3", "Q4"))

# 创建生存模型
sfit1 <- survfit(Surv(OS.time, OS) ~ gs.exp_level, data = sur_merge)

# 绘制生存曲线（使用四种不同颜色）
ggsurvplot(
  sfit1,
  pval = TRUE,                    # 显示P值
  data = sur_merge,               # 数据集
  palette = c("blue", "green", "orange", "red"),  # 四分位数颜色（Q1到Q4）
  risk.table = TRUE,              # 显示风险表
  legend.title = "Expression",     # 图例标题
  legend.labs = c("Q1 (low)", "Q2", "Q3", "Q4 (high)")  # 修改图例标签
)

#####################################################################################################
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

sample<-str_sub(rownames(gs.exp))
sample<-str_replace_all(sample,"\\.","-")
rownames(gs.exp)<-sample
gs.exp$Sample <- str_split(rownames( gs.exp), 
                           pattern =  '-' , simplify =  T )[,4]
#####剔除正常样本
unique(gs.exp$Sample)

for (i in c('01A' , '01B' ,"02A" , "01C" )) {
  
  gs.exp[which(gs.exp$Sample==i),]$Sample= "Tumor"
}

gs.exp_Tumor  <- gs.exp[which(gs.exp$Sample == 'Tumor'  ),]
gs.exp_Tumor <- dplyr::select( gs.exp , -Sample   )  


for (i in 1: length(  colnames( gs.exp_Tumor )  )) {
  gs.exp_Tumor[,i] <-  as.numeric( gs.exp_Tumor[,i]  )
}


###########----------------------------------------------------------------------------------------------------------------------

library("ggsci")
library("ggpubr")

colnames(gs.exp_Tumor)
ggscatter(gs.exp_Tumor, x = "TM9SF2_GENESET", y = "neutrophil_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", title = 'TCGA_NSCLC',
          add = "reg.line", conf.int = TRUE, add.params = list(color = "blue"),
          cor.coef = TRUE, cor.method = "pearson",
          color = "grey",
)




