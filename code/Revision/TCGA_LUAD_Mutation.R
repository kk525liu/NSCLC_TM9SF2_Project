




library(TCGAbiolinks)
library(maftools)


# 构造查询，获取最新的非同义突变和同义突变整合数据
query_mut <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

# 下载并读取数据（这会在本地生成文件，注意预留几百MB空间）
GDCdownload(query_mut)
maf_prep <- GDCprepare(query_mut)

# 使用 maftools 读取 MAF 对象
# maftools 默认会自动过滤掉同义突变 (Synonymous mutations)，只保留能引起氨基酸改变的突变
luad_maf <- read.maf(maf = maf_prep)


#########

# 查看 LUAD 队列总突变图谱摘要（可以用来输出 SCI 级别的 Oncoplot）
# oncoplot(maf = luad_maf, top = 10) # 前10个高频突变基因，KRAS通常排前三

# 获取所有基因的突变统计信息
gene_summary <- getGeneSummary(luad_maf)

# 筛选 KRAS 基因
kras_stat <- gene_summary[Hugo_Symbol == "KRAS"]
print(kras_stat)

# 计算具体的突变百分比
total_samples <- as.numeric(luad_maf@summary[3, 'summary']) # 获取有突变测序数据的总样本数，通常是 500 多例
kras_mut_samples <- kras_stat$MutatedSamples

kras_mut_rate <- (kras_mut_samples / total_samples) * 100
cat(sprintf("TCGA-LUAD 队列中 KRAS 的非同义突变率为: %.2f%%\n", kras_mut_rate))

# 提取具体的突变氨基酸变化 (Amino acid change)
kras_mut_details <- subsetMaf(maf = luad_maf, genes = "KRAS")
# 查看最常见的突变类型分布
table(kras_mut_details@data$HGVSp_Short) 
# 您将看到类似 p.G12C, p.G12V, p.G12D 的频次统计


#############   Plot   
library(ggplot2)
library(dplyr)
library(forcats)

# 2. 构造数据框 
kras_mut_data <- data.frame(
  Mutation = c("p.G12C", "p.G12V", "p.G12D", "p.G12A", "p.G12S", 
               "p.G13C", "p.L19F", "p.G13D", "p.Q61L", "p.A146P", 
               "p.D33E", "p.K88*", "p.Q61H"),
  Count = c(59, 46, 18, 17, 6, 6, 3, 2, 2, 1, 1, 1, 1)
)

# 3. 数据处理：计算百分比并根据突变频率排序
kras_mut_data <- kras_mut_data %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    # 区分 Codon 12 热点与其他突变，用于配色
    Highlight = ifelse(grepl("G12", Mutation), "Codon 12 Hotspots", "Other Mutations") 
  ) %>%
  # 按照 Count 从大到小排序因子水平
  mutate(Mutation = fct_reorder(Mutation, Count))

# 4. 绘制符合 SCI 审美的棒棒糖图 (Lollipop Plot)
ggplot(kras_mut_data, aes(x = Mutation, y = Count, color = Highlight)) +
  # 绘制线段和点
  geom_segment(aes(x = Mutation, xend = Mutation, y = 0, yend = Count), linewidth = 1.2) +
  geom_point(size = 5) +
  # 添加数值标签 (展示 Count 和 %)
  geom_text(aes(label = sprintf("%d (%.1f%%)", Count, Percentage)), 
            hjust = -0.3, size = 3.5, color = "black", fontface = "bold") +
  # 翻转坐标轴，更利于长文本展示
  coord_flip() +
  # 自定义高级学术配色
  scale_color_manual(values = c("Codon 12 Hotspots" = "#E64B35", "Other Mutations" = "#4DBBD5")) +
  # 调整 y 轴范围，为标签留出空间
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # 完善标签与标题
  labs(
    title = "Landscape of KRAS Mutations in TCGA-LUAD Cohort",
    subtitle = paste0("Total Non-synonymous Mutations: n = ", sum(kras_mut_data$Count)),
    x = "Amino Acid Change",
    y = "Frequency (Count)",
    color = "Mutation Region"
  ) +
  # 使用经典主题，去除多余网格线
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

# 提示：您可以使用 ggsave("KRAS_Mutation_Landscape.pdf", width = 8, height = 6) 导出为矢量图


library(ggplot2)
library(dplyr)


wt_samples <- total_samples - kras_mut_samples

# 构造绘图数据，将低频突变合并为 "Other" 以保持图形整洁
donut_data <- data.frame(
  Category = c("KRAS Wild-Type", "KRAS G12C", "KRAS G12V", "KRAS G12D", "Other KRAS Mut"),
  Count = c(wt_samples, 59, 46, 18, 40) # 40 是除了C/V/D之外其余突变的总和
)

# 计算比例和排版位置
donut_data <- donut_data %>%
  arrange(desc(Category)) %>%
  mutate(
    Fraction = Count / sum(Count),
    Ymax = cumsum(Fraction),
    Ymin = c(0, head(Ymax, n = -1)),
    Label_Pos = (Ymax + Ymin) / 2,
    # 格式化标签：名称 + 占比
    Label = paste0(Category, "\n", round(Fraction * 100, 1), "%")
  )



# 绘图
ggplot(donut_data, aes(ymax = Ymax, ymin = Ymin, xmax = 4, xmin = 3, fill = Category)) +
  geom_rect(color = "white", linewidth = 1) +
  # 添加文本标签，远离圆心
  geom_text(aes(x = 4.8, y = Label_Pos, label = Label), size = 4.5, fontface = "bold") +
  # 转换为极坐标，形成圆环
  coord_polar(theta = "y") +
  # 设定 x 轴范围，控制圆环中空的粗细 (xmin越小，环越粗)
  xlim(c(1, 5)) +
  # 高级学术配色：野生型用灰色，热点突变用同色系的渐变或对比色
  scale_fill_manual(values = c(
    "KRAS Wild-Type" = "#E0E0E0", 
    "KRAS G12C" = "#E64B35", 
    "KRAS G12V" = "#F39B7F", 
    "KRAS G12D" = "#4DBBD5", 
    "Other KRAS Mut" = "#91D1C2"
  )) +
  theme_void() + # 去除所有背景元素
  theme(
    legend.position = "none", # 标签已在图上，可省略图例
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = -20))
  ) +
  labs(title = "Proportion of KRAS Mutations in TCGA-LUAD")






