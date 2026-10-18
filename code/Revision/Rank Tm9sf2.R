


# 加载必要的包
library(readxl)
library(tidyverse)
library(ggrepel)

# 1. 读取数据
file_path <- "/home/kk/Lung_cancer_Data/scRNA-seq/data/Li_lab/LYX/B6_PBS_Rag1_median_rra.gene_summary.xlsx"
df <- read_excel(file_path)

# 2. 数据预处理与排名计算
plot_data <- df %>%
  # 剔除 pos|score 为空的无效行
  filter(!is.na(`pos|score`)) %>%
  # 按照正向选择得分从小到大(越小越显著)进行严格排序
  arrange(`pos|score`) %>%
  # 生成 Rank 排名列 (由于存在得分相同的基因，row_number() 会赋予它们连续的绝对排名)
  mutate(
    Rank = row_number(),
    # 计算 -log10(score) 用于 Y 轴，使图形呈现瀑布流曲线
    log_score = -log10(`pos|score`),
    
    # 设定分组变量，用于后续自定义颜色
    # 将 Tm9sf2 设为 Target，前 50 名设为 Top Hits（可自行修改），其余为 Background
    Gene_Group = case_when(
      id == "Tm9sf2" ~ "Target",
      Rank <= 50 ~ "Top Hits", 
      TRUE ~ "Background"
    )
  )

# 提取目标基因所在行，方便后续专门为其画点和添加标签
target_gene <- plot_data %>% filter(id == "Tm9sf2")

# 3. 绘制 SCI 风格图表
p <- ggplot(plot_data, aes(x = Rank, y = log_score)) +
  
  # 第1层：绘制背景基因 (浅灰色，透明度高，不抢眼)
  geom_point(data = filter(plot_data, Gene_Group == "Background"), 
             color = "#E0E0E0", size = 1.2, alpha = 0.6) +
  
  # 第2层：绘制 Top 显著基因 (例如用经典的 Nature 蓝)
  geom_point(data = filter(plot_data, Gene_Group == "Top Hits"), 
             color = "#4DBBD5", size = 1.8, alpha = 0.9) +
  
  # 第3层：高亮目标基因 Tm9sf2 (使用 Nature 红，尺寸放大)
  geom_point(data = target_gene, 
             color = "#E64B35", size = 3.5) +
  
  # 第4层：添加防重叠的基因名称标签 (ggrepel)
  geom_text_repel(data = target_gene,
                  # 标签内容：基因名 (斜体) + Rank 排名
                  aes(label = sprintf("italic('%s')~'(Rank: %d)'", id, Rank)),
                  parse = TRUE, # 允许解析表达式以实现部分斜体
                  color = "#E64B35",
                  size = 5,
                  box.padding = 1.5,
                  point.padding = 0.5,
                  nudge_x = 200,   # 将标签向右推开，避免挡住曲线
                  nudge_y = 0.3,   # 将标签向上推开
                  segment.color = "black", # 连接线颜色
                  segment.size = 0.5) +
  
  # 第5层：设置坐标轴标签
  labs(
    x = "Gene Rank (Positive Selection)",
    y = expression(bold(-log[10]("RRA Score"))) # Y轴使用标准的数学公式表达
  ) +
  
  # 第6层：SCI 经典主题美化 (黑白分明，加粗坐标轴，去除冗余网格)
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 0.8, color = "black"), # 粗坐标轴线
    axis.ticks = element_line(linewidth = 0.8, color = "black"), # 粗刻度线
    axis.text = element_text(color = "black", size = 12),       # 黑色刻度文字
    axis.title = element_text(face = "bold", size = 14),        # 加粗轴标题
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)        # 留白边缘
  )

# 打印图像
print(p)
