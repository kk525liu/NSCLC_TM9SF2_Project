


# 加载R包
library(readxl)
library(tidyverse)

# 1. 读取数据
file_path <- "/home/kk/Lung_cancer_Data/scRNA-seq/data/Li_lab/LYX/T cell Flow data.xlsx"

# 小技巧：设置 col_names = FALSE，把表头也当做数据读进来，方便后续处理合并单元格
raw_data <- read_excel(file_path, col_names = FALSE)

# 2. 提取并清洗多级表头
# 提取第一行 (Marker)，并使用简单的循环填补合并单元格留下的 NA 空白
markers <- as.character(raw_data[1, -1])
for (i in seq_along(markers)) {
  if (is.na(markers[i]) || markers[i] == "NA") {
    markers[i] <- markers[i - 1]
  }
}

# 提取第二行 (实验组 Group)
groups <- as.character(raw_data[2, -1])

# 拼接形成标准化的列名，例如 "Tcm_sgNC"
new_col_names <- c("Cell_Type", paste(markers, groups, sep = "_"))

# 3. 提取核心数据并格式化
df_data <- raw_data[-c(1, 2), ] # 剔除前两行表头
colnames(df_data) <- new_col_names

df_clean <- df_data %>%
  # 关键步骤：向下填充 CD4+ 和 CD8+ (解决原Excel首列合并单元格问题)
  fill(Cell_Type, .direction = "down") %>%
  # 将除了第一列之外的所有数据列转换为数值型 (numeric)
  mutate(across(-Cell_Type, as.numeric))

# 4. 宽表转长表 (Melt) 并计算均值
df_long <- df_clean %>%
  pivot_longer(
    cols = -Cell_Type,
    names_to = c("Marker", "Group"),
    names_sep = "_", # 利用下划线将 Tcm_sgNC 拆分为 Marker列 和 Group列
    values_to = "Value"
  )

df_mean <- df_long %>%
  # 按 细胞类型、实验组、Marker 聚合计算平均值
  group_by(Cell_Type, Group, Marker) %>%
  summarize(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# 设定 Y 轴的顺序，保证 sgNC 显示在最上方
# ggplot的Y轴默认是从下往上画，所以这里顺序反转排列一下
df_mean$Group <- factor(df_mean$Group, levels = c("sgTm9sf2#2", "sgTm9sf2#1", "sgNC"))

# 5. 绘制热图 (ggplot2)
p <- ggplot(df_mean, aes(x = Marker, y = Group, fill = Mean_Value)) +
  geom_tile(color = "white", linewidth = 0.5) + # 绘制热图方块并加白色边框
  geom_text(aes(label = sprintf("%.2f", Mean_Value)), 
            color = "white", size = 4) + # 在方块中心显示保留两位小数的数值
  facet_wrap(~ Cell_Type) + # 神仙函数：自动根据 CD4+ 和 CD8+ 分两张子图并排
  scale_fill_viridis_c(
    name = "Mean Percentage (%)", 
    option = "D",
    limits = c(0, 1),           # 设定颜色的上下限，你可以根据需要调整为 3 或 4
    oob = scales::squish        # 关键！将超出 limits 的值“挤压”到边界颜色，而不是显示为NA
  ) + 
  labs(
    title = "T Cell Subpopulations Heatmap",
    x = "Cell State / Marker",
    y = "Experimental Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 15, face = "bold"), # 调整分面标题 (CD4+/CD8+)
    panel.grid = element_blank(), # 去除背景网格线
    axis.text.x = element_text(angle = 45, hjust = 1), # X轴标签倾斜45度防止重叠
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.position = "right"
  )

# 显示图表
print(p)

