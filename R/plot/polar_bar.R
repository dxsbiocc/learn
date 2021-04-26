####################  径向柱状图  ####################
df <- read_delim("~/Downloads/coad_caseccc_2015/data_mutations_mskcc.txt", delim = "\t")

select(df, Tumor_Sample_Barcode, Hugo_Symbol) %>%
  count(df, Hugo_Symbol) %>%
  filter(n > 3) %>%
  arrange(n) %>%
  ggplot(aes(Hugo_Symbol, n, fill = Hugo_Symbol)) +
  # ggplot(aes(factor(Hugo_Symbol, levels = Hugo_Symbol), n, fill = Hugo_Symbol)) +
  geom_col() +
  geom_text(aes(y = n - 2, label = n), colour = "white") +
  coord_polar(start = 0) +
  ylim(c(-10, 35)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "grey80",size=.25),
    axis.text.x = element_text(size = 9, colour="black", angle = seq(-10, -350, length.out = 27)),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )
####################  分组径向柱状图  ####################

# 设置空白柱子的个数
empty_bar = 2
# 自定义突变类型
mut_type <- c("Ins", "Del", "Mismatch", "Silent")
# 构造数据
data <- tibble(
  gene=paste( "Gene ", seq(1,60), sep=""),
  group=c(rep('Ins', 10), rep('Mismatch', 30), rep('Del', 14), rep('Silent', 6)) ,
  value=sample(seq(10,100), 60, replace=T)
) %>%
  # 添加 NA 数据，用于在分组之间绘制空白柱形
  add_row(tibble(
    gene = rep(NA, empty_bar * length(mut_type)),
    group = rep(mut_type, empty_bar),
    value = gene
  )) %>%
  mutate(group = factor(group, levels = mut_type)) %>%
  # 排序，为了让统一分组绘制在一起
  arrange(group)

# 构造唯一标识，用作 x 轴，并按该顺序绘制
data$id = 1:nrow(data)
# 添加显示文本的角度
angle <- 90 - 360 * (data$id - 0.5) / nrow(data)
# 添加内圈注释
base_anno <- group_by(data, group) %>%
  summarise(start = min(id), end = max(id) - empty_bar) %>%
  mutate(mid = (start + end) / 2)
  
ggplot(data, aes(id, value, fill = group)) +
  geom_col(position = position_dodge2()) +
  geom_text(aes(y = value + 18, label = gene), size = 2.5, alpha = 0.6, 
            angle = ifelse(angle < -90, angle+180, angle)) +
  # 内圈注释
  geom_segment(data = base_anno, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "grey40") +
  geom_text(data = base_anno, aes(x = mid, y = -18, label = group), 
            angle = c(-26, -100, -50, 26), colour = "grey40") +
  coord_polar() +
  ylim(-100,120) +
  theme(
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
  )

####################  径向热图  ####################
bit_data <- read_csv("~/Downloads/bit_data.csv")

group_by(bit_data, year, month) %>%
  summarise(value = mean(High), .groups = "drop") %>%
  ggplot(aes(factor(month), year, fill = value)) +
  geom_tile(width = 1, colour = "white") +
  coord_polar() +
  ylim(c(2010, 2020)) +
  scale_fill_gradientn(colours = rainbow(10)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "grey80",size=.25),
    axis.text.x = element_text(size = 9, colour="black", angle = seq(-10, -350, length.out = 12)),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank()
  )
 
#################### 径向圈圈图 ####################
group_by(bit_data, year, month) %>%
  summarise(value = mean(High)) %>%
  mutate(
    xmin = month,
    xmax = month + 1,
    ymin = (year - 2015) * 10 + 1,
    ymax = ymin + sample(1:5, n(), replace = TRUE)
    ) %>%
  ggplot(aes(fill = value)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  scale_x_continuous(breaks = seq(1.5, 12.5, 1), labels = month.name) +
  scale_fill_gradientn(colours = rainbow(10)) +
  coord_polar() +
  ylim(c(-5, 40)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "grey80",size=.25),
    axis.text.x = element_text(size = 9, colour="black", angle = seq(-10, -350, length.out = 12)),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank()
  )
