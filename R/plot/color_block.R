library(tidyverse)
library(cowplot)

# 色块图
############## 华夫饼图 ##############
# 堆积型
df <- tibble(
  x = rep(1:10, 10),
  y = rep(1:10, each=10),
  class = sort(sample(mpg$class, 100))
)

sample_n(df, 67) %>% 
  arrange(x, y) %>% 
  group_by(x) %>% 
  mutate(y = 1:n()) %>%
  ggplot(aes(x, y, fill = class)) +
  geom_tile(colour = "white") +
  # geom_point(size = 12, shape = 21) +
  coord_fixed() +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
 
# 百分比堆积型
ggplot(df, aes(x, y, fill = class)) +
  geom_tile(colour = "white") +
  scale_y_continuous(trans = "reverse") +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

############## 马赛克图 ##############
df <- tibble(
  gene = rep(LETTERS[1:5], 4),
  status = rep(c("alpha", "beta", "gamma", "delta"), each = 5),
  value = sample(1:100, 20),
  percent = rep(c(10, 30, 30, 20, 10), 4)
)

df %>% 
  group_by(status) %>%
  mutate(xmax = cumsum(percent), xmin = xmax - percent) %>%
  group_by(gene) %>% 
  mutate(ytmp = value * 100 / sum(value), ymax = cumsum(ytmp), ymin = ymax - ytmp) %>%
  mutate(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2) %>%
  ggplot() + 
  geom_rect(aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = status),
            colour = "black") +
  geom_text(aes(x, y, label = paste0(round(ytmp, 2), "%"))) +
  geom_text(aes(x = x, y = 103, label = gene)) +
  theme(panel.background = element_blank())

############## 瀑布图 ##############
# > data
# A tibble: 2,432 x 3
#   sample   gene  MutFunc          
#   <chr>    <chr> <chr>            
# 1 C1803704 ATM   nonsynonymous SNV
# 2 C1803704 BRAF  nonsynonymous SNV
# 3 C1803704 BRCA1 synonymous SNV   
# 4 C1803704 BRCA1 nonsynonymous SNV
# 5 C1803704 BRCA2 nonsynonymous SNV
# 6 C1803704 BRCA2 synonymous SNV   
# 7 C1803704 BRCA2 stopgain         
# 8 C1803704 BRD4  nonsynonymous SNV
# 9 C1803704 EOMES nonsynonymous SNV
# 10 C1803704 EPCAM nonsynonymous SNV
# … with 2,422 more rows

genes <- count(data, gene) %>%
  top_n(n = 20, wt = n) %>%
  mutate(percent = round(n * 100 / sum(n), 1)) %>%
  arrange(desc(n))
  
samples <- subset(data, gene %in% genes$gene) %>%
  count(sample) %>% arrange(desc(n)) %>%
  rename(num = n)

df <- inner_join(data, genes) %>%
  mutate(gene = factor(gene, levels = rev(genes$gene)),
         sample = factor(sample, levels = samples$sample)) %>%
  inner_join(samples)

p1 <- ggplot(df) +
  geom_tile(aes(x = sample, y = gene, fill = MutFunc)) +
  # geom_text(aes(x = -4, y = gene, label = percent), data = genes) +
  # 图例行数的调整放到 fill，单独用 guides 无效
  scale_fill_discrete(guide = guide_legend(nrow = 3)) +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 4),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = 'white', colour = 'black')
        )
# 样本的突变基因数目条形图
p2 <- ggplot(df) +
  geom_bar(aes(x = sample, fill = MutFunc)) +
  # 使用 expand 删除数据与轴之间的空隙
  scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0, 25), 
                     expand = expansion(mult = 0, add = 0)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.length.y.left = unit(.25, "cm"),
    axis.line.y.left = element_line(colour = "black"),
  )
# 突变基因的突变频数条形图
p3 <- ggplot(df) +
  geom_bar(aes(y = gene, fill = MutFunc)) +
  scale_x_continuous(position = "top", breaks = seq(0, 280, 70),
                     limits = c(0, 280),
                     expand = expansion(mult = 0, add = 0)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.length.x.top = unit(.25, "cm"),
    axis.line.x.top = element_line(colour = "black")
  )

# 组合图形
p12 <- plot_grid(p2, p1, nrow = 2, align = "v",
          rel_heights = c(1, 5))

pr <- plot_grid(NULL, p3, NULL, nrow = 3, align = "v",
                rel_heights = c(0.95, 5, 1.05))

plot_grid(p12, pr, ncol = 2, align = "h",
          rel_widths = c(5, 1))
