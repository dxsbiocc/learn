# 带标签的并列柱状图
mpg %>% 
  group_by(class, drv) %>% 
  summarise(count = n()) %>%
  ggplot(aes(class, count)) +
  geom_col(aes(fill=drv), position = position_dodge2(preserve = 'single')) +
  geom_text(aes(label=count), 
            position = position_dodge2(width = 0.9, preserve = 'single'), 
            vjust = -0.2, hjust = 0.5)

# 带标签的堆叠柱状图
mpg %>% 
  group_by(class, drv) %>% 
  summarise(count = n()) %>%
  mutate(cumcount = cumsum(count)) %>%
  ggplot(aes(class, count)) +
  geom_col(aes(fill=drv), position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = cumcount), 
            position = position_stack(), # 可以不设置该参数
            vjust = 0.5, hjust=0.5)

# 标签居中的堆叠柱状图
mpg %>% 
  group_by(class, drv) %>% 
  summarise(count = n()) %>%
  mutate(cumcount = cumsum(count), midcount = cumcount - count/2) %>%
  ggplot(aes(class, count)) +
  geom_col(aes(fill = drv), position = position_stack(reverse = TRUE)) +
  geom_text(aes(y = midcount, label = cumcount), hjust=0.5)

# 金字塔图
df <- tibble(
  gene = factor(paste0("gene_", rep(1:16, 2)), levels = paste0("gene_", 16:1)),
  stat = c(seq(-10, -100, -10), seq(-90, -40, 10), seq(10, 100, 10), seq(90, 40, -10)),
  direct = rep(c("down", "up"), each=16)
)

ggplot(df, aes(gene, stat, fill = direct)) + 
  geom_col() +
  coord_flip() + 
  scale_y_continuous(breaks = seq(-100, 100, 20),
                     labels = c(seq(100, 0, -20), seq(20, 100, 20)))
                     
# 偏差图
df <- tibble(
  gene = factor(paste0("gene_", 1:20), levels = paste0("gene_", 20:1)),
  stat = c(seq(100, 10, -10), seq(-10, -100, -10)),
  direct = factor(rep(c("up", "down"), each=10), levels = c("up", "down"))
)

ggplot(df, aes(gene, stat, fill = direct)) + 
  geom_col() +
  coord_flip()
