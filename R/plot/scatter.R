library(tidyverse)

# 克利夫兰点图
group_by(mpg, manufacturer) %>%
  summarise(displ = max(displ)) %>%
  ggplot(aes(displ, manufacturer)) +
  geom_point(shape = 21, colour = "black", fill="#dd3497", size = 3)
  
# 棒棒糖图
group_by(mpg, manufacturer) %>%
  summarise(displ = max(displ)) %>%
  ggplot(aes(displ, manufacturer)) +
  geom_segment(aes(x = 0, xend = displ, 
                   y = manufacturer, yend = manufacturer)) +
  geom_point(shape = 21, colour = "black", fill="#dd3497", size = 3)
  
# 哑铃图
group_by(mpg, manufacturer, year) %>%
  summarise(displ = max(displ)) %>%
  ggplot(aes(displ, manufacturer)) +
  geom_line(aes(group = manufacturer)) +
  geom_point(aes(fill = factor(year)), shape = 21, colour = "black", size = 3)
  
# 残差分析图
ggplot(dsamp, aes(carat, price)) +
  geom_segment(aes(xend = carat, yend = predicted), alpha = 0.2) +
  geom_point(aes(fill = abs_rd, size = abs_rd), shape = 21, colour = "black") +
  geom_smooth(method = "lm", se = FALSE, colour = "lightgrey") +
  geom_point(aes(y = predicted), shape = 1) +
  scale_fill_continuous(low = "#fb9a99", high = "#b2df8a")
