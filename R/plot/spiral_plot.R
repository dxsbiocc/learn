library(tidyverse)

bit_data <- read_csv("~/Downloads/bit_data.csv")

#################### 螺旋柱状图 ####################
data <- group_by(bit_data, year, month) %>%
  summarise(value = mean(High)) %>%
  mutate(
    xmin = month,
    xmax = month + 1,
    ymin = (year - 2015) * 12 + month,
    ymax = ymin + runif(n(), min = 3, max = 10)
    # ymax = ymin + 12  # 螺旋热图
  ) %>%
  # 根据每行的四个坐标，构建成 4 行，代表 4 个点
  # 因为 geom_polygon 是逆时针连接起来的，所以点
  # 的顺序也要依次排列
  rowwise() %>%
  do(with(., tibble(
    year = year,
    month = month,
    value = value,
    x = c(xmin, xmax, xmax, xmin),
    y = c(ymin, ymin + 1, ymax + 1, ymax)
  )))

# 在这里，设置 group 参数非常重要。我们要每个月份
# 作为一组，绘制一个平行四边形，所以需要根据年份
# 与月份两列来进行分组
ggplot(data, aes(x, y, group = paste(year, month))) +
  geom_polygon(aes(fill = value), colour = "black") +
  coord_polar() +
  scale_x_continuous("month", breaks = seq(1.5, 12.5, 1), labels = month.name) +
  scale_fill_gradientn(values = seq(0,1,0.2), colours = c('cyan','blue','green','orange','red')) +
  ylim(c(-6, 60)) +
  theme(
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank()
  )

#################### 螺旋面积图 ####################
df <- tibble(
  date = seq(as.Date("2015-01-01"), as.Date("2019-12-31"), "days"),
  value = runif(length(date), min = 50, max = 300)
) %>%
  mutate(year = year(date), month = month(date), day = yday(date)) %>%
  # 去除闰年多出的一天
  filter(day != 366) %>%
  mutate(
    # 构造 y 轴梯度
    ymin = (year - 2015) * 364 + day, 
    ymax = ymin + value
  )

ggplot(df, aes(x = day, group = year)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "orange") +
  # geom_linerange(aes(ymin = ymin, ymax = ymax, colour = value)) # 渐变色
  geom_line(aes(y = ymin)) +
  geom_line(aes(y = ymax), colour = "grey40") +
  coord_polar(start = 1) +
  scale_x_continuous(
    breaks = c(1,31,59,90,120,151,181,212,243,273,304,334),
    labels = month.name
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major.x = element_line(colour = "grey20", size = 0.25),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank()
  )
