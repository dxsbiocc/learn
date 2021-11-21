library(tidyverse)
library(ggnewscale)


data <- read_delim('Downloads/gene_sig.txt')


# 根据配对列表生成上、下三角坐标
triangle <- function(pairs, type = "up") {
  # 默认的上三角坐标基
  x = c(0, 0, 1)
  y = c(0, 1, 1)
  # 下三角的坐标基
  if (type == "lower") {
    x = c(0, 1, 1)
    y = c(0, 0, 1)
  }
  # 生成三角矩阵
  mat = do.call(
    rbind,
    apply(pairs, 1, function (row) {
      a = row[1]
      b = row[2]
      data.frame(
        x = x + a,
        y = y + b,
        group = paste(a, b, sep = "-")
      )
    }))
  return(mat)
}

triangle_data <- function(data, row = 1, col = 2) {
  # 生成所有组合
  rows = length(unique(data[[row]]))
  cols = length(unique(data[[col]]))
  pairs = merge(1:rows, 1:cols)
  # 获取上三角坐标
  upper <- triangle(pairs)
  colnames(upper) <- c(paste0("upper.", colnames(upper)[1:2]), "group")
  # 获取下三角坐标
  lower <- triangle(pairs, type = "lower")[1:2]
  colnames(lower) <- paste0("lower.", colnames(lower))
  # 合并坐标
  upper_lower = cbind(upper, lower)
  # 根据分组信息将坐标连接到数据中
  data %>% transmute(across(where(is.factor), ~ as.character(as.numeric(.)))) %>%
    unite("group", row:col, sep = "-") %>%
    cbind(data, .) %>%
    right_join(upper_lower, by = "group")
}

data <- mutate_at(data, 1:2, ~ as.factor(.))

trian_data <- triangle_data(data, row = 2, col = 1)

df <- mutate(trian_data, 
             cor = replace_na(cor, 0),
             p = replace_na(p, 1))

tmp <- data %>% transmute(across(where(is.factor), as.numeric)) %>%
  `names<-`(c("y", "x")) %>%
  cbind(data, .) %>%
  as.data.frame()

points <- do.call(rbind, apply(tmp, 1, function(row) {
  p = as.numeric(row['p'])
  x = as.numeric(row['x'])
  y = as.numeric(row['y'])
  df = data.frame()
  if (p < 0.001) {
    df = rbind(df, data.frame(x = x + 0.9, y = y + 0.5))
  }
  if (p < 0.01) {
    df = rbind(df, data.frame(x = x + 0.9, y = y + 0.3))
  }
  if (p < 0.05) {
    df = rbind(df, data.frame(x = x + 0.9, y = y + 0.1))
  }
  df
}))


ggplot(trian_data) +
  geom_polygon(aes(upper.x, upper.y, fill = abs(cor), group = group), colour = "grey") +
  scale_fill_gradientn(colors = colorRampPalette(c("#1E3163", "#00C1D4", "#FFED99","#FF7600"))(10)) +
  new_scale("fill") +
  geom_polygon(aes(lower.x, lower.y, fill = p, group = group)) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu")) +
  geom_point(data = points, aes(x, y), size = 0.4) +
  scale_x_continuous(breaks = c(1:length(unique(data[[2]]))) + 0.5, expand = c(0,0),
                     labels = sort(unique(data[[2]]))) +
  scale_y_continuous(expand = c(0, 0), breaks = c(1:length(unique(data[[1]]))) + 0.5,
                     labels = sort(unique(data[[1]])), sec.axis = dup_axis()) +
  theme(
    plot.margin = margin(0.5,0.01,0.5,0.01, "cm"),
    axis.title = element_blank(),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5)
  )


