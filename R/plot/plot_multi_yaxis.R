library(ggplot2)
library(gtable)
library(grid)


hinvert_title_grob <- function(grob){
  # 交换宽度
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
  
  # 修改对齐
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
  grob
}

add_yaxis_left <- function(g1, g2) {
  # 将坐标轴添加到左侧
  # 添加坐标轴
  pos <- c(subset(g1$layout, name == "ylab-l", select = t:r))
  index <- which(g2$layout$name == "axis-l")
  yaxis <- g2$grobs[[index]]
  # 先添加 3mm 间距
  g <- gtable_add_cols(g1, unit(3, "mm"), pos$l - 1)
  # 再添加轴
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$l - 1)
  g <- gtable_add_grob(g, yaxis, pos$t, pos$l, pos$b, pos$l, clip = "off", name = "axis-l")
  # 添加轴标签
  index <- which(g2$layout$name == "ylab-l")
  ylab <- g2$grobs[[index]]
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$l - 1)
  g <- gtable_add_grob(g, ylab, pos$t, pos$l, pos$b, pos$l, clip = "off", name = "ylab-l")
  g
}

add_yaxis_right <- function(g1, g2, pos) {
  # 将坐标轴添加到右侧
  # ============ 2. 轴标签 ============ #
  index <- which(g2$layout$name == "ylab-l")
  ylab <- g2$grobs[[index]]
  ylab <- hinvert_title_grob(ylab)
  # 添加轴标签
  g <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pos$r)
  g <- gtable_add_grob(g, ylab, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "ylab-r")
  # ============ 3. 轴设置 ============ #
  index <- which(g2$layout$name == "axis-l")
  yaxis <- g2$grobs[[index]]
  # 将 Y 轴线移动到最左边
  yaxis$children[[1]]$x <- unit.c(unit(0, "npc"), unit(0, "npc"))
  # 交换刻度线和刻度标签
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  # 移动刻度线
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + unit(3, "pt")
  # 刻度标签位置转换和对齐
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  yaxis$children[[2]] <- ticks
  # 添加轴，unit(3, "mm") 增加轴间距
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l] + unit(3, "mm"), pos$r)
  g <- gtable_add_grob(g, yaxis, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "axis-r")
  g
}

add_yaxis <- function(g1, g2, offset = 0) {
  # ============ 1. 主绘图区 ============ #
  # 获取主绘图区域
  pos <- c(subset(g1$layout, name == "panel", select = t:r))
  # 添加图形
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], 
                       pos$t, pos$l, pos$b * ((offset - 2) * 0.00001 + 1), pos$l)
  if (offset > 3 && offset %% 2 == 0) {
    g <- add_yaxis_left(g, g2)
  } else {
    g <- add_yaxis_right(g, g2, pos)
  }
  g
}

# 接受可变参数，可添加多个 Y 轴
plot_multi_yaxis <- function(..., right_label_reverse = TRUE) {
  args <- list(...)
  my_theme <- theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA))
  len <- length(args)
  args[[1]] <- args[[1]] + my_theme
  g <- ggplotGrob(args[[1]])
  for (i in len:2) { 
    if (i < 4 || i %% 2 && right_label_reverse) {
      # 为轴标签添加旋转
      args[[i]] <- args[[i]] + 
        theme(axis.title.y = element_text(angle = 270))
    }
    args[[i]] <- args[[i]] + my_theme
    # 获取 gtable 对象
    g2 <- ggplotGrob(args[[i]])
    g <- add_yaxis(g, g2, offset = i)
  }
  # 绘制图形
  grid.newpage()
  grid.draw(g)
}

colors <- c('#5470C6', '#91CC75', '#EE6666')
data <- data.frame(
  category = factor(substr(month.name, 1, 3), levels = substr(month.name, 1, 3)),
  Evaporation = c(2.0, 4.9, 7.0, 23.2, 25.6, 76.7, 135.6, 162.2, 32.6, 20.0, 6.4, 3.3),
  Precipitation = c(2.6, 5.9, 9.0, 26.4, 28.7, 70.7, 175.6, 182.2, 48.7, 18.8, 6.0, 2.3),
  Temperature = c(2.0, 2.2, 3.3, 4.5, 6.3, 10.2, 20.3, 23.4, 23.0, 16.5, 12.0, 6.2)
)


p1 <- ggplot(data, aes(category, Evaporation)) + 
  geom_col(fill = colors[1], width = 0.3, position = position_nudge(x = -0.2)) + 
  labs(x = "month", y = "Evaporation(ml)") +
  scale_y_continuous(limits = c(0, 250), expand = c(0,0)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[1]), 
        axis.ticks.y = element_line(color = colors[1]), 
        axis.title.y = element_text(color = colors[1]), 
        axis.line.y = element_line(color = colors[1]), 
        axis.line.x = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p2 <- ggplot(data, aes(category, Precipitation)) + 
  geom_col(fill = colors[2], width = 0.3, position = position_nudge(x = 0.2)) + 
  labs(x = "month", y = "Precipitation(ml)") +
  scale_y_continuous(limits = c(0, 250), expand = c(0,0))  +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[2]), 
        axis.ticks.y = element_line(color = colors[2]), 
        axis.title.y = element_text(color = colors[2]), 
        axis.line.y = element_line(color = colors[2]), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p3 <- ggplot(data, aes(category, Temperature, group = 1)) + 
  geom_line(colour = colors[3]) + 
  geom_point(aes(colour = colors[3]), fill = "white", shape = 21, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 25), expand = c(0,0)) +
  labs(x = "month", y = expression(paste("Temperature (", degree, " C)"))) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[3]), 
        axis.ticks.y = element_line(color = colors[3]), 
        axis.title.y = element_text(color = colors[3]), 
        axis.line.y = element_line(color = colors[3]), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


plot_multi_yaxis(p1, p2, p3)
