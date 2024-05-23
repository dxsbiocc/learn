library(ggplot2)
library(gtable)
library(grid)


# 反转标题 grobs 的对齐方式和位置
hinvert_title_grob <- function(grob){
    # 交换宽度
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    # 在一些版本的ggplot2中，grob 对象可能不包含 vp 或 layout 组件
    if (!is.null(grob$vp) && !is.null(grob$vp[[1]]$layout)) {
        grob$vp[[1]]$layout$widths[1] <- widths[3]
        grob$vp[[1]]$layout$widths[3] <- widths[1]
    }
    # 修改对齐方式
    grob$children[[1]] <- editGrob(
        grob$children[[1]], 
        hjust = 1 - grob$children[[1]]$hjust, 
        vjust = 1 - grob$children[[1]]$vjust, 
        x = unit(1, "npc") - grob$children[[1]]$x
    )
    grob
}
# 左侧添加Y轴
add_yaxis_left <- function(g1, g2) {
    # 获取左侧Y轴标签的位置
    pos <- c(subset(g1$layout, name == "ylab-l", select = t:r))
    index <- which(g2$layout$name == "axis-l")
    yaxis <- g2$grobs[[index]]
    # 添加列并插入Y轴
    g <- gtable_add_cols(g1, unit(3, "mm"), pos$l - 1)
    g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$l - 1)
    g <- gtable_add_grob(g, yaxis, pos$t, pos$l, pos$b, pos$l, clip = "off")
    # 获取左侧Y轴标签
    index <- which(g2$layout$name == "ylab-l")
    ylab <- g2$grobs[[index]]
    # 插入Y轴标签
    g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$l - 1)
    g <- gtable_add_grob(g, ylab, pos$t, pos$l, pos$b, pos$l, clip = "off")
    g
}
# 右侧添加Y轴
add_yaxis_right <- function(g1, g2, pos) {
    # 获取右侧Y轴标签
    index <- which(g2$layout$name == "ylab-l")
    ylab <- g2$grobs[[index]]
    ylab <- hinvert_title_grob(ylab)
    # 添加列并插入Y轴标签
    g <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pos$r)
    g <- gtable_add_grob(g, ylab, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "ylab-r")
    # 获取右侧Y轴
    index <- which(g2$layout$name == "axis-l")
    yaxis <- g2$grobs[[index]]
    # 调整Y轴线位置
    yaxis$children[[1]]$x <- unit.c(unit(0, "npc"), unit(0, "npc"))
    # 获取刻度
    ticks <- yaxis$children[[2]]
    # 调整刻度线位置
    for (i in seq_along(ticks$grobs)) {
        # 刻度标签
        if (inherits(ticks$grobs[[i]], "titleGrob")) {
            ticks$grobs[[i]] <- hinvert_title_grob(ticks$grobs[[i]])
        # 刻度线
        } else if (inherits(ticks$grobs[[i]], "polyline")) {
            ticks$grobs[[i]]$x <- ticks$grobs[[i]]$x - unit(1, "npc") + unit(3, "pt")
        }
    }
    # 反转刻度线和刻度标签
    ticks$widths <- rev(ticks$widths)
    ticks$grobs <- rev(ticks$grobs)
    
    yaxis$children[[2]] <- ticks
    # 添加列并插入Y轴
    g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l] + unit(3, "mm"), pos$r)
    g <- gtable_add_grob(g, yaxis, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "axis-r")
    g
}
# 添加Y轴函数，根据偏移量判断添加在左侧或右侧
add_yaxis <- function(g1, g2, offset = 0) {
    pos <- c(subset(g1$layout, name == "panel", select = t:r))
    # 添加主绘图区域
    g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], 
                          pos$t, pos$l, pos$b * ((offset - 2) * 0.00001 + 1), pos$l)
    # 根据偏移量判断添加在左侧或右侧
    if (offset > 3 && offset %% 2 == 0) {
        g1 <- add_yaxis_left(g1, g2)
    } else {
        g1 <- add_yaxis_right(g1, g2, pos)
    }
    g1
}
# 绘制多Y轴图形的函数
plot_multi_yaxis <- function(..., right_label_reverse = TRUE) {
    args <- list(...)
    my_theme <- theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA))
    len <- length(args)
    args[[1]] <- args[[1]] + my_theme
    g <- ggplotGrob(args[[1]])
    for (i in len:2) { 
        if (i < 4 || i %% 2 && right_label_reverse) {
            # 旋转轴标签
            args[[i]] <- args[[i]] + 
                theme(axis.title.y = element_text(angle = 270))
        }
        args[[i]] <- args[[i]] + my_theme
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
