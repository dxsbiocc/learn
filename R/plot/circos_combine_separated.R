library(circlize)

# 获取当前图形参数，用于复原
op <- par(no.readonly = TRUE)

# 设置一行三列图形排布
par(mar = c(2, 2, 2, 2), mfrow = c(1, 3))

# 1. 第一幅图

plot_circos1 <- function() {
  circos.par("canvas.xlim" = c(-1, 1.5), "canvas.ylim" = c(-1, 1.5), start.degree = -45)
  circos.initialize(sectors = letters[1:4], xlim = c(0, 1))
  # 添加空轨迹
  circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
  # 绘制 a、b 两个扇形区域，并添加文本
  circos.updatePlotRegion(sector.index = "a")
  circos.text(0.5, 0.5, "first one", niceFacing = TRUE, 
              facing = "bending.outside")
  
  circos.updatePlotRegion(sector.index = "b")
  circos.text(0.5, 0.5, "first one", niceFacing = TRUE, 
              facing = "bending.outside")
  highlight.sector(
    c("a", "b"), track.index = 1, 
    col = "#5aae6180"
  )
  circos.clear()
}
plot_circos1()

# 添加外框和轴
box()
axis(side = 1)
axis(side = 2)

# 2. 第二幅图

plot_circos2 <- function() {
  circos.par("canvas.xlim" = c(-1.5, 1), "canvas.ylim" = c(-1.5, 1), start.degree = -45)
  circos.initialize(sectors = letters[1:4], xlim = c(0, 1))
  circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
  # 绘制 c、d 两个扇形区域，并添加文本
  circos.updatePlotRegion(sector.index = "d")
  circos.text(0.5, 0.5, "second one", niceFacing = TRUE,
              facing = "bending.outside")
  circos.updatePlotRegion(sector.index = "c")
  circos.text(0.5, 0.5, "second one", niceFacing = TRUE,
              facing = "bending.outside")
  highlight.sector(
    c("d", "c"), track.index = 1, 
    col = "#9970ab80"
  )
  circos.clear()
}
plot_circos2()

# 添加外框和轴
box()
axis(side = 1)
axis(side = 2)

# 3. 第三幅图

plot_circos1()
# 添加图层
par(new = TRUE)
plot_circos2()

# 复原参数
par(op)
