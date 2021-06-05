library(circlize)

# 获取当前图形参数，用于复原
op <- par(no.readonly = TRUE)

# 设置一行三列图形排布
par(mar = c(2, 2, 2, 2), mfrow = c(1, 3))
圆形
# 1. 第一幅图
plot_circos1 <- function() {
  sectors <- letters[1:4]
  circos.initialize(sectors = sectors, xlim = c(0, 1))
  circos.trackPlotRegion(
    ylim = c(0, 1), 
    panel.fun = function(x, y) {
      circos.text(
        0.5, 0.5, "inner circos", col = 6,
        niceFacing = TRUE, facing = "bending.outside"
      )
    }
  )
  circos.clear()
}
plot_circos1()

# 添加外框
box()
# 添加坐标轴
axis(side = 1)
axis(side = 2)

# 2. 第二幅图
# 设置更大的 canvas.xlim 和 canvas.ylim
# 同样的 xlim  = c(0, 1)，会绘制更小的圆形
plot_circos2 <- function() {
  circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2))
  sectors <- letters[1:3]
  circos.initialize(sectors = sectors, xlim = c(0, 1))
  circos.trackPlotRegion(
    ylim = c(0, 1), 
    panel.fun = function(x, y) {
      circos.text(
        0.5, 0.5, "inner circos", col = 7,
        niceFacing = TRUE, facing = "bending.outside"
        )
    }
  )
  circos.clear()
}
plot_circos2()

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
