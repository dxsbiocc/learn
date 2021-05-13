library(ComplexHeatmap)


########## decorate text ##########
group <- list(
  A = "Cell cycle",
  B = "DNA replication",
  C = "Mismatch repair",
  D = "MAPK signaling pathway"
)

ha = rowAnnotation(
  foo = anno_empty(
    border = FALSE, 
    # 计算空白注释的宽度
    width = max_text_width(unlist(group)) + unit(4, "mm"))
  )

Heatmap(matrix(
  rnorm(1000), nrow = 100), 
  name = "mat", 
  # 分 4 块
  row_km = 4, 
  right_annotation = ha
  )

for(i in 1:4) {
  decorate_annotation(
    "foo", 
    # 选择热图块
    slice = i, {
    # 添加颜色框
    grid.rect(
      x = 0, 
      width = unit(2, "mm"), 
      gp = gpar(
        fill = rainbow(4)[i], 
        col = NA
        ), 
      just = "left"
      )
    # 绘制文本
    grid.text(
      group[[i]], 
      x = unit(4, "mm"), 
      gp = gpar(
        col = rainbow(4)[i]
      ),
      just = "left")
  })
}

########## decorate points ##########

ha <- HeatmapAnnotation(
  foo = anno_empty(
    border = TRUE, 
    # 固定注释的高度
    height = unit(3, "cm"))
  )

ht <- Heatmap(
  matrix(rnorm(100), nrow = 10), 
  name = "mat", 
  top_annotation = ha
  )

# 先绘制热图
ht <- draw(ht)
# 获取热图的列顺序
co <- column_order(ht)
# 生成 10 个均匀分布的随机值
value <- runif(10)
# 将空白装饰一下
decorate_annotation("foo", {
  # 对应于矩阵的列
  x = 1:10
  # 根据列顺序重排这些随机值，用于设置 y 轴的值
  value = value[co]
  # 添加 viewport
  pushViewport(
    viewport(
      # x 轴范围
      xscale = c(0.5, 10.5), 
      # y 轴范围
      yscale = c(0, 1))
    )
  # 在中间添加一条水平虚线
  grid.lines(
    c(0.5, 10.5), c(0.5, 0.5), 
    gp = gpar(lty = 2),
    default.units = "native"
    )
  # 添加点
  grid.points(
    x, value, 
    pch = 16, 
    size = unit(2, "mm"),
    gp = gpar(
      # 虚线上下的值设置不同颜色
      col = ifelse(value > 0.5, "red", "blue")), 
    default.units = "native"
    )
  # 设置 y 轴断点
  grid.yaxis(at = c(0, 0.5, 1))
  popViewport()
})
