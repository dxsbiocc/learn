library(circlize)
library(ComplexHeatmap)

########## 单列热图 ########## 
# 设置热图配色
col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c("#8c510a", "white", "#01665e")
)

# 绘制主热图
exp <- matrix(rnorm(50*10), nrow = 50)
main <- Heatmap(
  exp, name = "main_matrix",
  col = col_fun
  )

# 绘制一列离散型热图
v <- sample(letters[1:2], 50, replace = TRUE)
lnc <- Heatmap(
  v, name = "mat1", 
  # 设置离散型颜色
  col = structure(c("#f46d43", "#66bd63"), name = letters[1:2]),
  top_annotation = HeatmapAnnotation(
    summary = anno_summary(
      height = unit(3, "cm")
    )
  ), 
  width = unit(1, "cm"),
  )

# 绘制一列连续型热图
v <- rnorm(50)
mi <- Heatmap(
  v, name = "mat2", 
  col = col_fun,
  top_annotation = HeatmapAnnotation(
    summary = anno_summary(
      gp = gpar(fill = 2:3), 
      height = unit(3, "cm"))
  ), 
  width = unit(1, "cm")
  )

# 按列添加多个热图
ht_list <- main + lnc + mi

split <- sample(letters[1:2], 50, replace = TRUE)

draw(ht_list, 
     row_split = split, 
     ht_gap = unit(5, "mm"), 
     heatmap_legend_list = list(lgd_boxplot)
     )
     
########## 缩放注释示例 ##########

# 生成表达谱
exp <- matrix(rnorm(100*10), nrow = 100)
# 设置分组
group <- sample(letters[1:3], 100, 
                replace = TRUE, 
                prob = c(1, 5, 10)
                )

panel_fun <- function(index, nm) {
  # 添加绘图 viewport
  pushViewport(viewport(xscale = range(exp), yscale = c(0, 2)))
  # 绘制区域外侧框线
  grid.rect()
  grid.xaxis(gp = gpar(fontsize = 8))
  # 添加箱线图
  grid.boxplot(
    exp[index, ], pos = 1, 
    gp = gpar(fill = index + 6),
    direction = "horizontal"
    )
  popViewport()
}

Heatmap(
  exp, name = "exp", 
  col = col_fun,
  row_split = group,
  height = unit(10, "cm"),
  width = unit(12, "cm"),
  right_annotation = rowAnnotation(
    foo = anno_zoom(
      align_to = group,
      which = "row",
      panel_fun = panel_fun,
      size = unit(2, "cm"),
      gap = unit(1, "cm"),
      width = unit(4, "cm")
    ))
  )
