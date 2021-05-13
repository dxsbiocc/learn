library(ComplexHeatmap)
library(glue)

pdf("~/Downloads/hm.pdf")

set.seed(123)
mat <- matrix(
  rnorm(50*50), 
  nrow = 50
  )

group <- c(
  "MAPK",
  "PI3K-Akt",
  "ErbB",
  "Cell cycle",
  "Apoptosis"
)

split = rep(1:5, each = 10)

ha <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(
    gp = gpar(fill = 2:6), 
    labels = group
    )
)

Heatmap(
  mat, name = "mat", 
  column_split = split, 
  top_annotation = ha, 
  column_title = NULL
  )

block_group_anno <- function(group, empty_anno, gp = gpar(),
                             label = NULL, label_gp = gpar()) {
  # 获取最左侧 viewport
  seekViewport(glue(
    'annotation_{anno}_{slice}', 
    slice = min(group),
    anno = empty_anno)
    )
  # 获取左下角坐标点
  loc1 <- deviceLoc(
    x = unit(0, "npc"), 
    y = unit(0, "npc")
  )
  # 获取最右侧 viewport
  seekViewport(glue(
    'annotation_{anno}_{slice}', 
    slice = max(group),
    anno = empty_anno)
  )
  # 获取右上角坐标点
  loc2 <- deviceLoc(
    x = unit(1, "npc"), 
    y = unit(1, "npc")
  )
  # 切换到全局 viewport
  seekViewport("global")
  # 绘制矩形
  grid.rect(
    loc1$x, loc1$y,
    width = loc2$x - loc1$x, 
    height = loc2$y - loc1$y, 
    just = c("left", "bottom"), 
    gp = gp
  )
  # 如果传递了标签，则添加标签
  if (!is.null(label)) {
    grid.text(
      label, 
      x = (loc1$x + loc2$x) * 0.5, 
      y = (loc1$y + loc2$y) * 0.5,
      gp = label_gp
    )
  }
}
# 将前三个热图块作为一组
block_group_anno(1:3, "empty", gp = gpar(fill = "red"), label = "Signal transduction")
# 后两个作为一组
block_group_anno(4:5, "empty", gp = gpar(fill = "green"), label = "Cellular Processes")
dev.off()
