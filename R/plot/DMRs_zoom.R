library(circlize)
library(RColorBrewer)

load(system.file(package = "circlize", "extdata", "tagments_WGBS_DMR.RData"))

# 设置颜色映射
chr_bg_color <- paste0(sample(
  c(brewer.pal(9, "Set1"), 
    brewer.pal(8, "Set2"), 
    brewer.pal(12, "Set3")), 
  size = 22), "80")

names(chr_bg_color) <- paste0("chr", 1:22)
# 绘制整个基因组
f1 <- function() {
  circos.par(gap.after = 2, start.degree = 90)
  circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), 
                                plotType = c("ideogram", "labels"), ideogram.height = 0.03)
}

f2 <- function() {
  circos.par(cell.padding = c(0, 0, 0, 0), gap.after = c(rep(1, nrow(tagments)-1), 10))
  circos.genomicInitialize(tagments, plotType = NULL)
  # 绘制检测区域中的 DMR 的数据
  circos.genomicTrack(
    DMR1, ylim = c(-0.6, 0.6), 
    panel.fun = function(region, value, ...) {
      # 添加虚线
      for (h in seq(-0.6, 0.6, by = 0.2)) {
        circos.lines(
          CELL_META$cell.xlim, c(h, h),
          lty = 3, col = "#AAAAAA"
        )
      }
      circos.lines(
        CELL_META$cell.xlim, c(0, 0), lty = 3, 
        col = "#888888"
      )
      # 根据 DMR 值的正负设置不同颜色的点
      circos.genomicPoints(
        region, value, pch = 16, cex = 0.5,
        col = ifelse(value[[1]] > 0, "#E41A1C", "#377EB8")
      )
    },
    # 设置背景色
    bg.col = chr_bg_color[tagments$chr],
    track.margin = c(0.02, 0)
  )
  # 添加 y 轴刻度和标签
  circos.yaxis(
    side = "left",
    at = seq(-0.6, 0.6, by = 0.3),
    sector.index = get.all.sector.index()[1],
    labels.cex = 0.4
  )
  # 添加最内层颜色圆环，标注检测区域所属染色体
  circos.track(
    ylim = c(0, 1),
    track.height = mm_h(2),
    bg.col = add_transparency(chr_bg_color[tagments$chr], 0)
  )
}

circos.nested(
  f1, f2, correspondance, 
  connection_col = chr_bg_color[correspondance[[1]]]
)
