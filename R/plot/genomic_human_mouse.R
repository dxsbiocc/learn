library(circlize)


# 获取人类和小鼠的包含条带信息的基因组数据
human_cytoband <- read.cytoband(species = "hg19")$df
mouse_cytoband <- read.cytoband(species = "mm10")$df
# 为染色体添加前缀
human_cytoband[ ,1] <- paste0("human_", human_cytoband[, 1])
mouse_cytoband[ ,1] <- paste0("mouse_", mouse_cytoband[, 1])
# 合并两份数据
cytoband <- rbind(human_cytoband, mouse_cytoband)

# 设置染色体顺序
chromosome.index <- c(
  paste0("human_chr", c(1:22, "X", "Y")), 
  rev(paste0("mouse_chr", c(1:19, "X", "Y")))
)

# 获取人类和小鼠的染色体长度信息
human_chromInfo <- read.chromInfo(species = "hg19")$df
mouse_chromInfo <- read.chromInfo(species = "mm10")$df
human_chromInfo[ ,1] <- paste0("human_", human_chromInfo[, 1])
mouse_chromInfo[ ,1] <- paste0("mouse_", mouse_chromInfo[, 1])
chromInfo <- rbind(human_chromInfo, mouse_chromInfo)
chromInfo[, 1] <- factor(chromInfo[ ,1], levels = chromosome.index)

# 设置间距
circos.par(gap.after = c(rep(1, 23), 5, rep(1, 20), 5))
# 初始化布局，不添加图形
circos.genomicInitialize(chromInfo, plotType = NULL)
# 添加数字染色体号
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2] + mm_y(2),
      gsub(".*chr", "", CELL_META$sector.index),
      cex = 0.6,
      niceFacing = TRUE
    )
  },
  track.height = mm_h(1),
  cell.padding = c(0, 0, 0, 0),
  bg.border = NA
)
# 添加分组颜色轨迹
highlight.chromosome(
  paste0("human_chr", c(1:22, "X", "Y")),
  col = "#66c2a5", track.index = 1
)
highlight.chromosome(
  paste0("mouse_chr", c(1:19, "X", "Y")),
  col = "#fc8d62", track.index = 1
)

# 添加 ideogram
circos.genomicIdeogram(cytoband)

# 创建随机数据
human_df <- generateRandomBed(200, species = "hg19")
mouse_df <- generateRandomBed(200, species = "mm10")
human_df[, 1] <- paste0("human_", human_df[, 1])
mouse_df[, 1] <- paste0("mouse_", mouse_df[, 1])
df <- rbind(human_df, mouse_df)
# 添加点图
circos.genomicTrack(
  df,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, col = rand_color(1), cex = 0.5, ...)
  }
)

# 添加人类与小鼠基因组之间的连接
human_mid <- data.frame(
  chr = paste0("human_chr", 1:19),
  mid = round((human_chromInfo[1:19, 2] + human_chromInfo[1:19, 3]) / 2)
)
mouse_mid <- data.frame(
  chr = paste0("mouse_chr", 1:19),
  mid = round((mouse_chromInfo[1:19, 2] + mouse_chromInfo[1:19, 3]) / 2)
)
circos.genomicLink(human_mid, mouse_mid, col = rand_color(19))
circos.clear()
# 添加注释
text(-0.9,-0.8, "Human\ngenome")
text(0.9, 0.8, "Mouse\ngenome")
