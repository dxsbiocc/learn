library(circlize)
library(gridBase)
library(RColorBrewer)
library(ComplexHeatmap)


# 读取数据
res_list <- readRDS("Downloads/meth.rds")
# 分配数据变量
type <- res_list$type
mat_meth <- res_list$mat_meth
mat_expr <- res_list$mat_expr
direction <- res_list$direction
cor_pvalue <- res_list$cor_pvalue
gene_type <- res_list$gene_type
anno_gene <- res_list$anno_gene
dist <- res_list$dist
anno_enhancer <- res_list$anno_enhancer
# k-means 聚类
km <- kmeans(mat_meth, centers = 5)$cluster
# 为数据分配颜色
col_meth <- colorRamp2(c(0, 0.5, 1), c("#a6611a", "#f5f5f5", "#018571"))
col_direction <- c("hyper" = "red", "hypo" = "blue")
col_expr <- colorRamp2(c(-2, 0, 2), c("#d01c8b", "#f7f7f7", "#4dac26"))
col_pvalue <- colorRamp2(c(0, 2, 4), c("#f1a340", "#f7f7f7", "#998ec3"))

col_gene_type <- structure(brewer.pal(length(unique(gene_type)), "Set3"), names = unique(gene_type))
col_anno_gene <- structure(brewer.pal(length(unique(anno_gene)), "Set1"), names = unique(anno_gene))

col_dist <- colorRamp2(c(0, 10000), c("#ef8a62", "#67a9cf"))
col_enhancer <- colorRamp2(c(0, 1), c("#fc8d59", "#99d594"))

# 创建连接数据
df_link <- data.frame(
  from_index = sample(nrow(mat_meth), 20),
  to_index = sample(nrow(mat_meth), 20)
)

# 绘制圆形热图
circlize_plot <- function() {
  circos.heatmap(mat_meth, split = km, col = col_meth, track.height = 0.12)
  circos.heatmap(direction, col = col_direction, track.height = 0.01)
  circos.heatmap(mat_expr, col = col_expr, track.height = 0.12)
  circos.heatmap(cor_pvalue, col = col_pvalue, track.height = 0.01)
  circos.heatmap(gene_type, col = col_gene_type, track.height = 0.01)
  circos.heatmap(anno_gene, col = col_anno_gene, track.height = 0.01) 
  circos.heatmap(dist, col = col_dist, track.height = 0.01)
  circos.heatmap(anno_enhancer, col = col_enhancer, track.height = 0.03)
  
  # 添加连接线
  for(i in seq_len(nrow(df_link))) {
    circos.heatmap.link(
      df_link$from_index[i], df_link$to_index[i], col = rand_color(1))
  }
  circos.clear()
}

# 设置图例
lgd_meth <- Legend(title = "Methylation", col_fun = col_meth)
lgd_direction <- Legend(
  title = "Direction", at = names(col_direction), 
  legend_gp = gpar(fill = col_direction)
)
lgd_expr <- Legend(title = "Expression", col_fun = col_expr)
lgd_pvalue <- Legend(
  title = "P-value", col_fun = col_pvalue, at = c(0, 2, 4), 
  labels = c(1, 0.01, 0.0001)
)
lgd_gene_type <- Legend(
  title = "Gene type", at = names(col_gene_type), 
  legend_gp = gpar(fill = col_gene_type)
)
lgd_anno_gene <- Legend(
  title = "Gene anno", at = names(col_anno_gene), 
  legend_gp = gpar(fill = col_anno_gene)
)
lgd_dist <- Legend(
  title = "Dist to TSS", col_fun = col_dist, 
  at = c(0, 5000, 10000), labels = c("0kb", "5kb", "10kb")
)
lgd_enhancer <- Legend(
  title = "Enhancer overlap", col_fun = col_enhancer, 
  at = c(0, 0.25, 0.5, 0.75, 1), 
  labels = c("0%", "25%", "50%", "75%", "100%")
)

# 创建 png 图形设备，并设置足够的大小
# 注意：如果图形设备的大小太小，会提示 "figure margins too large"
# 并且，gridOMI() 会返回负值
png(filename = "~/Downloads/a.png", width = 1000, height = 800)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region

pushViewport(viewport(
  x = 0, y = 0.5, width = circle_size, 
  height = circle_size, just = c("left", "center"))
)
# 设置 new = TRUE，避免重新创建图形
par(omi = gridOMI(), new = TRUE)
circlize_plot()
upViewport()

# 获取图形设备的高度
h <- dev.size()[2]
lgd_list <- packLegend(
  lgd_meth, lgd_direction, lgd_expr, 
  lgd_pvalue, lgd_gene_type,  lgd_anno_gene, 
  lgd_dist, lgd_enhancer, 
  max_height = unit(0.9*h, "inch")
)
draw(lgd_list, x = circle_size, just = "left")
dev.off()
circos.clear()
