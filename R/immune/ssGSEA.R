library(GSVA)
# 基因表达数据文件
expression <-"mRNA_expression.txt"
# 读取表型数据文件
pheno <- "mRNA_group.txt" 
# 背景基因集合文件
gene_set <- "gene_set.txt" 

# 读取基因表达数据
expression <- read.table(expression,sep="\t", header=T, row.names = 1)
# 读取背景基因集合
gene_set <- read.table(gene_set,header=T, sep="\t", stringsAsFactors = F)[,c(1:2)]
# 存储的是每个免疫细胞对应的基因,构建背景基因集合
bg_genes <- split(as.matrix(gene_set)[,1], gene_set[,2])

# 进行 ssGSEA 分析
gsva_matrix <- gsva(as.matrix(expression), bg_genes, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
# 输出结果
write.table(gsva_matrix,"gsva_matrix.txt", sep="\t", quote=FALSE, row.names = TRUE)

# 绘制图形
library(pheatmap)
# 样本的类型信息，为了后面画图的时候，加个注释条
pheno <- read.table(pheno,header=T,stringsAsFactors = F)
# 添加注释信息
annotation_col = data.frame(patient=pheno$type)
rownames(annotation_col)<-colnames(gsva_matrix)
# 保存图形结果
pdf(file="pheatmap_ssGSEA.pdf")
pheatmap(gsva_matrix,
         show_colnames = F,    # 不展示行名
         cluster_rows = F,     # 不对行聚类
         cluster_cols = F,     # 不对列聚类
         annotation_col = annotation_col,   # 加注释
         cellwidth=5,cellheight=5,          # 设置单元格的宽度和高度
         fontsize=5)           # 字体大小
dev.off()
