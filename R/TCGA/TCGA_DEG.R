library(TCGAbiolinks)
library(tidyverse)


get_count <- function(cancer) {
  query <- GDCquery(
    project = cancer,
    data.category = "Gene expression",
    data.type = "Gene expression quantification",
    platform = "Illumina HiSeq",
    file.type  = "results",
    sample.type = c("Primary Tumor"),
    legacy = TRUE
  )
  # 选择 20 个样本
  query$results[[1]] <-  query$results[[1]][1:20,]
  GDCdownload(query)
  
  exp.count <- GDCprepare(
    query,
    summarizedExperiment = TRUE,
  )
  return(exp.count)
}

gbm.exp <- get_count("TCGA-GBM")
lgg.exp <- get_count("TCGA-LGG")

dataPrep_GBM <- TCGAanalyze_Preprocessing(
  object = gbm.exp,
  cor.cut = 0.6,
  datatype = "raw_count"
)

dataPrep_LGG <- TCGAanalyze_Preprocessing(
  object = lgg.exp,
  cor.cut = 0.6,
  datatype = "raw_count"
)
# 使用 gcContent 方法进行标准化
dataNorm <- TCGAanalyze_Normalization(
    tabDF = cbind(dataPrep_LGG, dataPrep_GBM),
    geneInfo = TCGAbiolinks::geneInfo,
    method = "gcContent"
)
# 分位数过滤
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut =  0.25
) 

dataLGG <- subset(dataFilt, select = gbm.exp$barcode)
dataGBM <- subset(dataFilt, select = lgg.exp$barcode)

## edgeR
DEGs.edgeR <- TCGAanalyze_DEA(
  mat1 = dataLGG,
  mat2 = dataGBM,
  Cond1type = "LGG",
  Cond2type = "GBM",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "exactTest"
)

DEGs.limma <- TCGAanalyze_DEA(
  mat1 = dataLGG,
  mat2 = dataGBM,
  pipeline = "limma",
  Cond1type = "LGG",
  Cond2type = "GBM",
  fdr.cut = 0.01 ,
  logFC.cut = 1
)

##################################################
#               DESeq2 差异基因分析              #
##################################################
library(DESeq2)

counts <- cbind(dataLGG, dataGBM)
coldata <- data.frame(
  row.names = colnames(counts),
  group = factor(ifelse(colnames(counts) %in% gbm.exp$barcode, "gbm", "lgg"))
)
# 构造 DESeqDataSet 数据结构
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ group
)
# 差异分析
dds <- DESeq(dds)
# 
resultsNames(dds) 
# [1] "Intercept"        "group_lgg_vs_gbm"
# 获取结果
res <- results(dds, name="group_lgg_vs_gbm")
# 筛选差异基因
DEGs.DESeq2 <- res[res$padj < 0.01 & abs(res$log2FoldChange) >= 1, ]

##################################################
#        三种分析方法识别的差异基因比较          #
##################################################
library(UpSetR)
library(RColorBrewer)

g.edgeR <- rownames(DEGs.edgeR)
g.limma <- rownames(DEGs.limma)
g.DESeq2 <- rownames(DEGs.DESeq2)

set_list <- list(
  edgeR = g.edgeR,
  limma = g.limma,
  DESeq2 = g.DESeq2
)
upset(
  fromList(set_list),
  order.by = "freq", 
  scale.intersections = "log10",
  sets.bar.color = brewer.pal(7, "Set2")[1:3],
  matrix.color = brewer.pal(4, "Set1")[2],
  main.bar.color = brewer.pal(7, "Set2"),
)

# 韦恩图
library(VennDiagram)

grid.newpage()
venn_ploy <- venn.diagram(
  x = list(
    edgeR = g.edgeR,
    limma = g.limma,
    DESeq2 = g.DESeq2
  ),
  filename = NULL,
  fill = brewer.pal(3, "Set2")
)
grid.draw(venn_ploy)

DEGs.common <- intersect(g.edgeR, intersect(g.limma, g.DESeq2))
DEGs.exp <- subset(dataFilt, rownames(dataFilt) %in% DEGs.common)
save(DEGs.exp, DEGs.edgeR, file = "~/Downloads/gbm_lgg_deg.rda")
