load("~/Downloads/gbm_lgg_deg.rda")

gene_list <- rownames(DEGs.exp)

##################################################
#                  go 富集分析                   #
##################################################
system.time(
  ansEA <- TCGAanalyze_EAcomplete(
    TFname="gbm Vs lgg", gene_list
  )
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = gene_list,
  nBar = 10,
  filename = "~/Downloads/go_enrichment.pdf"
)
# clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# symbol to ID
gene.id <- bitr(
  gene_list, fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) 

# GSEA

gene_info <- DEGs.edgeR %>%
  rownames_to_column(var = "SYMBOL") %>%
  filter(SYMBOL %in% gene_list) %>%
  inner_join(., gene.id[,1:2], by = "SYMBOL") %>%
  arrange(desc(logFC))

geneList <- gene_info$logFC
names(geneList) <- as.character(gene_info$ENTREZID)


# 富集分析
go <- enrichGO(
  gene = gene.id$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)

go2 <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  ont          = "ALL",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
# 网络图
ego <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(ego, showCategory = 2, foldChange = geneList)
# 设置分类的大小，可以是 pvalue 或 geneNum
p2 <-
  cnetplot(
    ego,
    showCategory = 2,
    categorySize = "geneNum",
    foldChange = geneList
  )
# 设置圆形布局
p3 <-
  cnetplot(
    ego,
    showCategory = 3,
    foldChange = geneList,
    circular = TRUE,
    colorEdge = TRUE
  )
cowplot::plot_grid(
  p1, p2, p3, ncol = 3,
  labels = LETTERS[1:3],
  rel_widths = c(.8, .8, 1.2)
)

p1 <- heatplot(ego)
p2 <- heatplot(ego, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

edo <- pairwise_termsim(kegg)
emapplot_cluster(edo, node_scale=1.5, layout="kk") 
##################################################
#                 KEGG 富集分析                  #
##################################################
kegg <- enrichKEGG(
  gene = gene.id$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

kegg2 <- gseKEGG(
  geneList = geneList,
  organism     = 'hsa',
  minGSSize    = 120,
  pvalueCutoff = 0.25,
  verbose      = FALSE
)
gseaplot(kegg2, geneSetID = "hsa05033")
