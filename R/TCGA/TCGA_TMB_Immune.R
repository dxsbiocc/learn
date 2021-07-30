library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(maftools)
library(tidyverse)
library(survival)
library(survminer)
library(immunedeconv)

# -------------------------------------- #
#               1.TMB 分析               #
# -------------------------------------- #
paad.maf <- GDCquery_Maf(
  tumor = "PAAD",
  save.csv = TRUE,
  directory = "~/Downloads/TCGA",
  pipelines = "mutect2"
)

paad.clin <- GDCquery_clinic("TCGA-PAAD")
clin <- paad.clin %>% 
  dplyr::select(c(
    "submitter_id", "days_to_last_follow_up", "vital_status",
    "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
    "ajcc_pathologic_stage", "gender"
    )) %>%
  `names<-`(c("Tumor_Sample_Barcode", "time", "status", "T", "N", "M", "stage", "gender"))

paad.maf <- read.maf(
  "~/Downloads/TCGA/TCGA-PAAD/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/fea333b5-78e0-43c8-bf76-4c78dd3fac92/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz",
  clinicalData = clin,
  isTCGA = TRUE
)

getSampleSummary(paad.maf)
#Shows gene summary.
getGeneSummary(paad.maf)
#shows clinical data associated with samples
getClinicalData(paad.maf)
#Shows all fields in MAF
getFields(paad.maf)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = paad.maf, basename = 'paad.maf')

plotmafSummary(maf = paad.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
oncoplot(maf = paad.maf, colors = vc_cols, top = 15)

par(oma = c(3, 4, 5, 1))
somaticInteractions(maf = paad.maf, top = 25, pvalue = c(0.05, 0.1))

paad.maf@data[["TumorVAF"]] <- paad.maf@data$t_alt_count / paad.maf@data$t_depth
plotVaf(maf = paad.maf, vafCol = 'TumorVAF')

maf.TMB <- tcgaCompare(
  maf = paad.maf, cohortName = 'PAAD-maf', 
  logscale = TRUE, capture_size = 50
)

get_TMB <- function(file) {
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  
  df <- read_csv(file, col_select = use_cols)
  data <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data)
}

paad.csv <- get_TMB('~/Downloads/TCGA/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.csv')
paad.tmb <- paad.csv %>% 
  filter(TMB != max(TMB)) %>%
  mutate(
    label = if_else(TMB >= median(TMB), "High", "Low"), 
    .before = "MaxVAF"
  ) %>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 12))

# 生存分析
data.surv <- clin %>% filter(!is.na(time)) %>%
  mutate(
    time = time / 30, 
    status = if_else(status == "Alive", 0, 1)
  ) %>% 
  inner_join(paad.tmb, by = "Tumor_Sample_Barcode") %>%
  dplyr::select(c("Tumor_Sample_Barcode", "time", "status", "label"))

km_fit <- survfit(Surv(time, status) ~ label, data = data.surv)

ggsurvplot(
  km_fit, data = data.surv,
  pval = TRUE, surv.median.line = "hv",
  legend.labs=c("TMB-High","TMB-Low"),
  legend.title="Group",
  title="Overall survival",
  xlab = "Time(month)",
  risk.table = TRUE
)

## 其他临床特征的相关性
data.cate <- inner_join(paad.tmb, clin, by = "Tumor_Sample_Barcode") %>%
  filter_at(vars(T, N, M), all_vars(!endsWith(., "X") & !is.na(.))) %>%
  mutate(
    N = gsub("(N\\d).*", "\\1", N, perl = TRUE),
    stage = if_else(startsWith(stage, c("Stage III", "Stage IV")), "Stage III-IV", "Stage I-II")
  )
# T 分期
ggboxplot(data.cate, x = "T", y = "TMB",
          fill = "T") +
  stat_compare_means(comparisons = list(
    c("T2", "T3"), c("T2", "T4"), c("T3", "T4")
  )) +
  stat_compare_means(label.y = 4.5) 
# 性别
ggboxplot(data.cate, x = "gender", y = "TMB",
          fill = "gender") +
  stat_compare_means(comparisons = list(
    c("female", "male")
  )) +
  stat_compare_means(label.y = 4.5) 

# -------------------------------------- #
#             2. 差异表达分析            #
# -------------------------------------- #
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
  # 获取 read count
  exp.count <- GDCprepare(
    query,
    summarizedExperiment = TRUE,
  )
  return(exp.count)
}
paad.exp <- get_count("TCGA-PAAD")

# 预处理
prepData <- TCGAanalyze_Preprocessing(
  object = paad.exp,
  cor.cut = 0.6,
  datatype = "raw_count"
)

# 标准化
dataNorm <- TCGAanalyze_Normalization(
  tabDF = prepData,
  geneInfo = TCGAbiolinks::geneInfo,
)
# 分位数过滤
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut =  0.25
)

colnames(dataFilt) <- substr(colnames(dataFilt), 1, 12)
sample.high <- subset(paad.tmb, label == "High")$Tumor_Sample_Barcode
sample.low <- subset(paad.tmb, label == "Low")$Tumor_Sample_Barcode

dataHigh <- dataFilt[,colnames(dataFilt) %in% sample.high]
dataLow <- dataFilt[,colnames(dataFilt) %in% sample.low]

DEGs <- TCGAanalyze_DEA(
  mat1 = dataLow,
  mat2 = dataHigh,
  metadata = FALSE,
  Cond1type = "Low",
  Cond2type = "High",
  fdr.cut = 0.01,
  logFC.cut = 1,
  method = "glmLRT"
)

DEGsFiltLevel <- TCGAanalyze_LevelTab(
  DEGs, "TMB High", "TMB Low",
  dataHigh, dataLow
)

## ---------------------- ##
##       2.2 富集分析     ##
## ---------------------- ##
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

gene.id <- bitr(
  rownames(DEGs), fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) 
# go 
go <- enrichGO(
  gene = gene.id,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
# kegg
kegg <- enrichKEGG(
  gene = gene.id$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# gsea
gene_info <- DEGs %>%
  rownames_to_column(var = "SYMBOL") %>%
  filter(SYMBOL %in% gene_list) %>%
  inner_join(., gene.id[,1:2], by = "SYMBOL") %>%
  arrange(desc(logFC))

geneList <- gene_info$logFC
names(geneList) <- gene_info$ENTREZID

go2 <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  ont          = "ALL",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

kegg2 <- gseKEGG(
  geneList = geneList,
  organism     = 'hsa',
  minGSSize    = 120,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
## ---------------------- ##
##     2.3 免疫相关基因   ##
## ---------------------- ##
immune_gene <- read_delim("~/Downloads/GeneList.txt", delim = "\t")

gene_list <- intersect(rownames(DEGs), immune_gene$Symbol)

exp.count <- cbind(dataHigh, dataLow)

# 基因表达
log(exp.count[gene_list,]) %>% 
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "sample", values_to = "expression") %>%
  mutate(TMB = if_else(sample %in% sample.low, "Low", "High")) %>%
  ggplot(aes(Gene, expression, fill = TMB)) +
  geom_boxplot()

gene.sur <- exp.count[c("CALCA", "ERAP2"),] %>%
  t() %>% as.data.frame() %>%
  mutate(
    CALCA = if_else(CALCA < median(CALCA), 0, 1),
    ERAP2 = if_else(ERAP2 < median(ERAP2), 0, 1)
  ) %>%
  rownames_to_column("Tumor_Sample_Barcode") %>%
  inner_join(data.surv, by="Tumor_Sample_Barcode")

km.CALCA <- survfit(Surv(time, status) ~ CALCA, data = gene.sur)
km.ERAP2 <- survfit(Surv(time, status) ~ ERAP2, data = gene.sur)

ggsurvplot(
  km.CALCA, data = gene.sur,
  pval = TRUE, surv.median.line = "hv",
  legend.labs=c("High expression", "Low expression"),
  legend.title="",
  title="Overall survival",
  xlab = "Time(month)",
  risk.table = TRUE
)

# -------------------------------------- #
#             2. 免疫浸润分析            #
# -------------------------------------- #
library(GenomicFeatures)

# 计算基因长度（太麻烦了，弃）
txdb <- makeTxDbFromGFF("~/Downloads/gencode.v38.annotation.gtf.gz", format = "gtf")
exons.list <- exonsBy(txdb, by = "gene")
gene.exon.size <- lapply(exons.list, function (x) {
  sum(IRanges::width(IRanges::reduce(x)))
})

gene.length <- do.call(rbind, lapply(gene.exon.size, data.frame))

gene.length %<>%
  rownames_to_column("ENSEMBL") %>%
  mutate(ENSEMBL = substr(ENSEMBL, 1, 15))

symbol.length <- gene.length %>% {
    unimap <- mapIds(
    org.Hs.eg.db, keys = .$ENSEMBL, keytype = "ENSEMBL", 
    column = "SYMBOL", multiVals = "filter")
    filter(., ENSEMBL %in% names(unimap))
  } %>%
  mutate(Symbol = AnnotationDbi::select(
    org.Hs.eg.db, keys = .$ENSEMBL, keytype = "ENSEMBL",
    columns = "SYMBOL")$SYMBOL, .before = ENSEMBL) %>%
  dplyr::select(-ENSEMBL) %>%
  `names<-`(c("Symbol", 'lemgth'))

# 获取 FPKM
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM-UQ",
  barcode = colnames(paad.exp)
)
GDCdownload(query)
data <- GDCprepare(query)

exp.fpkm <- assay(data) %>% as.data.frame() %>%
  rownames_to_column(var = "Ensembl_ID") %>% {
    # 如果一个 Ensembl ID 映射到多个 gene symbols，要将它删除
    unimap <- mapIds(
      org.Hs.eg.db, keys = .$Ensembl_ID, keytype = "ENSEMBL", 
      column = "SYMBOL", multiVals = "filter")
    filter(., Ensembl_ID %in% names(unimap))
  } %>%
  # 将 Ensembl ID 对应到基因 symbol
  mutate(Symbol = AnnotationDbi::select(
    org.Hs.eg.db, keys = .$Ensembl_ID, keytype = "ENSEMBL",
    columns = "SYMBOL")$SYMBOL, .before = Ensembl_ID) %>% 
  # 删除 Ensembl_ID 列和未配对到 symbol 的 NA 行
  dplyr::select(!Ensembl_ID) %>%
  filter(!is.na(Symbol)) %>%
  # 对基因进行分组取平均值
  group_by(Symbol) %>%
  summarise_all(mean) %>% 
  column_to_rownames(var = "Symbol") %>%
  `names<-`(substr(colnames(.), 1, 12))

# 转换为 TPM
exp.tpm <- exp(log(exp.fpkm) - log(sum(exp.fpkm)) + log(1e6))

# TIMER 分析
immune.timer <- deconvolute(
  exp.tpm, "timer", 
  indications = rep("paad", 20)
)

immune.timer %>%
  pivot_longer(cols = -cell_type, names_to = "sample", values_to = "percent") %>%
  ggplot(aes(sample, percent, fill = cell_type)) +
  geom_bar(stat='identity') +
  coord_flip()

# cox
immune.os <- immune.timer %>%
  column_to_rownames("cell_type") %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column("Tumor_Sample_Barcode") %>%
  inner_join(data.surv, by = "Tumor_Sample_Barcode")
  
cox <- coxph(Surv(time, status) ~ `B cell` + `T cell CD4+` + `T cell CD8+` + 
               Neutrophil + Macrophage + `Myeloid dendritic cell`, 
             data = immune.os)
cox_fit <- survfit(cox)

# CIBERSORT 分析
write.table(exp.tpm, file = "~/Downloads/TCGA/PAAD_TPM.txt", sep = "\t")

source("~/Documents/WorkSpace/RStudio/Cibersort.R")

immune_cibersort <- CIBERSORT(
  sig_matrix = "~/Documents/WorkSpace/RStudio/LM22.txt", 
  mixture_file = "~/Downloads/TCGA/PAAD_TPM.txt",
  perm = 10, QN = F
)

immune_cibersort[, 1:22] %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -sample, names_to = "cell_type", values_to = "percent") %>%
  ggplot(aes(sample, percent, fill = cell_type)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1)) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

# TMB 与免疫细胞相关性
immune.tmb <- immune_cibersort[, 1:22] %>% as.data.frame() %>%
  mutate(label = if_else(rownames(immune_cibersort) %in% sample.high, 1, 0))

res.wilcox <- data.frame()
for (i in colnames(immune.tmb)[1:22]) {
  p = wilcox.test(immune.tmb[[i]], immune.tmb$label)$p.value
  res.wilcox <- rbind(res.wilcox, c(i, p))
}
colnames(res.wilcox) <- c("cell_type", "p-value")
  
immune.cate <- immune.tmb %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = c(-sample, -label), names_to = "cell_type", values_to = "percent")

label_high <- group_by(immune.cate, cell_type) %>%
  summarise(y = max(percent))

ggplot(immune.cate, aes(cell_type, percent, fill = factor(label))) +
  geom_boxplot(
    position = position_dodge2(preserve = "single", width = 0.5)) +
  stat_compare_means(aes(group = label, label = paste0("p = ", ..p.format..)), 
                     label.y = label_high$y, method = "wilcox.test") +
  scale_fill_discrete(labels = c("Low", "High"), guide_legend(title = "TMB")) +
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1))
