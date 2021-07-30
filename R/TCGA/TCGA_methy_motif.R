library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

#------------------------------------
# 获取 DNA 同时检测甲基化和表达的样本
#------------------------------------
# 结肠和直肠数据
lgg.samples <- matchedMetExp("TCGA-LGG", n = 10)
gbm.samples <- matchedMetExp("TCGA-GBM", n = 10)
samples <- c(lgg.samples,gbm.samples)

#-----------------------------------
# 1 - Methylation
# ----------------------------------
query <- GDCquery(
  project = c("TCGA-LGG","TCGA-GBM"),
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  legacy = TRUE,
  barcode = samples
)
GDCdownload(query)
met <- GDCprepare(query, save = FALSE)

# 我们以 chr9 为例
met <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
# 删除 NA 值
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
# 去除重复样本
met <- met[, substr(colnames(met), 14, 16) != "01B"]

#----------------------------
# Mean methylation
#----------------------------
TCGAvisualize_meanMethylation(
  met,
  groupCol = "project_id",
  group.legend  = "Groups",
  filename = NULL,
  print.pvalue = TRUE
)

#-------   识别差异甲基化位点   ----------
res <- TCGAanalyze_DMC(
  met,
  # colData 函数获取的矩阵中分组列名
  groupCol = "project_id",
  group1 = "TCGA-GBM",
  group2 = "TCGA-LGG",
  p.cut = 0.05,
  diffmean.cut = 0.15,
  save = FALSE,
  legend = "State",
  plot.filename = "~/Downloads/COAD_READ_metvolcano.png",
  cores =  1
)

#--------------------------
# DNA Methylation heatmap
#-------------------------
library(ComplexHeatmap)

coad_clin <- GDCquery_clinic(project = "TCGA-COAD", type = "Clinical")
read_clin <- GDCquery_clinic(project = "TCGA-READ", type = "Clinical")

use_cols <- c("bcr_patient_barcode", "disease","gender","vital_status","race")
clinical <- coad_clin %>% 
  dplyr::select(use_cols) %>%
  add_row(dplyr::select(read_clin, use_cols)) %>%
  subset(bcr_patient_barcode %in% substr(samples, 1, 12))


# 获取 Hypermethylated 和 Hypomethylated 的探针
sig_met <- filter(res, status != "Not Significant")

res_data <- subset(met,subset = (rownames(met) %in% rownames(sig_met)))

ta <- HeatmapAnnotation(
  df = clinical[, c("disease", "gender", "vital_status", "race")],
  col = list(
    disease = c("COAD" = "grey", "READ" = "black"),
    gender = c("male" = "blue", "female" = "pink")
  ))


ra <- rowAnnotation(
  df = sig_met$status,
  col = list(
    "status" =
      c("Hypomethylated" = "orange",
        "Hypermethylated" = "darkgreen")
  ),
  width = unit(1, "cm")
)

heatmap  <- Heatmap(
  assay(res_data),
  name = "DNA methylation",
  col = matlab::jet.colors(200),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  bottom_annotation = ta,
  column_title = "DNA Methylation"
) 
# Save to pdf
png("~/Downloads/heatmap.png",width = 600, height = 400)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()

save(sig_met, res_data, file = "~/Downloads/CRC.rda")
#---------------------------
#          motif 分析
#---------------------------
library(rGADEM)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifStack)

probes <- rowRanges(res_data)

sequence <- GRanges(
  seqnames = as.character(seqnames(probes)),
  ranges = IRanges(start = start(ranges(probes)) - 100,
                   end = start(ranges(probes)) + 100),
  strand = "*"
)
#look for motifs
gadem <- GADEM(sequence, verbose = FALSE, genome = Hsapiens)

nMotifs(gadem)

# 打印模体
pwm <- getPWM(gadem)
pfm  <- new("pfm",mat=pwm[[1]],name="Novel Site 1")
plotMotifLogo(pfm)

# 配对分析
library(MotIV)

analysis.jaspar <- motifMatch(pwm)

summary(analysis.jaspar)

alignment <- viewAlignments(analysis.jaspar)
print(alignment[[1]])
