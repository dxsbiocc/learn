library(TCGAbiolinks)
library(tidyverse)
library(gaia)

# Select common CN technology available for GBM and LGG
#############################
## CNV data pre-processing ##
#############################
query.gbm.nocnv <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Copy number variation",
  legacy = TRUE,
  file.type = "nocnv_hg19.seg",
  sample.type = c("Primary Tumor")
)
# 只挑选 20 个样本
query.gbm.nocnv$results[[1]] <- query.gbm.nocnv$results[[1]][1:20,]

GDCdownload(query.gbm.nocnv)

cnvMatrix <- GDCprepare(query.gbm.nocnv)

cnvMatrix %<>%
  filter(abs(Segment_Mean) > 0.3) %>%
  mutate(label = if_else(Segment_Mean < -0.3, 0, 1)) %>%
  dplyr::select(-Segment_Mean) %>%
  `names<-`(c("Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")) %>%
  mutate(Chromosome = ifelse(Chromosome == "X", 23, Chromosome),
         Chromosome = ifelse(Chromosome == "Y", 24, Chromosome),
         Chromosome = as.integer(Chromosome)) %>%
  as.data.frame()

# 获取 common CNV
file <- "Downloads/CNV.hg19.bypos.111213.txt"
commonCNV <- readr::read_tsv(
  file, 
  progress = FALSE
) %>%
  mutate(markerID = paste(Chromosome, Start, sep = ":"))

# 获取探针信息
file <- "~/Downloads/genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt"
markersMatrix <-  readr::read_tsv(
  file, col_names = FALSE, 
  col_types = "ccn", progress = FALSE
)
# 处理探针数据
markersMatrix_filtered <- markersMatrix %>%
  `names<-`(c("Probe.Name", "Chromosome", "Start")) %>%
  mutate(Chromosome = ifelse(Chromosome == "X", 23, Chromosome),
         Chromosome = ifelse(Chromosome == "Y", 24, Chromosome),
         Chromosome = as.integer(Chromosome)) %>%
  mutate(markerID = paste(Chromosome, Start, sep = ":")) %>%
  filter(!duplicated(markerID)) %>%
  filter(!markerID %in% commonCNV$markerID) %>%
  dplyr::select(-markerID)

# 使用 gaia 识别 recurrent CNV
set.seed(200)
markers_obj <- load_markers(as.data.frame(markersMatrix_filtered))
nbsamples <- length(unique(cnvMatrix$Sample.Name))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
suppressWarnings({
  cancer <- "GBM"
  results <- runGAIA(
    cnv_obj, markers_obj,
    output_file_name = paste0("~/Downloads/GAIA_", cancer, "_flt.txt"),
    # 指定需要分析的变异类型，-1 分析所有的变异
    aberrations = -1,
    # 指定需要分析的染色体，默认为 -1，表示所有染色体
    chromosomes = 9,
    # 使用近似的方式可以加快计算速度
    approximation = TRUE,
    # 设置迭代次数
    num_iterations = 5000,
    threshold = 0.25
  )
})

###############################
## visualizing recurrent CNV ##
###############################
RecCNV <- as_tibble(results) %>%
  mutate_at(vars("q-value"), as.numeric) %>%
  mutate_at(vars(-"q-value"), as.integer) %>% {
    minval = min(.$`q-value`[which(.$`q-value` != 0)])
    mutate(., `q-value` = ifelse(`q-value` == 0, minval, `q-value`))
  } %>%
  mutate(score = -log10(`q-value`)) %>%
  as.data.frame()

threshold <- 0.3
gaiaCNVplot(RecCNV,threshold)

save(results, RecCNV, threshold, file = paste0("~/Downloads/", cancer,"_CNV_results.rda"))
