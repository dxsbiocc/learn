library(ComplexHeatmap)
library(TCGAbiolinks)
library(tidyverse)
library(circlize)
library(glue)

##################################################
## Genomic aberration overview - oncoprint plot ## 
##################################################
LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
GBMmut <- GDCquery_Maf(tumor = "GBM", pipelines = "mutect2")
mut <- LGGmut %>% add_row(GBMmut)
save(mut, file ="~/Downloads/mafMutect2LGGGBM.rda")

# 提取胶质瘤通路中基因的突变信息
Glioma_signaling <- as_tibble(TCGAbiolinks:::listEA_pathways) %>%
  filter(grepl("glioma", Pathway, ignore.case = TRUE)) %>%
  subset(Pathway == "Glioma Signaling")

Glioma_signaling_genes <- unlist(str_split(Glioma_signaling$Molecules, ","))[1:10]
mut <- mut[mut$Hugo_Symbol %in% Glioma_signaling_genes,]

gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")
clinical <- gbm_clin %>% add_row(lgg_clin)


TCGAvisualize_oncoprint(
  mut = mut,
  gene = Glioma_signaling_genes,
  annotation = clinical[, c("bcr_patient_barcode", "disease", "vital_status", "ethnicity")],
  annotation.position = "bottom",
  color=c("background"="#CCCCCC", "DEL"="#fb8072", "INS"="#80b1d3", "SNP"="#fdb462"),
  rm.empty.columns = FALSE
)

###############################################
## Genomic aberration overview - Circos plot ## 
###############################################

LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")

mut_type <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins"
)
mut <- filter(LGGmut, Variant_Classification %in% mut_type) %>%
  mutate(
    mut.id = glue("{Chromosome}:{Start_Position}-{End_Position}|{Reference_Allele}"),
    .before = Hugo_Symbol
  ) %>%
  filter(duplicated(mut.id)) %>%
  dplyr::select(c("Chromosome","Start_Position","End_Position","Variant_Classification","Hugo_Symbol")) %>%
  distinct_all() %>% {
    # 将 typeNames 作为全局变量，后面还要使用
    typeNames <<- unique(.$Variant_Classification)
    type = structure(1:length(typeNames), names = typeNames)
    mutate(., type = type[Variant_Classification])
  } %>%
  dplyr::select(c(1:3,6,4,5))

load("~/Downloads/GBM_CNV_results.rda")

threshold <- 0.3
cnv <- as_tibble(RecCNV) %>%
  # 值考虑小于阈值的区域，threshold = 0.3
  filter(`q-value` <= threshold) %>%
  select(c(1, 3, 4, 2)) %>%
  # 将染色体名称转换为 chr 开头的字符串
  mutate(
    Chromosome = ifelse(
      Chromosome == 23, "X", ifelse(Chromosome == 24, "Y", Chromosome)
    ),
    Chromosome = paste0("chr", Chromosome)
  ) %>%
  mutate(CNV = 1) %>%
  `names<-`(c("Chromosome","Start_position","End_position","Aberration_Kind","CNV"))

# 绘制图形
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
# 添加 CNV 结果
colors <- c("forestgreen","firebrick")
circos.genomicTrackPlotRegion(
  cnv, ylim = c(0, 1.2),
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region, value, ytop.column = 2,
      ybottom = 0, col = colors[value[[1]] + 1],
      border = "white"
    )
    cell.xlim = get.cell.meta.data("cell.xlim")
    circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
  }
)
# 添加突变结果
colors <- structure(c("blue","green","red", "gold"), names = typeNames[1:4])
circos.genomicTrackPlotRegion(
  mut, ylim = c(1.2, 4.2),
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(
      region, value, cex = 0.8, pch = 16, 
      col = colors[value[[2]]], ...)
  }
)
circos.clear()

# 添加图例
legend(
  -0.2, 0.2, bty = "n", y.intersp = 1,
  c("Amp", "Del"), pch = 15,
  col = c("firebrick", "forestgreen"),
  title = "CNVs", text.font = 1,
  cex = 0.8, title.adj = 0
)
legend(
  -0.2, 0, bty = "n", y.intersp = 1,
  names(colors), pch = 20, col = colors,
  title = "Mutations", text.font = 1,
  cex = 0.8, title.adj = 0
)

###############################################
## Chr 17 Genomic aberration  - Circos plot  ## 
###############################################
par(mar=c(.1, .1, .1, .1),cex=1.5)
circos.par(
  "start.degree" = 90, canvas.xlim = c(0, 1),
  canvas.ylim = c(0, 1), gap.degree = 270,
  cell.padding = c(0, 0, 0, 0), 
  track.margin = c(0.005, 0.005)
)
circos.initializeWithIdeogram(chromosome.index = "chr17")
circos.par(cell.padding = c(0, 0, 0, 0))
# 添加 CNV 结果
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(
  cnv, ylim = c(0, 1.2),
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region, value, ytop.column = 2,
      ybottom = 0, col = colors[value[[1]]],
      border = "white"
    )
    cell.xlim = get.cell.meta.data("cell.xlim")
    circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
  }
)

# 添加单基因突变结果
gene.mut <- mutate(mut, genes.mut = glue("{Hugo_Symbol}-{type}")) %>%
  filter(!duplicated(genes.mut)) %>% {
    n.mut <- table(.$genes.mut)
    mutate(., num = n.mut[.$genes.mut])
  } %>%
  mutate(
    genes.num = as.character(glue("{Hugo_Symbol}({num})")),
    type = type / 2
  ) %>%
  dplyr::select(-(6:8))
  
s.mutt$num.Freq <- NULL

colors <- c("blue","green","red","gold")
names(colors)  <- typeNames[1:4]
circos.genomicTrackPlotRegion(
  gene.mut, ylim = c(0.3, 2.2), track.height = 0.05,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(
      region, value, cex = 0.4, pch = 16,
      col = colors[value[[2]]], ...
    )
  }
)

circos.genomicTrackPlotRegion(
  gene.mut, ylim = c(0, 1), track.height = 0.1, bg.border = NA
)
i_track = get.cell.meta.data("track.index")

circos.genomicTrackPlotRegion(
  gene.mut, ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    circos.genomicText(
      region, value, y = 1, labels.column = 3,
      col = colors[value[[2]]], facing = "clockwise",
      adj = c(1, 0.5), posTransform = posTransform.text,
      cex = 0.4, niceFacing = TRUE
    )
  },
  track.height = 0.1,
  bg.border = NA
)

circos.genomicPosTransformLines(
  gene.mut,
  posTransform = function(region, value)
    posTransform.text( 
      region,  y = 1, labels = value[[3]],
      cex = 0.4, track.index = i_track + 1
    ),
  direction = "inside",
  track.index = i_track
)

circos.clear()

legend(
  0, 0.38, bty = "n", y.intersp = 1, c("Amp", "Del"),
  pch = 15, col = c("firebrick", "forestgreen"),
  title = "CNVs", text.font = 1, cex = 0.7, title.adj = 0
)
legend(
  0, 0.2, bty = "n", y.intersp = 1, names(colors),
  pch = 16, col = colors, title = "Mutations",
  text.font = 1, cex = 0.7, title.adj = 0
)
