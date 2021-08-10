library(org.Hs.eg.db)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(gelnet)
library(tidyverse)
library(synapser)

# 1. RNA expression
# ================================================== #
#                 data processing                    #
# ================================================== #
synLogin(email = "dxsbiocc@gmail.com", password = "dxs920466915")
synRNA <- synGet("syn2701943", downloadLocation = "~/Downloads/PCBC")

exp <- read_delim(file = synRNA$path) %>%
  # cut the Ensembl suffix
  separate(col = "tracking_id", sep = "\\.", into = c("Ensembl_ID", "suffix")) %>%
  dplyr::select(-suffix) %>%
  column_to_rownames("Ensembl_ID") %>%
  as.matrix()

# ENSEMBL -> Symbol
unimap <- mapIds(
  org.Hs.eg.db, keys = rownames(exp), keytype = "ENSEMBL", 
  column = "SYMBOL", multiVals = "filter"
)

data.exp <- exp[names(unimap),]
rownames(data.exp) <- unimap

# label
synMeta <- synTableQuery("SELECT UID, Diffname_short FROM syn3156503")

metaInfo <- synMeta$asDataFrame() %>%
  dplyr::select(UID, Diffname_short) %>%
  column_to_rownames("UID") %>%
  filter(!is.na(Diffname_short))

y <- metaInfo[colnames(X), ]
names(y) <- colnames(X)

# ================================================== #
#                 training model                     #
# ================================================== #
# Mean-center the data
X <- data.exp
m <- apply(X, 1, mean)
X <- X - m

sc <- which(y == "SC")
X.sc <- X[, sc]
X.or <- X[, -sc]

model.RNA <- gelnet(t(X.sc), NULL, 0, 1)
save(X, y, model.RNA, file = "~/Downloads/PCBC/model.rda")

# ================================================== #
#                      predict                       #
# ================================================== #
load('~/Downloads/PCBC/model.rda')

get_expression <- function(proj, n = 20) {
  query <- GDCquery(
    project = proj,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "HTSeq - FPKM"
  )
  query$results[[1]] <- query$results[[1]][1:n,]
  GDCdownload(query)
  data <- GDCprepare(query)
  exp <- assay(data) %>% as.data.frame() %>%
    rownames_to_column(var = "Ensembl_ID") %>% {
      unimap <- mapIds(
        org.Hs.eg.db, keys = .$Ensembl_ID, keytype = "ENSEMBL", 
        column = "SYMBOL", multiVals = "filter")
      filter(., Ensembl_ID %in% names(unimap))
    } %>%
    mutate(Symbol = AnnotationDbi::select(
      org.Hs.eg.db, keys = .$Ensembl_ID, keytype = "ENSEMBL",
      columns = "SYMBOL")$SYMBOL, .before = Ensembl_ID) %>%
    dplyr::select(!Ensembl_ID) %>%
    filter(!is.na(Symbol)) %>%
    group_by(Symbol) %>%
    summarise_all(mean) %>% 
    column_to_rownames(var = "Symbol")
  return(exp)
}

exp <- get_expression("TCGA-BRCA")
# 
predict.mRNAsi <- function(exp, modelPath='model.rda') {
  load(modelPath)
  
  common <- intersect(names(model.RNA$w), rownames(exp))
  X <- exp[common, ]
  w <- model.RNA$w[common]
  
  score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
  score <- score - min(score)
  score <- score / max(score)
}
score <- predict.mRNAsi(exp, '~/Downloads/PCBC/model.rda')
# 2. DNA methylation
# ================================================== #
#                      load data                     #
# ================================================== #
load('~/Downloads/PCBC/pcbc.data.Rda')
load('~/Downloads/PCBC/pcbc.pd.f.Rda')

m <- apply(pcbc.data, 1, mean)
pcbc.data.norm <- pcbc.data - m

SC <- pcbc.pd.f[pcbc.pd.f$Diffname_short %in% 'SC',]
X <- pcbc.data.norm[, SC$UID]
model.DNA <- gelnet(t(X), NULL, 0, 1)

save(model.DNA, model.RNA, file = "~/Downloads/PCBC/model-weight.rda")

# ================================================== #
#                      predict                       #
# ================================================== #
coad.samples <- matchedMetExp("TCGA-COAD", n = 10)
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  legacy = TRUE,
  barcode = coad.samples
)
GDCdownload(query)
met <- GDCprepare(query, save = FALSE)

### replace NA
replace.NA <-function(data, by = "mean"){
  # Do we have NAs?
  if(is.na(table(is.na(data))["TRUE"])){
    message("No NAs were found")
    return(data)
  }
  # get NAs index 
  idx <- which(is.na(data) == TRUE,arr.ind=TRUE)
  count <- table(rownames(idx))
  message("======= Status Number of NA in probes ========")
  message("--------------------- Summary------------------")
  print(summary(as.numeric(count)))
  message("\n----------- Probes with more nb of NAs -----------")
  print(head(sort(count,decreasing = T)))
  message("===============================================")
  
  idx <- cbind(idx, mean = NA, median = NA)
  
  # For each NA value calculate the mean for the same probe for the samples
  # where it belongs
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    probe <- rownames(idx)[line]
    
    # get the probe value for all samples in the same group 
    aux <- data[rownames(data) %in% probe, ] 
    
    idx[line,3] <- mean(as.numeric(aux),na.rm = TRUE)
    idx[line,4] <- median(as.numeric(aux),na.rm = TRUE)
  }
  # Step 2 replace
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    if(by == "mean"){
      data[idx[line,1],idx[line,2]] <- idx[line,3]  
    } else if(by == "median") { 
      data[idx[line,1],idx[line,2]] <- idx[line,4]
    }
  }
  return(data)
}

data.met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))

predict.mDNAsi <- function(met, modelPath='model.rda') {
  load(modelPath)
  
  common <- intersect(names(model.DNA$w), rownames(met))
  X <- met[common, ]
  w <- model.DNA$w[common]
  
  score <- t(w) %*% X
  score <- score - min(score)
  score <- score / max(score)
}

score.meth <- predict.mDNAsi(data.met, "~/Downloads/PCBC/model-weight.rda")

# ================================================== #
#                     all models                     #
# ================================================== #
library(readxl)

models.file <- "~/Downloads/DNAmethylation_and_RNAexpression_Stemness_Signatures.xlsx"

sheets <- excel_sheets(models.file)
for (sheet in sheets) {
  n <- 0
  if (sheet != "mRNAsi")
    n <- 1
  assign(sheet, read_excel(models.file, sheet = sheet, skip = n))
}

save(list = sheets, file = "~/Downloads/Stemness_index.rda")

load("~/Downloads/Stemness_index.rda")
