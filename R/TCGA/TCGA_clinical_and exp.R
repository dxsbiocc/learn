library(SummarizedExperiment)
library(TCGAbiolinks)
library(org.Hs.eg.db)
library(tidyverse)

#### ============================================================ ####
####     Step 1. Querying and processing clinical informations    ####
#### ============================================================ ####

get_clinical <- function(proj) {
  query <- GDCquery(
    project = proj,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
  )
  GDCdownload(query)
  clinical <- GDCprepare(query)
  return(clinical)
}

clinical.brca <- get_clinical("TCGA-BRCA")
clinical.patient <- clinical.brca$clinical_patient_brca

# define the columns which you want to use
select_cols <- c("bcr_patient_barcode", "age_at_diagnosis", "menopause_status",
             "ajcc_tumor_pathologic_pt", "ajcc_nodes_pathologic_pn", 
             "ajcc_metastasis_pathologic_pm", "ajcc_pathologic_tumor_stage",
             "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc")
# set abbreviation
my_cols <- c("barcode", "age", "menopause", "T", "N", "M", "stage", "ER", "PR", "HER2")
# NULL values
filter_row <- c("[Unknown]", "[Not Available]", "[Not Evaluated]")

# processing
my_info <- clinical.patient %>%
  dplyr::select(select_cols) %>%
  `names<-`(my_cols) %>%
  # filter rows which contains NULL value
  filter_all(all_vars(! . %in% filter_row)) %>%
  # filter rows which status is indeterminate
  filter_at(vars(ER, PR, HER2), all_vars(. %in% c("Negative", "Positive"))) %>%
  filter_at(vars(T, N, M), all_vars(!endsWith(., "X"))) %>%
  filter(
    startsWith(menopause, c("Pre", "Post")) & startsWith(stage, "Stage")
  ) %>% 
  # Custom classification
  mutate(
    age = if_else(age < 50, "<50", "≥50"),
    menopause = if_else(startsWith(menopause, "Pre"), "No", "Yes"),
    N = gsub("(N\\d).*", "\\1", N, perl = TRUE),
    T = if_else(gsub("(N\\d).*", "\\1", T, perl = TRUE) == "T1", "T1", "T2+T3+T4"),
    stage = if_else(startsWith(stage, c("Stage I", "Stage II")), "I-II", "III-IV")
  ) %>% 
  # according to ER, PR and HER2 status to classify subtype
  mutate(
    HR = if_else(ER == "Positive" | PR == "Positive", "HR+", "HR-"),
    Subtype = if_else(HER2 == "Positive", paste(HR, "HER2+", sep = "/"), paste(HR, "HER2-", sep = "/")),
    Subtype = if_else(Subtype == "HR-/HER2-", "TNBC", Subtype)
  )

#### ============================================================ ####
####        Step 2. Obtain and process expression profiles        ####
#### ============================================================ ####

get_expression <- function(proj) {
  query <- GDCquery(
    project = proj,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "HTSeq - FPKM"
  )
  
  GDCdownload(query)
  data <- GDCprepare(query)
  
  exp <- assay(data) %>% as.data.frame() %>%
    rownames_to_column(var = "Ensembl_ID") %>% {
      # if one Ensembl ID mapping to many gene symbols, remove it
      unimap <- mapIds(
        org.Hs.eg.db, keys = .$Ensembl_ID, keytype = "ENSEMBL", 
        column = "SYMBOL", multiVals = "filter")
      filter(., Ensembl_ID %in% names(unimap))
    } %>%
    # convert Ensembl to gene symbol
    mutate(Symbol = AnnotationDbi::select(
      org.Hs.eg.db, keys = .$Ensembl_ID, keytype = "ENSEMBL",
      columns = "SYMBOL")$SYMBOL, .before = Ensembl_ID) %>% 
    # remove 'Ensembl_ID' column
    dplyr::select(!Ensembl_ID) %>%
    filter(!is.na(Symbol)) %>%
    # Take the average as expression values, if a gene has multiple rows
    group_by(Symbol) %>%
    summarise_all(mean) %>% 
    column_to_rownames(var = "Symbol")
  return(exp)
}
  
exp <- get_expression("TCGA-BRCA")
# Tumor and Normal
label <- colnames(exp) %>%
  regexPipes::gsub("\\w+-\\w+-\\w+-(\\d+).*", "\\1") %>%
  as.numeric()
label <- if_else(label < 10, 1, 0)
# t-test
t_res <- t.test(exp["HAUS5", label == 1], exp["HAUS5", label == 0], alternative = "two.sided")

t_res$statistic
# 8.966298 
t_res$p.value
# 1.480187e-16

#### ============================================================ ####
#### Step 3. ANOVA for HAUS5 expression and clinical information  ####
#### ============================================================ ####

exp.haus5 <- exp[, label == 1] %>%
  `names<-`(gsub("(\\w+-\\w+-\\w+)-.*", "\\1", colnames(.))) %>%
  dplyr::select(any_of(my_info$barcode)) %>%
  filter(rownames(.) == "HAUS5") %>% t()
  
cut_off <- mean(exp.haus5)

my_table <- subset(my_info, barcode %in% rownames(exp.haus5)) %>%
  add_column(exp = exp.haus5[,"HAUS5"]) %>%
  mutate(haus5 = if_else(exp < cut_off, "Low", "High"))

table(my_table$menopause, my_table$haus5)

my_number <- my_table %>% mutate(
  age = if_else(age == "<50", 0, 1),
  menopause = if_else(menopause == "No", 0, 1),
  T = if_else(T == "T1", 0, 1),
  N = as.numeric(gsub("N(\\d)", "\\1", N)),
  stage = if_else(stage == "I-II", 0, 1),
  Subtype = case_when(
    Subtype == "HR+/HER2+" ~ 4,
    Subtype == "HR+/HER2-" ~ 3,
    Subtype == "HR-/HER2+" ~ 2,
    Subtype == "TNBC" ~ 1
  )
)

kt <- kruskal.test(exp ~ N, my_table)
kt$p.value

save.image(file = "~/Downloads/TCGA-BRCA_exp_and_clinical.RData")
