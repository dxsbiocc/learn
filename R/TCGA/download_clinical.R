library(TCGAbiolinks)
library(tidyverse)
library(regexPipes)
library(xlsx)

getclinical <- function(proj, filename){
  message(proj)
  while(1){
    result = tryCatch({
      query <- GDCquery(project = proj, data.category = "Clinical",file.type = "xml")
      GDCdownload(query)
      clinical <- GDCprepare_clinic(query, clinical.info = "patient")
      for(i in c("admin","radiation","follow_up","drug", "stage_event", "new_tumor_event")){
        message(i)
        aux <- GDCprepare_clinic(query, clinical.info = i)
        if(is.null(aux) || nrow(aux) == 0) next
        # 为列名添加后缀
        replicated <- which(grep("bcr_patient_barcode", colnames(aux), value = T,invert = T) %in% colnames(clinical))
        colnames(aux)[replicated] <- paste0(colnames(aux)[replicated], ".", i)
        # 合并两个表
        if(!is.null(aux)) 
          clinical <- full_join(clinical, aux, by = "bcr_patient_barcode")
      }
      write.xlsx(clinical, file = filename, sheetName = proj)
      # return(clinical)
      return(TRUE)
    }, error = function(e) {
      message(paste0("Error clinical: ", proj))
    })
    message("try again!")
  }
}

filename <- "~/Downloads/TCGA_clinical.xlsx"
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% sort %>%
  map_chr(getclinical, filename)
