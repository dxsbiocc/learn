# install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

# library
library('org.Hs.eg.db')

# 获取所有可用的表
columns(org.Hs.eg.db)

#  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
#  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
# [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"        
# [19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
# [25] "UNIGENE"      "UNIPROT" 


# keytype 配合 keys 使用，在 select 函数中匹配 keys 参数指定的 id
keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
#  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
# [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"        
# [19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
# [25] "UNIGENE"      "UNIPROT"   

#  keys 返回数据库或表的键
head(keys(org.Hs.eg.db))
# [1] "1"  "2"  "3"  "9"  "10" "11"
head(keys(org.Hs.eg.db, keytype = 'SYMBOL'))
# [1] "A1BG"  "A2M"   "A2MP1" "NAT1"  "NAT2"  "NATP"

# read gene symbol
symbol <- read.table(file = '~/Downloads/symbol.txt', sep = '\t', header = FALSE)
symbol <- as.character(unique(symbol$V1))

# 将 symbol 对应到 entrezid
entrezid <- select(org.Hs.eg.db, keys=symbol, columns = 'ENTREZID', keytype = 'SYMBOL')
# 'select()' returned 1:1 mapping between keys and columns
# 是否存在未匹配的 SYMBOL
no_map <- sort(as.character(entrezid[is.na(entrezid$ENTREZID),'SYMBOL']))
# 进一步查看是否是基因别名 alias
alias <- select(org.Hs.eg.db, keys=no_map, columns = c('SYMBOL', 'ENTREZID'), keytype = 'ALIAS')

# 'select()' returned 1:many mapping between keys and columns
# 
# >alias
# 
#       ALIAS  SYMBOL ENTREZID
# 1    FAM63A  MINDY1    55793
# 2   FAM129B  NIBAN2    64855
# 3    MB21D1    CGAS   115004
# 4      AIM1  CRYBG1      202
# 5      AIM1   AURKB     9212
# 6      AIM1 SLC45A2    51151
# 7    TMEM57   MACO1    55219
# 8     WISP1    CCN4     8840
# 9     PYCRL   PYCR3    65263
# 10 C16orf59   TEDC2    80178
# 11  SDCCAG3   ENTR1    10807
# 12   GATSL3 CASTOR1   652968
# 13 C11orf84 SPINDOC   144097
# 14   DOPEY2   DOP1B     9980
# 15    AIM1L  CRYBG2    55057
# 16  FAM109A  PHETA1   144717
# 17    TMEM2  CEMIP2    23670
# 18 KIAA1524   CIP2A    57650
# 19   FAM64A  PIMREG    54478
# 20     GSG2  HASPIN    83903
# 21 KIAA1468   RELCH    57614
# 22     MURC  CAVIN4   347273
# 23    H2AFX    H2AX     3014
# 24 HIST1H1T    H1-6     3010
# 25 C14orf80   TEDC1   283643


# 删除多重配对的结果
uni_alias <- mapIds(org.Hs.eg.db, keys = no_map, column = 'SYMBOL', keytype = 'ALIAS', multiVals = 'filter')
# 重新匹配到 id
alias_symbol_id <- select(org.Hs.eg.db, keys = uni_alias, columns = 'ENTREZID', keytype = 'SYMBOL')
# 合并结果
res <- rbind(entrezid[!is.na(entrezid$ENTREZID),], alias_symbol_id)
# 输出结果
write.table(res, file = '~/Downloads/symbol_id.txt', sep = '\t', row.names = FALSE)


