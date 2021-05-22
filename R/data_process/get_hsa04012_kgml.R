library(rvest)
library(tidyverse)

# 解析 ErbB signaling pathway 网络结构
kmgl <- read_html('http://rest.kegg.jp/get/hsa04012/kgml')

node <- kmgl %>% html_elements(xpath = '//entry[@type="gene"]') %>%
  html_attrs() %>%
  lapply(function (x) c(x['id'], x['name'])) %>%
  do.call(rbind, .) %>%
  as.data.frame()

edge <- kmgl %>% html_elements(xpath = '//relation') %>%
  html_attrs() %>%
  lapply(function (x) c(x['entry1'], x['entry2'])) %>%
  do.call(rbind, .) %>%
  as.data.frame()

net <- inner_join(node, edge, by = c('id' = 'entry1')) %>%
  inner_join(node, by = c('entry2' = 'id')) %>%
  dplyr::select(name.x, name.y) %>%
  dplyr::rename(source = name.x, target = name.y)

edges_list <- list()
for (i in seq_along(net[,1])) {
  gids <- str_split(gsub("hsa:", '', net[i,]), ' ')
  s <- data.frame('source' = gids[[1]])
  t <- data.frame('target' = gids[[2]])
  edges_list[[i]] <- crossing(s, t)
}

all_edges <- do.call(rbind, edges_list)

# EGF-ERBB2-RAS-ERK signaling pathway genes
sub_path <- c(
  1950, 2064, 1956, 2885,
  6654, 6655, 3265, 3845,
  4893, 369, 673, 5894, 5604,
  5605, 5594, 5595
  )

path_id <- all_edges %>%
  filter(source %in% sub_path & target %in% sub_path)

# 将 ENTREZID 转换为 SYMBOL
library(org.Hs.eg.db)
s <- select(org.Hs.eg.db, keys = path_id$source, columns = 'SYMBOL', 'ENTREZID')
t <- select(org.Hs.eg.db, keys = path_id$target, columns = 'SYMBOL', 'ENTREZID')

sub_path <- tibble(
  source = s$SYMBOL,
  target = t$SYMBOL
) %>%
  distinct()

# 设置基因的类型，是否癌基因或抑癌基因
genes <- unique(union(sub_path$source, sub_path$target))
gtype <- rep("other", length(genes))
gtype[grep("(RAS)|(RAF)", genes)] <- "proto-oncogene"

gene_type <- data.frame(
  genes = genes,
  type = gtype
)
# 获取突变信息
mut <- read_delim('data_mutations_mskcc.txt', delim = '\t') %>%
  dplyr::select(Hugo_Symbol, Variant_Classification) %>%
  filter(Hugo_Symbol %in% genes)

# 每个基因的突变频数，即突变类型
mut_cnt <- group_by(mut, Hugo_Symbol) %>%
  summarise(mut_count = n())
mut_info <- distinct(mut) %>%
  group_by(Hugo_Symbol) %>%
  summarise(mut_class = n()) %>%
  left_join(mut_cnt)
# 将基因的所有信息整合在一起
node_attr <- left_join(gene_type, mut_info, by = c("genes" = "Hugo_Symbol")) %>%
  replace_na(list(mut_class = 0, mut_count = 0))
