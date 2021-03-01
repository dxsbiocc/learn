library(ggplot2)
library(immunedeconv)
library(tidyverse)

# 读入表达数据
exprMatrix <- read.table("expression.txt",header=TRUE,row.names=1, as.is=TRUE)
# 使用 quantiseq 方法评估免疫浸润，method 可根据需求设置
res  <- deconvolute(exprMatrix, method="quantiseq") 
write.table(res, "quantiseq.txt", sep="\t", col.names=T, row.names=F, quote=F)
save(res, "quantiseq.Rdata")

# 画图
pdf("quantiseq.pdf")
res %>%
  gather(sample, fraction, -cell_type) %>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(res)))
dev.off()
