# 可以从官网下载
source('Cibersort.R')

# 其中 nperm 指定的是置换的次数，QN 分位数归一化,如果是芯片设置为 T，如果是测序就设置为 F，测序数据最好是TPM
res_cibersort <- CIBERSORT('LM22.txt', 'Data.txt', perm = 10, QN = F)

write.table(res_cibersort, file="res_cibersort.txt", sep="\t", col.names=T, row.names=F, quote=F)
save(res_cibersort, file="res_cibersort.Rdata")

# 画图
library(ggplot2)
library(tidyverse)
colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

my_theme <- function(){
  theme(panel.grid = element_blank(),       # 网格线
        panel.border = element_blank(),     # 面板边框
        legend.position="right",            # legend位置
        legend.text = element_text(size=8), # legend内容大小
        legend.title = element_text(size=8),# legend标题大小
        axis.line = element_line(size=1),   # 坐标轴线
        text = element_text(family="Times"),# 文本字体
        axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
        axis.title = element_text(size=10,face="bold"),  # 轴标题
        plot.title = element_text(hjust=0.5,size=10))    # 距，设置绘图区域距离边的据类，上、右、下、左
}  

p1 <- res_cibersort[,1:22] %>% reshape2::melt() %>%
  ggplot(aes(x=Var1,y=value,fill=Var2)) +
  geom_bar(stat='identity') +
  coord_flip()  +
  scale_fill_manual(values =colour ) +
  theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()

pdf("cibersort.pdf")
p1
dev.off()
