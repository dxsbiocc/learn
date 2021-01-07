# 导入模块
library(rjson)

# 读取数据
data <- fromJSON(file = 'Downloads/data.json')

# 查看数据格式
# > str(data)
# List of 2
# $ :List of 4
# ..$ Name      : chr "Mario"
# ..$ Age       : num 32
# ..$ Occupation: chr "Plumber"
# ..$ Rank      : num 3
# $ :List of 4
# ..$ Name      : chr "Peach"
# ..$ Age       : num 21
# ..$ Occupation: chr "Princess"
# ..$ Rank      : num 1
mario <- data[[1]]
mario$Age <- 45
peach <- data[[2]]
peach$Rank <- 9

# 转换为 json string
outJson <- toJSON(data)
# 保存为 new_data.json
write(outJson, file = "Downloads/new_data.json")

# 安装模块
install.packages("XML")
# 导入模块
library(XML)
# 解析 xml 文件
hsa <- xmlParse("Downloads/hsa05130.xml")
# 提取根节点
oot <- xmlRoot(hsa)
# 查看根节点名称
xmlName(root)
# [1] "pathway"
# 查看根节点的子节点数目
xmlSize(root)
# [1] 293
# 查看第一个子节点
root[[1]]
# <entry id="4" name="path:hsa04810" type="map" link="https://www.kegg.jp/dbget-bin/www_bget?hsa04810">
#   <graphics name="Regulation of actin cytoskeleton" fgcolor="#000000" bgcolor="#FFFFFF" type="roundrectangle" x="1237" y="777" width="119" height="34"/>
#   </entry>

root[[1]][[1]]  # 查看第一个子节点的第一个子节点
# <graphics name="Regulation of actin cytoskeleton" fgcolor="#000000" bgcolor="#FFFFFF" type="roundrectangle" x="1237" y="777" width="119" height="34"/> 

xmlSApply(root, xmlName)  # 根节点的所有子节点名称
xmlSApply(root[[1]], xmlAttrs)  # 子节点 1 的所有子节点属性
# graphics                          
# name    "Regulation of actin cytoskeleton"
# fgcolor "#000000"                         
# bgcolor "#FFFFFF"                         
# type    "roundrectangle"                  
# x       "1237"                            
# y       "777"                             
# width   "119"                             
# height  "34" 
xmlSApply(root, xmlSize)  # 所有子节点大小

# xpath 语法获取节点属性 id=4 的 entry
getNodeSet(root, "//entry[@id=4]")
# [[1]]
# <entry id="4" name="path:hsa04810" type="map" link="https://www.kegg.jp/dbget-bin/www_bget?hsa04810">
#   <graphics name="Regulation of actin cytoskeleton" fgcolor="#000000" bgcolor="#FFFFFF" type="roundrectangle" x="1237" y="777" width="119" height="34"/>
#   </entry> 
#   
#   attr(,"class")
# [1] "XMLNodeSet"

# 转换为 list
hsa_list <- xmlToList(root)
# 转换为 dataframe
hsa_df <- xmlToDataFrame(root)
# 保存
saveXML(root, file="hsa05130.xml",encoding="UTF-8")
