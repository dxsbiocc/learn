library(tidyverse)
library(RColorBrewer)

mat1 <- as.data.frame(round(cor(mtcars), 2))

for (i in 1:11) {
  for (j in i:11) {
    mat1[i,j] <- NA
  }
}

mat2 <- as.data.frame(round(cor(mtcars), 2))

for (i in 1:11) {
  for (j in 1:i) {
    mat2[i,j] <- NA
  }
}

var_name <- data1 %>% 
  filter(var1 == var2)

mat1$var1 <- rownames(mat1)
data1 <- gather(mat1, key = "var2", value = "corr", -var1) %>%
  mutate(var1 = factor(var1, levels = rownames(mat1)),
         var2 = factor(var2, levels = rownames(mat1)))

mat2$var1 <- rownames(mat2)
data2 <- gather(mat2, key = "var2", value = "corr", -var1) %>%
  mutate(var1 = factor(var1, levels = rownames(mat2)),
         var2 = factor(var2, levels = rownames(mat2)))

my_color <- brewer.pal(5, "Spectral")

ggplot(data1, aes(var1, var2)) +
  geom_tile(data = data2, aes(fill = corr)) +
  geom_text(data = data2, aes(label = corr), colour = "black", size = 5, na.rm = TRUE) +
  geom_point(aes(fill = corr, size = corr), shape = 21, na.rm = TRUE) +
  geom_text(data = var_name, aes(label = var1), size = 5) +
  scale_fill_gradientn(colours = my_color, na.value = "white") +
  scale_colour_gradientn(colours = my_color) +
  scale_size_area(max_size = 15, guide = FALSE) +
  scale_x_discrete(position = 't') +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    axis.title = element_blank()
    )
