library(tidyverse)
library(UpSetR)

movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")

# 自定义散点图
my_scatter <- function(data, x, y) {
  p <- ggplot(data, aes_string(x, y, colour = "color")) +
    geom_point() +
    scale_colour_identity() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  p
}
# 自定义密度图
my_density <- function(data, x, y) {
  data$decades <- data[, y] %/% 10 * 10
  data <- data[which(data$decades >= 1970), ]
  p <- ggplot(data, aes_string(x)) +
    geom_density(aes(fill = factor(decades)), alpha = 0.3) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"), 
      legend.key.size = unit(0.4, "cm")
    )
  p
}

# 添加元数据，评分
sets <- names(movies[3:19])
avgRottenTomatoesScore <- round(runif(17, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets, avgRottenTomatoesScore))
names(metadata) <- c("sets", "avgRottenTomatoesScore")

metadata$avgRottenTomatoesScore <- as.numeric(as.character(metadata$avgRottenTomatoesScore))

# 城市信息
Cities <- sample(c("Boston", "NYC", "LA"), 17, replace = T)
metadata <- cbind(metadata, Cities)
metadata$Cities <- as.character(metadata$Cities)

# 接收状态
accepted <- round(runif(17, min = 0, max = 1))
metadata <- cbind(metadata, accepted)

upset(movies, 
      # 查询
      queries = list(
        list(
          query = intersects, 
          params = list("Drama"), 
          color = "red", 
          active = F), 
        list(
          query = intersects, 
          params = list("Action", "Drama"), 
          active = T), 
        list(
          query = intersects,
          params = list("Drama", "Comedy", "Action"), 
          color = "orange", 
          active = T)), 
      # 元数据图
      set.metadata = list(
        data = metadata, 
        plots = list(
          list(
            type = "hist", 
            column = "avgRottenTomatoesScore", 
            assign = 20), 
          list(
            type = "bool", 
            column = "accepted",
            assign = 5, 
            colors = c("#FF3333", "#006400")), 
          list(
            type = "text", 
            column = "Cities",
            assign = 5, 
            colors = c(
              Boston = "green", 
              NYC = "navy", 
              LA = "purple")), 
          list(
            type = "matrix_rows", 
            column = "Cities", 
            colors = c(
              Boston = "green", 
              NYC = "navy", 
              LA = "purple"), 
            alpha = 0.5)
          )
        ), 
      # 属性图
      attribute.plots = list(
        gridrows = 45, 
        plots = list(
          list(
            plot = my_scatter, 
            x = "ReleaseDate", 
            y = "AvgRating", 
            queries = T), 
          list(plot = my_density, 
               x = "AvgRating", 
               y = "ReleaseDate", 
               queries = F)), 
        ncols = 2), 
      query.legend = "bottom"
      )
