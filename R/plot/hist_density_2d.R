library(cowplot)

N <- 500

df <- tibble(
  x1 = rnorm(n = N, mean = 2),
  x2 = rnorm(n = N, mean = 2),
  
  y1 = rnorm(n = N, mean = 2),
  y2 = rnorm(n = N, mean = 2)
)

top_hist <- ggplot(df, aes(y1)) +
  # geom_histogram(bins = 35, fill = "#1f78b4", colour = "black") +
  geom_density(fill = "#1f78b4", colour = "black") +
  theme_void()

right_hist <- ggplot(df, aes(y2)) +
  # geom_histogram(bins = 35, fill = "#1f78b4", colour = "black") +
  geom_density(fill = "#1f78b4", colour = "black") +
  coord_flip() +
  theme_void()

center <- ggplot(df, aes(y1, y2)) +
  # geom_hex(colour = "black") +
  # scale_fill_gradientn(colours = rainbow(10)) +
  geom_density2d(colour = "black") +
  geom_density2d_filled() +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.background=element_rect(fill="white",colour="black",size=0.25),
    axis.line=element_line(colour="black",size=0.25),
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    legend.position = "none"
  )

p1 <- plot_grid(top_hist, center, align = "v", 
          nrow = 2, rel_heights = c(1, 4))

p2 <- plot_grid(NULL, right_hist, align = "v",
          nrow = 2, rel_heights = c(1, 4))

plot_grid(p1, p2, ncol = 2, 
          rel_widths = c(4, 1))

glist <- list(top_hist, center, right_hist)


lay <- rbind(c(1, NA),
             c(2, 3))

grid.arrange(grobs = glist, layout_matrix = lay,
             widths=c(4,1), heights=c(1,4))
