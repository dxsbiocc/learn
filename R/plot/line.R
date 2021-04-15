# 坡度图
library(ggrepel)

mpg %>% 
  group_by(year, manufacturer) %>%
  summarise(value = sum(displ)) %>%
  pivot_wider(names_from = year, values_from = value) %>%
  mutate(class = if_else((`1999` - `2008`) > 0, "#8dd3c7", "#bebada")) %>%
  ggplot() +
  geom_segment(aes(x = 1, xend = 2, y = `1999`, yend = `2008`, colour = class),
               size = .75, show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "solid", size = 1, colour = "#ff7f00") +
  geom_vline(xintercept = 2, linetype = "solid", size = 1, colour = "#1f78b4") +
  geom_point(aes(x = 1, y = `1999`), size = 3, shape = 21, fill = "green") +
  geom_point(aes(x = 2, y = `2008`), size = 3, shape = 21, fill = "red") +
  scale_colour_manual(labels = c("Up", "Down"), values = c("#8dd3c7", "#bebada")) +
  xlim(.5, 2.5) +
  
  geom_text_repel(aes(x = 1, y = `1999`, label = `1999`), 
                  hjust = "left", size = 3.5) +
  geom_text_repel(aes(x = 2, y = `2008`, label = `2008`), 
                  hjust = "right", size = 3.5) +
  geom_text(aes(y = 1.03*max(max(`1999`), max(`2008`))), label = "1999", x = 1,
            size = 5, hjust = 1.2) +
  geom_text(aes(y = 1.03*max(max(`1999`), max(`2008`))), label = "2008", x = 2,
            size = 5, hjust = -.2) +
  theme_void()
  
  # 路线图
  sample_n(mtcars, 10) %>%
  ggplot(aes(mpg, disp)) +
  geom_point(colour = "#69b3a2", na.rm = TRUE) +
  geom_segment(aes(xend = c(tail(mpg, n=-1), NA),
                   yend = c(tail(disp, n=-1), NA)),
               arrow = arrow(length=unit(0.3,"cm")),
               colour = "#69b3a2") +
  geom_text(aes(label = disp), hjust = 1.2) +
  theme_bw()
