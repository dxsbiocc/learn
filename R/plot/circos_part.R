library(circlize)

df <- data.frame(
  sectors = rep("a", 100),
  x = runif(100),
  y = runif(100)
  )
sectors <- letters[1:4]

par(mfcol = c(1, 2))
circos.par("start.degree" = 90, "gap.degree" = 0)

circos.initialize(sectors, xlim = c(0, 1))

circos.track(
  df$sectors, x = df$x, y = df$y,
  panel.fun = function(x, y) {
    circos.points(x, y, pch = 16, cex = 0.5, col = 2)
  }
)
circos.track(df$sectors, x = df$x, y = df$y,
             panel.fun = function(x, y) {
               circos.lines(1:100/100, y, col = 3)
             })
circos.clear()
rect(0, 0, 1, 1)
text(0, 0, 0, col = "red", adj = c(0.5, 1))
text(1, 0, 1, col = "red", adj = c(0.5, 1))
text(0, 1, 1, col = "red", adj = c(0.5, 0))

par(mar = c(1, 1, 1, 1))
circos.par("canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1),
           "start.degree" = 90, "gap.after" = 270)

circos.initialize(sectors = df$sectors, xlim = c(0, 1))
circos.track(
  df$sectors, x = df$x, y = df$y,
  panel.fun = function(x, y) {
    circos.points(x, y, pch = 16, cex = 0.5, col = 2)
  }
)
circos.track(
  df$sectors, x = df$x, y = df$y,
  panel.fun = function(x, y) {
    circos.lines(1:100/100, y, col = 3)
  }
)
circos.clear()
box()
par(xpd = NA)
text(0, 0, 0, col = "red", adj = c(0.5, 1))
text(1, 0, 1, col = "red", adj = c(0.5, 1))
text(0, 1, 1, col = "red", adj = c(0.5, 0))
