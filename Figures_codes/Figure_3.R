## Figure 3: the Regime B designs, one replication each, drawn with the generators of
## Simulation_codes/Designs.R. Each title gives the available gain from RegimeB.csv.
source("Simulation_codes/Designs.R")
gain <- read.csv("Simulation_results/RegimeB.csv")
gain <- setNames(gain$available_gain, gain$design)
COL <- c(fit = "#4C72B0", cal = "#DD8452", target = "#C44E52")
star <- function(x, y, r = 0.045) {
  a <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]; rr <- rep(c(r, r * 0.42), 5)
  for (i in seq_along(x)) polygon(x[i] + rr * cos(a), y[i] + rr * sin(a),
                                  col = COL[["target"]], border = "white") }
panel <- function(g, title, ylab = "") {
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = ylab,
       main = title, font.main = 1, cex.main = 2.1, cex.lab = 1.8)
  abline(h = c(0, 0.5, 1), v = c(0, 0.5, 1), col = "grey93", lwd = 0.6); box(col = "grey55")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), col = "grey55", cex.axis = 1.9, padj = 0.3)
  axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1, col = "grey55", cex.axis = 1.9)
  if (!is.null(g$extra)) g$extra()
  points(g$S[g$fit, ], pch = 16, col = COL[["fit"]], cex = 0.9)
  points(g$S[g$cal, ], pch = 17, col = COL[["cal"]], cex = 1.05)
  star(g$S[g$te, 1], g$S[g$te, 2])
}
P <- list(list("halo10", "Halo 10"), list("halo20", "Halo 20"), list("halo40", "Halo 40"),
          list("transect", "Transect"), list("boundary", "Boundary"),
          list("expansion", "Monitoring network"))
dir.create("Figures", showWarnings = FALSE)
png("Figures/Figure3.png", width = 2370, height = 1700, res = 185)
par(mfrow = c(2, 3), mar = c(3.2, 6.4, 4.2, 0.8), oma = c(4.2, 0, 0, 0), mgp = c(4.2, 1.1, 0))
for (i in seq_along(P)) {
  set.seed(20 + i); g <- design_B(P[[i]][[1]])
  if (P[[i]][[1]] == "boundary")
    g$extra <- function() rect(0.12, 0.12, 0.88, 0.88, lty = 2, border = "grey45")
  if (P[[i]][[1]] == "expansion") g$extra <- function() {
    th <- seq(0, 2 * pi, length.out = 200)
    for (c0 in list(c(0.22, 0.22), c(0.78, 0.25), c(0.5, 0.8)))
      lines(c0[1] + 0.07 * cos(th), c0[2] + 0.07 * sin(th), lty = 2, col = "grey45") }
  panel(g, sprintf("%s: gain %.1f%%", P[[i]][[2]], gain[[P[[i]][[1]]]]),
        ylab = c("Nearest observations moved to calibration", "", "",
                 "Sampling patterns", "", "")[i])
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
lg <- legend("bottom", c("Fitting", "Calibration", "Prediction site"), pch = c(16, 17, NA),
             col = COL, horiz = TRUE, bty = "n", cex = 2.2, pt.cex = c(2.4, 2.7, 1))
r <- 0.0105; a <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]; rr <- rep(c(r, r * 0.42), 5)
polygon(lg$text$x[3] - 0.02 + rr * cos(a), lg$text$y[3] + rr * sin(a) * 2370 / 1700,
        col = COL[["target"]], border = "white")
dev.off()
