## Figure 4: the Regime C designs, one replication each, drawn with the generators of
## Simulation_codes/Designs.R.
source("Simulation_codes/Designs.R")
COL <- c(fit = "#4C72B0", cal = "#DD8452", target = "#C44E52")
star <- function(x, y, r = 0.045) {
  a <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]; rr <- rep(c(r, r * 0.42), 5)
  for (i in seq_along(x)) polygon(x[i] + rr * cos(a), y[i] + rr * sin(a),
                                  col = COL[["target"]], border = "white") }
circle <- function(x, y, r, ...) { th <- seq(0, 2 * pi, length.out = 200)
  lines(x + r * cos(th), y + r * sin(th), ...) }
panel <- function(g, title, ylab = "") {
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = ylab,
       main = title, font.main = 1, cex.main = 2.1, cex.lab = 1.8)
  abline(h = c(0, 0.5, 1), v = c(0, 0.5, 1), col = "grey93", lwd = 0.6); box(col = "grey55")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), col = "grey55", cex.axis = 1.9, padj = 0.3)
  axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1, col = "grey55", cex.axis = 1.9)
  if (!is.null(g$extra)) g$extra(g)
  points(g$S[g$fit, ], pch = 16, col = COL[["fit"]], cex = 0.9)
  points(g$S[g$cal, ], pch = 17, col = COL[["cal"]], cex = 1.05)
  star(g$S[g$te, 1], g$S[g$te, 2])
  if (!is.null(g$label)) g$label(g)
}
X <- list(
  quadrant = function(g) { rect(0.5, 0.5, 1, 1, col = adjustcolor("grey45", 0.07), border = NA)
                           rect(0.5, 0.5, 1, 1, lty = 2, border = "grey45") },
  clustered_cal = function(g) circle(0.3, 0.3, 0.30, lty = 2, col = "grey45"),
  patch = function(g) circle(0.13, 0.13, 0.10, lty = 2, col = "grey40"),
  block_holdout = function(g) { rect(0.5, 0.5, 1, 1, col = adjustcolor(COL[["target"]], 0.07), border = NA)
                                rect(0.5, 0.5, 1, 1, lty = 2, border = COL[["target"]]) },
  varying_range = function(g) abline(v = 0.5, lty = 2, col = "grey30"))
LAB <- list(
  ## the label goes as high in the withheld block as it can without covering a site
  block_holdout = function(g) {
    y <- Find(function(v) !any(g$S[g$te, 2] > v - 0.12 & g$S[g$te, 2] < v + 0.07),
              seq(0.965, 0.55, by = -0.005))
    rect(0.535, y - 0.065, 0.965, y + 0.02, col = adjustcolor("white", 0.85), border = NA)
    text(0.75, y, "withheld from both", adj = c(0.5, 1), col = COL[["target"]], cex = 1.5) },
  varying_range = function(g) {
    rect(c(0.13, 0.62), 0.88, c(0.37, 0.88), 0.975, col = adjustcolor("white", 0.85), border = NA)
    text(0.25, 0.965, expression(phi == 1), adj = c(0.5, 1), col = "grey30", cex = 1.9)
    text(0.75, 0.965, expression(phi == 10), adj = c(0.5, 1), col = "grey30", cex = 1.9) })
P <- list(list("disjoint", "Disjoint halves"), list("quadrant", "Quadrant"),
          list("clustered_cal", "Clustered calibration"), list("patch", "Distant patch"),
          list("block_holdout", "Block holdout"), list("varying_range", "Varying range"))
dir.create("Figures", showWarnings = FALSE)
png("Figures/Figure4.png", width = 2370, height = 1700, res = 185)
par(mfrow = c(2, 3), mar = c(3.2, 6.4, 4.2, 0.8), oma = c(4.2, 0, 0, 0), mgp = c(4.2, 1.1, 0))
for (i in seq_along(P)) {
  set.seed(10 + i); g <- design_C(P[[i]][[1]])
  g$extra <- X[[P[[i]][[1]]]]; g$label <- LAB[[P[[i]][[1]]]]
  panel(g, P[[i]][[2]], ylab = c("Roles assigned by region or distance", "", "",
                                 "Where a method loses coverage", "", "")[i])
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
lg <- legend("bottom", c("Fitting", "Calibration", "Prediction site"), pch = c(16, 17, NA),
             col = COL, horiz = TRUE, bty = "n", cex = 2.2, pt.cex = c(2.4, 2.7, 1))
r <- 0.0105; a <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]; rr <- rep(c(r, r * 0.42), 5)
polygon(lg$text$x[3] - 0.02 + rr * cos(a), lg$text$y[3] + rr * sin(a) * 2370 / 1700,
        col = COL[["target"]], border = "white")
dev.off()
