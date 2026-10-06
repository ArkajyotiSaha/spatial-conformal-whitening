## Figure 2: rank of the target score among the calibration scores in the two-cluster
## example, under nearby and distant calibration, for GSCP and the whitened interval.
S  <- as.matrix(read.csv("Simulation_codes/two_cluster_locations.csv"))
rk <- read.csv("Simulation_results/two_cluster_ranks.csv")
tc <- read.csv("Simulation_results/two_cluster.csv")
m <- 50; thr <- ceiling(0.9 * (m + 1)); NREP <- sum(rk$design == "nearby")
COL <- c(cal = "#DD8452", target = "#C44E52", gscp = "#4C72B0", white = "#55A868")
star <- function(x, y, r = 0.045, yscale = 1) {
  a <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]; rr <- rep(c(r, r * 0.42), 5)
  polygon(x + rr * cos(a), y + yscale * rr * sin(a), col = COL[["target"]], border = "white") }
dir.create("Figures", showWarnings = FALSE)
png("Figures/Figure2.png", width = 2400, height = 1080, res = 190)
layout(matrix(1:6, 2, byrow = TRUE), widths = c(0.85, 1.25, 1.25))
par(oma = c(3.5, 0, 0, 0), mgp = c(3, 0.9, 0))
for (d in c("nearby", "distant")) {
  cal <- if (d == "nearby") 51:100 else 1:50
  par(mar = c(3.4, 5, 3, 1), pty = "s")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "",
       ylab = paste(tools::toTitleCase(d), "calibration"), cex.lab = 1.75,
       main = if (d == "nearby") "Design" else "", font.main = 1, cex.main = 1.75)
  abline(h = c(0, 0.5, 1), v = c(0, 0.5, 1), col = "grey92", lwd = 0.6); box(col = "grey55")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), col = "grey55", cex.axis = 1.45)
  axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1, col = "grey55", cex.axis = 1.45)
  th <- seq(0, 2 * pi, length.out = 200)
  for (c0 in list(c(0.25, 0.5), c(0.75, 0.5)))
    lines(c0[1] + 0.18 * cos(th), c0[2] + 0.18 * sin(th), lty = 3, col = "grey60")
  points(S[cal, ], pch = 17, col = COL[["cal"]], cex = 1.1); star(S[101, 1], S[101, 2], r = 0.05)
  for (meth in c("gscp", "whitened")) {
    x <- rk[rk$design == d, meth]
    cv <- tc$coverage[tc$design == d & tc$method == if (meth == "gscp") "GSCP" else "Proposed (Whitened)"]
    par(mar = c(3.4, 1, 3, 1), pty = "m")
    h <- hist(x, breaks = seq(0.5, m + 2.5, 3), plot = FALSE)
    plot(NA, xlim = c(0.5, m + 1.5), ylim = c(0, 2200), axes = FALSE, xlab = "", ylab = "",
         main = sprintf("%s: coverage %.3f", if (meth == "gscp") "GSCP" else "Proposed (Whitened)", cv),
         font.main = 1, cex.main = 1.75)
    rect(thr + 0.5, 0, m + 1.5, 2200, col = adjustcolor(COL[["target"]], 0.13), border = NA)
    rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts,
         col = COL[[if (meth == "gscp") "gscp" else "white"]], border = "white")
    abline(h = NREP / ((m + 1) / 3), lty = 2, col = "grey25")
    axis(1, at = c(1, 13, 25, 38, 50), col = "grey55", cex.axis = 1.45); box(col = "grey55")
    if (d == "distant") mtext("Rank of the target score", 1, line = 2.9, cex = 1.15)
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
lg <- legend("bottomleft", c("Calibration", "Target location"), pch = c(17, NA), horiz = TRUE,
             col = COL[c("cal", "target")], bty = "n", cex = 1.5, inset = c(0.02, 0))
star(lg$text$x[2] - 0.014, lg$text$y[2], r = 0.009, yscale = 2400 / 1080)
dev.off()
