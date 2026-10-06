## Figure 1: the three fitting-calibration splits of the two-cluster example, and the
## table below it (Tables/Figure1_table.tex).
S <- as.matrix(read.csv("Simulation_codes/two_cluster_locations.csv"))
r <- read.csv("Simulation_results/two_cluster.csv")
COL <- c(fit = "#4C72B0", cal = "#DD8452", target = "#C44E52")
star <- function(x, y, r = 0.045, col = COL[["target"]]) {
  a <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]; rr <- rep(c(r, r * 0.42), 5)
  polygon(x + rr * cos(a), y + rr * sin(a), col = col, border = "white", lwd = 1.2) }
panel <- function(cal, title, legend = FALSE) {
  fit <- setdiff(1:100, cal)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "",
       main = title, font.main = 1, cex.main = 1.7)
  abline(h = seq(0, 1, 0.25), v = seq(0, 1, 0.25), col = "grey92", lwd = 0.6)
  box(col = "grey55")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), col = "grey55", cex.axis = 1.35)
  axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), col = "grey55", las = 1, cex.axis = 1.35)
  th <- seq(0, 2 * pi, length.out = 200)
  for (c0 in list(c(0.25, 0.5), c(0.75, 0.5)))
    lines(c0[1] + 0.18 * cos(th), c0[2] + 0.18 * sin(th), lty = 3, col = "grey60")
  points(S[fit, ], pch = 16, col = COL[["fit"]], cex = 0.85)
  points(S[cal, ], pch = 17, col = COL[["cal"]], cex = 0.9)
  star(S[101, 1], S[101, 2])
  if (legend) {
    lg <- legend("topleft", c("Fitting", "Calibration", "Target location"),
                 pch = c(16, 17, NA), col = COL, bg = "white", box.col = "grey80",
                 cex = 1.2, pt.cex = c(1, 1.1, 1))
    star(mean(c(lg$rect$left, lg$text$x[3])) - 0.008, lg$text$y[3], r = 0.034) }
}
dir.create("Figures", showWarnings = FALSE)
png("Figures/Figure1.png", width = 1500, height = 560, res = 200)
par(mfrow = c(1, 3), mar = c(2.8, 3.2, 2.8, 0.6), mgp = c(2, 0.8, 0))
panel(51:100, "(a) Nearby calibration", legend = TRUE)
panel(c(1:25, 51:75), "(b) Mixed calibration")
panel(1:50, "(c) Distant calibration")
dev.off()

## the table under the figure: coverage; mean width (CVR)
M <- c("Naive", "GSCP", paste0("LSCP, h=", c("0.05", "0.2", "0.4", "1", "2")), "Proposed (Whitened)")
lab <- c("Naive", "GSCP", paste0("LSCP, $h=", c("0.05", "0.20", "0.40", "1", "2"), "$"),
         "Proposed (Whitened)")
rows <- unlist(lapply(seq_along(M), function(i) {
  e <- sapply(c("nearby", "mixed", "distant"), function(d) {
    z <- r[r$method == M[i] & r$design == d, ]
    if (!is.finite(z$width)) sprintf("$%.3f$; $\\infty$", z$coverage)
    else sprintf("$%.3f$; $%.2f$ $(%.2f)$", z$coverage, z$width, z$cvr) })
  c(if (i == length(M)) "\\midrule", paste0(lab[i], "\n& ", paste(e, collapse = "\n& "), "\\\\")) }))
writeLines(c("\\begin{tabular}{lccc}", "\\toprule", "Method & Nearby & Mixed & Distant\\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "Tables/Figure1_table.tex")
