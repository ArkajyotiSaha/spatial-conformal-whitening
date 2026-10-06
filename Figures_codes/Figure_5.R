## Figure 5: coverage and mean width by quintile of the distance from each held-out
## monitor to the nearest fitting monitor, under the twenty random folds.
rf <- read.csv("Real_data_results/random_folds.csv")
KM_PER_UNIT <- 4547
M <- c(naive = "Naive", gscp_rf_oob = "GSCP", lscp_rf_oob = "LSCP", w_fc = "Whitened")
br <- quantile(rf$dist_fit[rf$method == "w_fc"], seq(0, 1, 0.2))
rf$quintile <- cut(rf$dist_fit, br, include.lowest = TRUE, labels = FALSE)
cv <- tapply(rf$covered, rf[c("quintile", "method")], mean)
wd <- tapply(rf$width, rf[c("quintile", "method")], mean)
km <- tapply(rf$dist_fit[rf$method == "w_fc"] * KM_PER_UNIT,
             rf$quintile[rf$method == "w_fc"], range)
lab <- sprintf("Q%d\n%.0f-%.0f", 1:5, sapply(km, `[`, 1), sapply(km, `[`, 2))
ST <- list(naive = list(col = "#9AA0A6", lty = 2, pch = 16, lwd = 1.6),
           gscp_rf_oob = list(col = "#4C72B0", lty = 1, pch = 15, lwd = 1.9),
           lscp_rf_oob = list(col = "#55A868", lty = 1, pch = 17, lwd = 1.9),
           w_fc = list(col = "#C44E52", lty = 1, pch = 18, lwd = 2.6))
panel <- function(v, ylim, title, nominal = FALSE) {
  plot(NA, xlim = c(0.8, 5.2), ylim = ylim, axes = FALSE, xlab = "", ylab = "",
       main = title, font.main = 1, cex.main = 1.6)
  abline(h = pretty(ylim), v = 1:5, col = "grey93", lwd = 0.6); box(col = "grey55")
  if (nominal) { abline(h = 0.90, lty = 3, col = "grey25")
                 text(3.5, 0.8985, "nominal", adj = c(0.5, 1), col = "grey35", cex = 1.15) }
  for (m in names(M)) with(ST[[m]], lines(1:5, v[, m], col = col, lty = lty, lwd = lwd,
                                          type = "o", pch = pch, cex = 1.3))
  axis(1, at = 1:5, labels = lab, padj = 0.6, col = "grey55", cex.axis = 1.25)
  axis(2, las = 1, col = "grey55", cex.axis = 1.25)
  mtext("Distance to the nearest fitting monitor (km)", 1, line = 4.1, cex = 1.2)
}
dir.create("Figures", showWarnings = FALSE)
png("Figures/Figure5.png", width = 2130, height = 980, res = 190)
par(mfrow = c(1, 2), mar = c(5.6, 4.6, 3, 1), oma = c(2.8, 0, 0, 0))
panel(cv, c(0.81, 0.945), "Coverage", nominal = TRUE)
panel(wd, c(3.9, 6.35), expression("Mean width (" * mu * "g/m"^3 * ")"))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
legend("bottom", M, col = sapply(ST, `[[`, "col"), lty = sapply(ST, `[[`, "lty"),
       pch = sapply(ST, `[[`, "pch"), lwd = sapply(ST, `[[`, "lwd"), horiz = TRUE,
       bty = "n", cex = 1.35)
dev.off()

## the numbers quoted in the text
cat("Coverage by quintile:\n"); print(round(cv[, names(M)], 3))
cat("Width relative to naive:\n"); print(round(1 - wd[, names(M)] / wd[, "naive"], 3))
