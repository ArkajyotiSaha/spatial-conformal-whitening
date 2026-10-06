## Table 4: coverage/mean width in the PM2.5 application. Random folds are pooled over
## the held-out monitors; whole-state results are averaged over the 49 regions with
## equal weight.
rf <- read.csv("Real_data_results/random_folds.csv")
ws <- read.csv("Real_data_results/whole_state.csv")
M  <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
cov3 <- function(x) sub("^0", "", sprintf("%.3f", x))
pooled <- t(sapply(M, function(m) { z <- rf[rf$method == m, ]
  c(mean(z$covered), mean(z$width)) }))
region <- aggregate(cbind(covered, width) ~ state + method, ws, mean)
state  <- t(sapply(M, function(m) colMeans(region[region$method == m, c("covered", "width")])))
row <- function(x) paste(sprintf("$%s/%.2f$", cov3(x[, 1]), x[, 2]), collapse = " & ")
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule",
             "Scheme & Naive & GSCP & LSCP & Whitened \\\\", "\\midrule",
             paste0("Random folds        & ", row(pooled), " \\\\"),
             paste0("Whole-state holdout & ", row(state), " \\\\"),
             "\\bottomrule", "\\end{tabular}"), "Tables/Table4.tex")
w <- pooled[, 2]
cat(sprintf("Random folds: whitened narrower than GSCP by %.0f%%, LSCP by %.0f%%, naive by %.0f%%\n",
            100 * (1 - w["w_fc"] / w["gscp_rf_oob"]), 100 * (1 - w["w_fc"] / w["lscp_rf_oob"]),
            100 * (1 - w["w_fc"] / w["naive"])))
