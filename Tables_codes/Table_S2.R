## Table S2: factorization and ordering, coverage and mean width averaged over the nine
## Regime A designs.
r <- rbind(read.csv("Simulation_results/RegimeA.csv")[, 1:7],
           read.csv("Simulation_results/S2_ordering.csv")[, 1:7])
O <- list(c("w_fc", "Dense Cholesky", "as supplied by the split"),
          c("w_fc_coordord", "Dense Cholesky", "coordinate sum"),
          c("w_fc_maximin", "Dense Cholesky", "maximin"),
          c("w_fc_randord", "Dense Cholesky", "random"),
          c("w_fc_nngp", "NNGP, $K=25$", "as supplied by the split"))
rows <- sapply(O, function(o) { z <- r[r$method == o[1], ]
  sprintf("%s & %s & $%.3f$ & $%.3f$\\\\", o[2], o[3], mean(z$coverage), mean(z$width)) })
writeLines(c("\\begin{tabular}{llcc}", "\\toprule",
             "Factor & Ordering & Coverage & Mean width\\\\", "\\midrule", rows,
             "\\bottomrule", "\\end{tabular}"), "Tables/TableS2.tex")
