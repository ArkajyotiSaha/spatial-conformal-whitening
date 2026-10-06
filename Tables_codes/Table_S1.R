## Table S1: LSCP with the fixed bandwidth h = 0.40 and with the selected bandwidth.
## Unbounded intervals count as covering and are left out of the mean width.
fx <- rbind(read.csv("Simulation_results/RegimeA.csv")[, 1:7],
            read.csv("Simulation_results/RegimeB.csv")[, 1:7],
            read.csv("Simulation_results/RegimeC.csv")[, 1:7])
fx <- fx[fx$method == "lscp_rf_oob", ]
se <- read.csv("Simulation_results/S1_bandwidth.csv")
D <- c(exp0.5 = "Exponential, $\\phi=0.5$", exp1 = "Exponential, $\\phi=1$",
       exp2 = "Exponential, $\\phi=2$", exp6 = "Exponential, $\\phi=6$",
       exp8 = "Exponential, $\\phi=8$", exp30 = "Exponential, $\\phi=30$",
       mat1.5 = "\\Matern{}, $\\nu=3/2$", mat2.5 = "\\Matern{}, $\\nu=5/2$", gauss = "Gaussian",
       halo10 = "Halo $10$", halo20 = "Halo $20$", halo40 = "Halo $40$",
       transect = "Transect", boundary = "Boundary", expansion = "Monitoring network",
       disjoint = "Disjoint halves", quadrant = "Quadrant",
       clustered_cal = "Clustered calibration", patch = "Distant patch",
       block_holdout = "Block holdout", varying_range = "Varying range")
cov3 <- function(x) sub("^0", "", sprintf("%.3f", x))

rows <- sapply(names(D), function(d) {
  f <- fx[fx$design == d, ]; s <- se[se$design == d, ]
  ## fixed h: coverage on the bounded intervals where some are unbounded
  ub_f <- f$n_unbounded - 5 * (f$reps - f$n); nb <- 5 * f$n - ub_f
  cf <- if (f$n_unbounded == 0) f$coverage else (f$coverage * 5 * f$n - ub_f) / nb
  ## selected h: over all targets, unbounded ones counted as covered
  sel <- if (is.na(s$n)) "unbounded" else sprintf("$%s/%.2f$",
           cov3((s$coverage * 5 * s$n + 5 * (s$reps - s$n)) / (5 * s$reps)), s$width)
  sprintf("%s & $%s/%.2f$ & %s & $%.1f$ & $%.2f$\\\\", D[[d]], cov3(cf), f$width, sel,
          100 * s$n_unbounded / (5 * s$reps), s$median_h) })
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule",
             "Design & Fixed $h=0.40$ & Selected $h$ & Selected unbounded (\\%) & Median selected $h$\\\\",
             "\\midrule", rows[1:9], "\\midrule", rows[10:15], "\\midrule", rows[16:21],
             "\\bottomrule", "\\end{tabular}"), "Tables/TableS1.tex")
f <- fx[fx$design == "patch", ]
cat(sprintf("Fixed h unbounded under distant patch at %.1f%% of targets\n",
            100 * f$n_unbounded / (5 * f$reps)))
