## Table 2: coverage/mean width, available and achieved gain (Regime B).
r <- read.csv("Simulation_results/RegimeB.csv")
D <- c(halo10 = "Halo $10$", halo20 = "Halo $20$", halo40 = "Halo $40$",
       transect = "Transect", boundary = "Boundary", expansion = "Monitoring network")
M <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
cov3 <- function(x) sub("^0", "", sprintf("%.3f", x))
bold <- function(x, on) ifelse(on, sprintf("\\mathbf{%s}", x), x)

available <- achieved <- c()
rows <- sapply(names(D), function(d) {
  z <- r[r$design == d, ]; z <- z[match(M, z$method), ]
  w <- round(z$width, 2)
  cells <- sprintf("$%s/%s$", cov3(z$coverage), bold(sprintf("%.2f", z$width), w == min(w)))
  available[d] <<- z$available_gain[1]
  achieved[d]  <<- 100 * (1 - z$width[4] / z$width[2])
  sprintf("%s & %s & $%.1f$ & $%.1f$\\\\", D[[d]], paste(cells, collapse = " & "),
          available[d], achieved[d]) })
writeLines(c("\\begin{tabular}{lcccccc}", "\\toprule",
             "Design & Naive & GSCP & LSCP & Whitened & Available gain & Achieved gain\\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "Tables/Table2.tex")
cat(sprintf("Correlation of achieved with available gain: %.2f\n", cor(available, achieved)))
