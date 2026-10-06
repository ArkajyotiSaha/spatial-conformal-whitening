## Table S3: coverage/mean width under a spatial t process.
r <- read.csv("Simulation_results/S3_t_process.csv")
D <- list(c("exp2_t10", "Exponential, $\\phi=2$", "10"), c("exp2_t5", "Exponential, $\\phi=2$", "5"),
          c("gauss_t10", "Gaussian, $\\phi=6$", "10"), c("gauss_t5", "Gaussian, $\\phi=6$", "5"))
M <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
cov3 <- function(x) sub("^0", "", sprintf("%.3f", x))
bold <- function(x, on) ifelse(on, sprintf("\\mathbf{%s}", x), x)
rows <- sapply(D, function(d) {
  z <- r[r$design == d[1], ]; z <- z[match(M, z$method), ]; w <- round(z$width, 2)
  cells <- sprintf("$%s/%s$", cov3(z$coverage), bold(sprintf("%.2f", z$width), w == min(w)))
  sprintf("%s & $%s$ & %s\\\\", d[2], d[3], paste(cells, collapse = " & ")) })
writeLines(c("\\begin{tabular}{llcccc}", "\\toprule",
             "Correlation & $\\nu_t$ & Naive & GSCP & LSCP & Whitened\\\\", "\\midrule", rows,
             "\\bottomrule", "\\end{tabular}"), "Tables/TableS3.tex")
