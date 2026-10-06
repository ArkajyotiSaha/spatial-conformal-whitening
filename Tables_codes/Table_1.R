## Table 1: coverage/mean width and CVR under random sampling (Regime A).
r <- read.csv("Simulation_results/RegimeA.csv")
D <- c(exp0.5 = "Exponential, $\\phi=0.5$", exp1 = "Exponential, $\\phi=1$",
       exp2 = "Exponential, $\\phi=2$", exp6 = "Exponential, $\\phi=6$",
       exp8 = "Exponential, $\\phi=8$", exp30 = "Exponential, $\\phi=30$",
       mat1.5 = "\\Matern{}, $\\nu=3/2$", mat2.5 = "\\Matern{}, $\\nu=5/2$",
       gauss = "Gaussian")
M <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
cov3 <- function(x) sub("^0", "", sprintf("%.3f", x))
bold <- function(x, on) ifelse(on, sprintf("\\mathbf{%s}", x), x)

rows <- sapply(names(D), function(d) {
  z <- r[r$design == d, ]; z <- z[match(M, z$method), ]
  w <- round(z$width, 2); cv <- round(z$cvr, 2)
  cells <- sprintf("$%s/%s\\ (%s)$", cov3(z$coverage),
                   bold(sprintf("%.2f", z$width), w == min(w)),
                   bold(sprintf("%.2f", z$cvr), cv == min(cv)))
  paste0(D[[d]], " & ", paste(cells, collapse = " & "), "\\\\") })
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule",
             "Design & Naive & GSCP & LSCP & Whitened\\\\", "\\midrule", rows,
             "\\bottomrule", "\\end{tabular}"), "Tables/Table1.tex")

## numbers quoted in the text
W <- sapply(M, function(m) r$width[r$method == m][match(names(D), r$design[r$method == m])])
red <- function(m) 100 * (1 - W[, "w_fc"] / W[, m])
cat(sprintf("Whitened narrower than GSCP by %.0f%% and LSCP by %.0f%% on average\n",
            mean(red("gscp_rf_oob")), mean(red("lscp_rf_oob"))))
cat(sprintf("Whitened narrower than naive by %.0f%% to %.0f%%, %.0f%% on average\n",
            min(red("naive")), max(red("naive")), mean(red("naive"))))
