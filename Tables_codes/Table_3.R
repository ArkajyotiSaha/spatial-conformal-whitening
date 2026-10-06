## Table 3: coverage/mean width and predicted coverage (Regime C).
r <- read.csv("Simulation_results/RegimeC.csv")
D <- c(disjoint = "Disjoint halves", quadrant = "Quadrant",
       clustered_cal = "Clustered calibration", patch = "Distant patch",
       block_holdout = "Block holdout", varying_range = "Varying range")
M <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
cov3 <- function(x) sub("^0", "", sprintf("%.3f", x))
bold <- function(x, on) ifelse(on, sprintf("\\mathbf{%s}", x), x)

rows <- sapply(names(D), function(d) {
  z <- r[r$design == d, ]; z <- z[match(M, z$method), ]
  w <- round(z$width, 2)
  ## an interval that is unbounded at some target has infinite mean width
  cells <- ifelse(z$n_unbounded > 0, sprintf("$%s/\\infty$", cov3(z$coverage)),
                  sprintf("$%s/%s$", cov3(z$coverage),
                          bold(sprintf("%.2f", z$width), w == min(w[z$n_unbounded == 0]))))
  sprintf("%s & %s & $%s$\\\\", D[[d]], paste(cells, collapse = " & "), cov3(z$predicted[4])) })
writeLines(c("\\begin{tabular}{lccccr}", "\\toprule",
             "Design & Naive & GSCP & LSCP & Whitened & Predicted coverage\\\\", "\\midrule",
             "\\multicolumn{6}{l}{\\emph{Role-dependent sampling}}\\\\", rows[1:3], "\\midrule",
             "\\multicolumn{6}{l}{\\emph{Coverage stress tests}}\\\\", rows[4:6],
             "\\bottomrule", "\\end{tabular}"), "Tables/Table3.tex")

## numbers quoted in the text and the caption
W <- function(m, d) r$width[r$method == m & r$design == d]
up <- c("disjoint", "quadrant", "clustered_cal")
rn <- sapply(up, function(d) 100 * (1 - W("w_fc", d) / W("naive", d)))
rc <- c(sapply(up, function(d) 100 * (1 - W("w_fc", d) / W("gscp_rf_oob", d))),
        sapply(up, function(d) 100 * (1 - W("w_fc", d) / W("lscp_rf_oob", d))))
cat(sprintf("Upper row: whitened %.0f%% to %.0f%% narrower than naive, %.0f%% to %.0f%% than GSCP and LSCP\n",
            min(rn), max(rn), min(rc), max(rc)))
z <- r[r$method == "lscp_rf_oob" & r$design == "patch", ]
nb <- 5 * z$n - (z$n_unbounded - 5 * (z$reps - z$n))
cat(sprintf("Distant patch: LSCP unbounded at %.0f%% of targets; bounded intervals %s/%.2f\n",
            100 * z$n_unbounded / (5 * z$reps),
            cov3((z$coverage * 5 * z$n - (5 * z$n - nb)) / nb), z$width))
w <- r[r$method == "w_fc", ]
cat("Predicted minus realized coverage of the whitened interval:\n")
print(data.frame(design = w$design, gap = round(w$predicted - w$coverage, 3)), row.names = FALSE)
