## Table 5: the regional diagnostic under whole-state holdout. For each region and split,
## krel_hat = xi * RMS(e_diag) / q_hat is estimated from the diagnostic half and the
## coverage is measured on the evaluation half; both are averaged over the 20 splits
## within each region and then over the regions of each tercile with equal weight.
d  <- read.csv("Real_data_results/diagnostic.csv")
xi <- qnorm(0.95); G <- function(t) 2 * pnorm(t) - 1
d$krel <- xi * d$kappa_diag / d$q; d$estimated <- G(xi / d$krel)
r <- aggregate(cbind(krel, estimated, coverage) ~ state, d, mean)
r <- r[order(r$krel), ]
cat(sprintf("%d regions; Spearman correlation of krel_hat with realized coverage %.2f\n",
            nrow(r), cor(r$krel, r$coverage, method = "spearman")))
nr <- nrow(r); size <- c(nr - 2 * (nr %/% 3), nr %/% 3, nr %/% 3)
r$tercile <- factor(rep(c("Lowest", "Middle", "Highest"), size), c("Lowest", "Middle", "Highest"))
tt <- aggregate(cbind(krel, estimated, coverage) ~ tercile, r, mean)
tt$regions <- as.vector(table(r$tercile))
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule",
  "Tercile of $\\widehat\\krel$ & Regions & Mean $\\widehat\\krel$ & Estimated coverage & Realized coverage \\\\",
  "\\midrule",
  sprintf("%s & $%d$ & $%.2f$ & $%.3f$ & $%.3f$ \\\\", tt$tercile, tt$regions, tt$krel,
          tt$estimated, tt$coverage),
  "\\bottomrule", "\\end{tabular}"), "Tables/Table5.tex")
