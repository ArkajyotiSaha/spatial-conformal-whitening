## Table S1: LSCP with the bandwidth selected by the interval score, on every design of
## Regimes A, B and C. The fixed-bandwidth entries are the LSCP results of Regimes A to C.
source("Simulation_codes/Utils.R")
source("Simulation_codes/Designs.R")
NREP <- 500
designs <- c(paste0("A_", names(EXCH)),
             paste0("B_", c("halo10", "halo20", "halo40", "transect", "boundary", "expansion")),
             paste0("C_", c("disjoint", "quadrant", "clustered_cal", "patch",
                            "block_holdout", "varying_range")))

pr <- per_rep(run_designs(designs, "lscp_rf_oob_adapt", NREP))
h  <- aggregate(pr$eta_10, pr["design"], median, na.rm = TRUE)
names(h)[2] <- "median_h"

result <- merge(summarise_reps(pr), h)
write.csv(result, "Simulation_results/S1_bandwidth.csv", row.names = FALSE)
