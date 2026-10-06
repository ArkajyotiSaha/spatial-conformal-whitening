## Regime C (Table 3): designs that assign the roles by region or break
## calibration-target comparability. Each replication carries 200 diagnostic
## observations from the target region; the predicted coverage of the whitened interval
## is the mean over replications of G(xi / krel_hat), krel_hat = xi * RMS(e_diag) / q_hat.
source("Simulation_codes/Utils.R")
source("Simulation_codes/Designs.R")
NREP <- 500
methods <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
designs <- c("disjoint", "quadrant", "clustered_cal", "patch", "block_holdout", "varying_range")

pr <- per_rep(run_designs(paste0("C_", designs), methods, NREP))
w  <- pr[pr$method == "w_fc", ]
pred <- aggregate(2 * pnorm(w$q / w$kappa_diag) - 1, w["design"], mean)
names(pred)[2] <- "predicted"; pred$method <- "w_fc"

result <- merge(summarise_reps(pr), pred, all.x = TRUE)
write.csv(result, "Simulation_results/RegimeC.csv", row.names = FALSE)
