## Regime B (Table 2): designs in which the calibration set is informative about the
## targets, with the gain available under the true covariance.
source("Simulation_codes/Utils.R")
source("Simulation_codes/Designs.R")
NREP <- 500
methods <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
designs <- c("halo10", "halo20", "halo40", "transect", "boundary", "expansion")

pr <- per_rep(run_designs(paste0("B_", designs), methods, NREP))
gain <- aggregate(pr$gain[pr$method == "w_fc"], pr[pr$method == "w_fc", "design", drop = FALSE], mean)
names(gain)[2] <- "available_gain"

result <- merge(summarise_reps(pr), gain)
write.csv(result, "Simulation_results/RegimeB.csv", row.names = FALSE)
