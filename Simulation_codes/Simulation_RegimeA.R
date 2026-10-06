## Regime A (Table 1): nine stationary designs under random sampling. The CVR comes
## from the fixed-design counterpart of each design.
source("Simulation_codes/Utils.R")
source("Simulation_codes/Designs.R")
NREP <- 500
methods <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")

random <- per_rep(run_designs(paste0("A_", names(EXCH)), methods, NREP))
fixed  <- per_rep(run_designs(paste0("D_", names(EXCH), "_fixedS"), methods, NREP))
fixed$design <- sub("_fixedS$", "", fixed$design)

result <- merge(summarise_reps(random), summarise_cvr(fixed)[c("design", "method", "cvr")])
write.csv(result, "Simulation_results/RegimeA.csv", row.names = FALSE)
