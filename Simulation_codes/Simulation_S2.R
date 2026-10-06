## Table S2: the whitened interval under other calibration orderings and with an NNGP
## factor, on the Regime A designs. The split order with the dense factor is in RegimeA.csv.
source("Simulation_codes/Utils.R")
source("Simulation_codes/Designs.R")
NREP <- 500
methods <- c("w_fc_coordord", "w_fc_maximin", "w_fc_randord", "w_fc_nngp")

pr <- per_rep(run_designs(paste0("A_", names(EXCH)), methods, NREP))
write.csv(summarise_reps(pr), "Simulation_results/S2_ordering.csv", row.names = FALSE)
