## Table S3: a spatial t process with nu_t = 5 or 10, under the exponential (phi = 2)
## and Gaussian (phi = 6) correlations of Regime A.
source("Simulation_codes/Utils.R")
source("Simulation_codes/Designs.R")
NREP <- 500
methods <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")
designs <- c("G_exp2_t10", "G_exp2_t5", "G_gauss_t10", "G_gauss_t5")

pr <- per_rep(run_designs(designs, methods, NREP))
write.csv(summarise_reps(pr), "Simulation_results/S3_t_process.csv", row.names = FALSE)
