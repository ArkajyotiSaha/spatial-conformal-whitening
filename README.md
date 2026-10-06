# spatial-conformal-whitening
Code to reproduce the simulation results and the PM2.5 data analysis in
"Conformal Prediction for Spatially Dependent Data via Sequential Whitening"
by Ayush Baran Sen and Arkajyoti Saha.

# Requirements
R with the packages `ranger`, `BRISC`, `RANN`, `readr` and `scp`
(`remotes::install_github("mhuiying/scp")`). Rebuilding the monitor data additionally
needs `dplyr`, `sf` and `terra`. Run every script from the repository root.

# Simulation
Simulation codes and results are stored in the folders Simulation_codes and
Simulation_results, respectively. Utils.R contains the methods and Designs.R the
simulated designs. The results provided are from the runs reported in the article. Each simulation
script runs 500 replications in parallel over the available cores (one on Windows) and
takes up to a few hours on a single core. Set NREP at the top of a script to a smaller
value for a quick check.

# PM2.5 data analysis
Real_data_codes/EPA_data/monitors.csv contains the annual mean PM2.5 concentrations and
covariates at the 723 monitors. Prepare_data.R builds it from the public sources listed
in the script. Running Real_data_analysis.R produces the results stored in
Real_data_results.

# Tables and figures
Table and figure codes are stored in Tables_codes and Figures_codes, and their outputs in
Tables and Figures. The table codes also print the numbers quoted in the text.

1. Figures 1 and 2: Simulation_two_cluster.R produces the simulation results. Figure_1.R produces Figure 1 and its table, and Figure_2.R produces Figure 2.
2. Table 1: Simulation_RegimeA.R produces the simulation results. Table_1.R produces Table 1.
3. Table 2 and Figure 3: Simulation_RegimeB.R produces the simulation results. Table_2.R and Figure_3.R produce Table 2 and Figure 3.
4. Table 3 and Figure 4: Simulation_RegimeC.R produces the simulation results. Table_3.R and Figure_4.R produce Table 3 and Figure 4.
5. Section 6.1: Network_summary.R produces the description of the monitor network.
6. Table 4 and Figure 5: Real_data_analysis.R produces the results. Table_4.R and Figure_5.R produce Table 4 and Figure 5.
7. Table 5: Real_data_analysis.R produces the results. Table_5.R produces Table 5.
8. Table S1: Simulation_S1.R produces the simulation results, which Table_S1.R combines with those of Regimes A to C.
9. Table S2: Simulation_S2.R produces the simulation results. Table_S2.R produces Table S2.
10. Table S3: Simulation_S3.R produces the simulation results. Table_S3.R produces Table S3.
