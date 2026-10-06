## The numbers describing the monitor network in Section 6.1: out-of-bag R^2 of the
## forest, the spatial share sigma^2 / (sigma^2 + tau^2) of the residual variance, the
## practical range 3 / phi, and the median number of other monitors within that range.
## The forest and the covariance are fitted to all 723 monitors, over twenty seeds.
source("Simulation_codes/Utils.R")
KM_PER_UNIT <- 4547          # the coordinates are Albers kilometres divided by 4547
dat <- read.csv("Real_data_codes/EPA_data/monitors.csv")
FEAT <- c("Elevation", "commercial", "residential", "industrial", "desert",
          "forest", "mobile", "suburban", "rural", "log_pop", "wind_ms")
y <- dat$pm25; X <- as.matrix(dat[, FEAT]); S <- as.matrix(dat[, c("V11", "V12")])
all <- seq_len(nrow(dat))
fits <- t(sapply(1:20, function(sd) {
  set.seed(sd)
  m  <- fit_mean(X, y, all)
  th <- suppressWarnings(fit_covariance(S, all, y - m$oob))
  c(r2 = 1 - m$model$prediction.error / var(y),
    share = th$sigma2 / (th$sigma2 + th$tau2), phi = th$phi) }))
med <- apply(fits, 2, median); rng <- 3 / med[["phi"]]
cat(sprintf("Monitors: %d\nOut-of-bag R^2: %.2f\nSpatial share of residual variance: %.0f%%\n",
            nrow(dat), med[["r2"]], 100 * med[["share"]]))
cat(sprintf("Practical range: %.0f km\nMedian number of monitors within the range: %.0f\n",
            rng * KM_PER_UNIT, median(rowSums(as.matrix(dist(S)) < rng) - 1)))
