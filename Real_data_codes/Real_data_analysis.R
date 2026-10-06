## The PM2.5 application: twenty random folds and whole-state holdout for the four
## methods, and the regional diagnostic. Writes Real_data_results/.
source("Simulation_codes/Utils.R")
## read with readr, as in the reported runs (its parser differs from read.csv in the last
## binary digit of a few covariate values, which can move a forest split)
dat <- as.data.frame(readr::read_csv("Real_data_codes/EPA_data/monitors.csv", show_col_types = FALSE))
FEAT <- c("Elevation", "commercial", "residential", "industrial", "desert",
          "forest", "mobile", "suburban", "rural", "log_pop", "wind_ms")
y <- dat$pm25; X <- as.matrix(dat[, FEAT]); S <- as.matrix(dat[, c("V11", "V12")])
n <- nrow(dat)
CAL_FRAC <- 0.33333                       # a third of the non-test monitors calibrate
K <- 20                                   # random folds
methods <- c("naive", "gscp_rf_oob", "lscp_rf_oob", "w_fc")

## One method on one fold (scheme "random") or one held-out state (scheme "state").
## The split depends on the scheme and the fold but not on the method. With diag_rep > 0
## the held-out monitors are split at random into a diagnostic half, which only
## estimates the relative predictive scale, and an evaluation half.
run_one <- function(method, scheme, idx, diag_rep = 0) {
  if (scheme == "state") {
    test <- which(dat$State.Code == idx)
    set.seed(2026); rest <- sample(setdiff(seq_len(n), test))
  } else {
    set.seed(2026); fold <- sample(rep_len(seq_len(K), n))
    test <- which(fold == idx)
    set.seed(2026 + idx); rest <- sample(setdiff(seq_len(n), test))
  }
  nf <- floor((1 - CAL_FRAC) * length(rest))
  fit <- rest[seq_len(nf)]; cal <- rest[(nf + 1):length(rest)]
  yy <- y
  if (diag_rep > 0) {
    set.seed(90000 + 97 * idx + diag_rep)
    p <- sample(test); h <- length(test) %/% 2
    attr(yy, "diag_idx") <- p[seq_len(h)]; test <- p[(h + 1):length(p)]
  }
  pt <- run_method(method, yy, X, S, fit, cal, test)
  data.frame(scheme = scheme, region = idx, split = diag_rep,
             station = dat$station_id[test][pt$point], pt)
}
run_all <- function(tasks)
  rbind_fill(parallel::mclapply(tasks, function(a) do.call(run_one, a), mc.cores = N_CORES))

## random folds
folds <- run_all(do.call(c, lapply(seq_len(K), function(k)
  lapply(methods, function(m) list(m, "random", k)))))
names(folds)[names(folds) == "region"] <- "fold"
write.csv(folds[c("fold", "method", "station", "covered", "width", "dist_fit")],
          "Real_data_results/random_folds.csv", row.names = FALSE)

## whole-state holdout over the 49 contiguous state-level regions (code 80 is Mexico)
states <- sort(setdiff(unique(dat$State.Code), 80))
held <- run_all(do.call(c, lapply(states, function(s)
  lapply(methods, function(m) list(m, "state", s)))))
names(held)[names(held) == "region"] <- "state"
write.csv(held[c("state", "method", "station", "covered", "width")],
          "Real_data_results/whole_state.csv", row.names = FALSE)

## the regional diagnostic: regions with at least sixteen monitors, twenty splits each
big <- states[table(dat$State.Code)[as.character(states)] >= 16]
dg <- run_all(do.call(c, lapply(big, function(s)
  lapply(1:20, function(r) list("w_fc", "state", s, r)))))
dg <- aggregate(cbind(coverage = dg$covered, q = dg$q, kappa_diag = dg$kappa_diag),
                dg[c("region", "split")], mean)
names(dg)[1] <- "state"
write.csv(dg[order(dg$state, dg$split), ], "Real_data_results/diagnostic.csv", row.names = FALSE)
