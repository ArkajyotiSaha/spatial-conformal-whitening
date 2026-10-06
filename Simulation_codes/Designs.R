## The simulated designs. make_rep(scen, b) returns replication b of design scen:
## the response y, covariates x, locations s, the fitting, calibration and test
## indices, the diagnostic indices (Regime C), and the available gain.
##
## Design names: A_<cov> (Regime A), D_<cov>_fixedS (fixed-design counterpart of A_<cov>,
## used for the CVR), B_<name> (Regime B), C_<name> (Regime C), G_<cov>_t<df> (spatial
## t process).

N_TOTAL <- 501; N_TEST <- 5; M_CAL <- 166; CAL_FRAC <- 1/3
N_DIAG  <- 200                       # diagnostic observations in Regime C
SIG2 <- 1.0; PHI <- 2.0; TAU2 <- 0.1

friedman <- function(Z) (10 * sin(pi * Z[,1] * Z[,2]) + 20 * (Z[,3] - 0.5)^2 +
                         10 * Z[,4] + 5 * Z[,5]) / 6
true_cov <- function(S, phi = PHI, nu = 0.5) {
  D <- as.matrix(dist(S))
  R <- if (is.infinite(nu)) exp(-(phi * D)^2)
       else if (abs(nu - 0.5) < 1e-8) exp(-phi * D)
       else if (abs(nu - 1.5) < 1e-8) { z <- sqrt(3) * phi * D; (1 + z) * exp(-z) }
       else { z <- sqrt(5) * phi * D; (1 + z + z^2 / 3) * exp(-z) }
  SIG2 * R + TAU2 * diag(nrow(S))
}
## tdf gives a multivariate t field: one inverse-gamma scale for the whole vector
sim_resid <- function(S, ..., nugget = TRUE, tdf = NA) {
  C <- true_cov(S, ...)
  if (!nugget) diag(C) <- diag(C) - TAU2
  z <- as.numeric(t(chol(C + 1e-9 * diag(nrow(S)))) %*% rnorm(nrow(S)))
  if (is.na(tdf)) return(z)
  z * sqrt(1 / rgamma(1, shape = tdf / 2, rate = tdf / 2))
}
unif_locs <- function() cbind(runif(N_TOTAL), runif(N_TOTAL))

## five prediction sites, the remainder permuted and cut at 1 - cal_frac
split_random <- function(S, cal_frac = CAL_FRAC) {
  force(S)                                    # draw the locations before the split
  test <- sample(seq_len(N_TOTAL), N_TEST)
  rest <- sample(setdiff(seq_len(N_TOTAL), test))
  nf <- floor((1 - cal_frac) * length(rest))
  list(S = S, fit = rest[seq_len(nf)], cal = rest[(nf + 1):length(rest)], te = test)
}
nn_of <- function(S, s0, k)
  order(sqrt(rowSums((S - matrix(s0, nrow(S), 2, byrow = TRUE))^2)))[seq_len(k)]

## ---- Regime B ---------------------------------------------------------------------
design_B <- function(name) {
  if (grepl("^halo", name)) {
    k <- as.integer(sub("halo", "", name))
    S <- unif_locs(); te <- sample(seq_len(N_TOTAL), N_TEST)
    near <- setdiff(unique(unlist(lapply(te, function(i) nn_of(S, S[i, ], k + 1)))), te)
    if (length(near) > M_CAL) near <- near[seq_len(M_CAL)]
    pool <- setdiff(seq_len(N_TOTAL), c(te, near))
    cal <- c(near, sample(pool, M_CAL - length(near)))
    list(S = S, fit = setdiff(seq_len(N_TOTAL), c(te, cal)), cal = cal, te = te)
  } else if (name == "expansion") {                       # monitoring network
    ctr <- rbind(c(0.22, 0.22), c(0.78, 0.25), c(0.5, 0.8))
    nf <- N_TOTAL - M_CAL - N_TEST
    Sf <- t(vapply(seq_len(nf), function(i) { c0 <- ctr[1 + (i %% 3), ]
      th <- runif(1, 0, 2 * pi); r <- 0.07 * sqrt(runif(1))
      c(c0[1] + r * cos(th), c0[2] + r * sin(th)) }, numeric(2)))
    far <- function(z) apply(z, 1, function(s)
      min(sqrt(rowSums((ctr - matrix(s, 3, 2, byrow = TRUE))^2))))
    Sr <- cbind(runif(M_CAL + N_TEST), runif(M_CAL + N_TEST))
    while (any(far(Sr) < 0.18)) { b <- far(Sr) < 0.18
      Sr[b, ] <- cbind(runif(sum(b)), runif(sum(b))) }
    S <- rbind(Sf, Sr); rest <- nf + seq_len(M_CAL + N_TEST); te <- rest[seq_len(N_TEST)]
    list(S = S, fit = seq_len(nf), cal = setdiff(rest, te), te = te)
  } else if (name == "transect") {
    nf <- N_TOTAL - M_CAL - N_TEST; ly <- seq(0.1, 0.9, length.out = 5)
    Sf <- cbind(runif(nf), ly[1 + (seq_len(nf) %% 5)] + rnorm(nf, 0, 0.012))
    mid <- ly[-5] + diff(ly)[1] / 2
    Sr <- cbind(runif(M_CAL + N_TEST),
                mid[1 + (seq_len(M_CAL + N_TEST) %% 4)] + rnorm(M_CAL + N_TEST, 0, 0.02))
    S <- rbind(Sf, Sr); rest <- nf + seq_len(M_CAL + N_TEST); te <- sample(rest, N_TEST)
    list(S = S, fit = seq_len(nf), cal = setdiff(rest, te), te = te)
  } else if (name == "boundary") {
    mg <- 0.12; nf <- N_TOTAL - M_CAL - N_TEST
    Sf <- cbind(runif(nf, mg, 1 - mg), runif(nf, mg, 1 - mg))
    Sr <- matrix(NA_real_, M_CAL + N_TEST, 2); k <- 0
    while (k < nrow(Sr)) { z <- c(runif(1), runif(1))
      if (min(z) < mg || max(z) > 1 - mg) { k <- k + 1; Sr[k, ] <- z } }
    S <- rbind(Sf, Sr); rest <- nf + seq_len(M_CAL + N_TEST); te <- sample(rest, N_TEST)
    list(S = S, fit = seq_len(nf), cal = setdiff(rest, te), te = te)
  } else stop("unknown design: ", name)
}

## ---- Regime C ---------------------------------------------------------------------
design_C <- function(name) {
  if (name == "disjoint") {
    nf <- N_TOTAL - M_CAL - N_TEST
    S <- rbind(cbind(runif(nf, 0, 0.45), runif(nf)),
               cbind(runif(M_CAL + N_TEST, 0.55, 1), runif(M_CAL + N_TEST)))
    rest <- nf + seq_len(M_CAL + N_TEST); te <- sample(rest, N_TEST)
    list(S = S, fit = seq_len(nf), cal = setdiff(rest, te), te = te)
  } else if (name == "quadrant") {
    S <- unif_locs(); q <- S[,1] > 0.5 & S[,2] > 0.5
    te <- sample(which(q), N_TEST)
    cal <- c(setdiff(which(q), te), sample(which(!q), max(0, M_CAL - sum(q) + N_TEST)))
    list(S = S, fit = setdiff(seq_len(N_TOTAL), c(te, cal)), cal = cal, te = te)
  } else if (name == "clustered_cal") {
    S <- unif_locs()
    d <- sqrt(rowSums((S - matrix(c(0.3, 0.3), N_TOTAL, 2, byrow = TRUE))^2))
    cal <- order(d)[seq_len(M_CAL)]
    te <- sample(setdiff(order(d, decreasing = TRUE)[seq_len(120)], cal), N_TEST)
    list(S = S, fit = setdiff(seq_len(N_TOTAL), c(cal, te)), cal = cal, te = te)
  } else if (name == "patch") {                           # distant patch
    th <- runif(M_CAL, 0, 2 * pi); rr <- 0.10 * sqrt(runif(M_CAL))
    Sc <- cbind(0.13 + rr * cos(th), 0.13 + rr * sin(th))
    nf <- N_TOTAL - M_CAL - N_TEST
    Sr <- cbind(runif(nf + N_TEST, 0.42, 1), runif(nf + N_TEST, 0.42, 1))
    S  <- rbind(Sc, Sr)
    rest <- M_CAL + seq_len(nf + N_TEST); te <- sample(rest, N_TEST)
    list(S = S, fit = setdiff(rest, te), cal = seq_len(M_CAL), te = te)
  } else if (name == "block_holdout") {
    S <- unif_locs(); inside <- S[,1] > 0.5 & S[,2] > 0.5
    te <- sample(which(inside), N_TEST)
    rest <- sample(which(!inside)); n_cal <- min(M_CAL, length(rest) - 20)
    list(S = S, fit = rest[-seq_len(n_cal)], cal = rest[seq_len(n_cal)], te = te)
  } else if (name == "varying_range") {
    ## two independent exponential fields, phi = 1 left and 10 right of s1 = 0.5
    S <- unif_locs(); right <- S[,1] >= 0.5
    te <- sample(which(right), N_TEST)
    rest <- sample(setdiff(seq_len(N_TOTAL), te))
    list(S = S, fit = rest[-seq_len(M_CAL)], cal = rest[seq_len(M_CAL)], te = te,
         stitch = list(left = !right, phi = c(1, 10), nu = 0.5, tau2 = 0.10))
  } else stop("unknown design: ", name)
}

## the diagnostic observations, drawn from the region the targets occupy
diag_locs <- function(name, S, te, k = N_DIAG) {
  draw <- function(f) { Z <- matrix(NA_real_, k, 2); j <- 0
    while (j < k) { z <- c(runif(1), runif(1))
      if (f(z)) { j <- j + 1; Z[j, ] <- z } }
    Z }
  if (name == "disjoint")           draw(function(z) z[1] > 0.55)
  else if (name == "quadrant")      draw(function(z) z[1] > 0.5 && z[2] > 0.5)
  else if (name == "patch")         draw(function(z) all(z > 0.42))
  else if (name == "block_holdout") draw(function(z) z[1] > 0.5 && z[2] > 0.5)
  else if (name == "varying_range") draw(function(z) z[1] >= 0.5)
  else {                                                  # clustered_cal: the far field
    d0 <- sqrt(rowSums((S[te, , drop = FALSE] -
                        matrix(c(0.3, 0.3), length(te), 2, byrow = TRUE))^2))
    draw(function(z) sqrt(sum((z - c(0.3, 0.3))^2)) >= min(d0))
  }
}

## the gain the design makes available under the true covariance, in percent
available_gain <- function(S, fit, cal, te, ...) {
  C <- true_cov(S, ...); c00 <- C[1, 1]
  nu <- function(idx, i) {
    L <- chol(C[idx, idx] + 1e-9 * diag(length(idx)))
    max(c00 - sum(backsolve(L, C[idx, i], transpose = TRUE)^2), 1e-12) }
  100 * mean(vapply(te, function(i) 1 - sqrt(nu(c(fit, cal), i) / nu(fit, i)), numeric(1)))
}

## ---- the designs -------------------------------------------------------------------
EXCH <- list(exp0.5 = c(0.5, 0.5), exp1 = c(1, 0.5), exp2 = c(2, 0.5),
             exp6 = c(6, 0.5), exp8 = c(8, 0.5), exp30 = c(30, 0.5),
             mat1.5 = c(6, 1.5), mat2.5 = c(6, 2.5), gauss = c(6, Inf))

## one set of locations and one split shared by every fixed-design counterpart
fixed_design <- local({ g <- NULL; function() {
  if (is.null(g)) { set.seed(20260824); g <<- split_random(unif_locs()) }
  g } })

design_spec <- function(scen) {
  parts <- strsplit(scen, "_")[[1]]; grp <- parts[1]
  name <- sub("^[A-Z]_", "", scen)
  if (grp == "A") { p <- EXCH[[name]]
    return(list(gen = function() split_random(unif_locs()),
                cov_par = list(phi = p[1], nu = p[2]))) }
  if (grp == "D") { p <- EXCH[[sub("_fixedS$", "", name)]]
    return(list(gen = fixed_design, fixed = TRUE, cov_par = list(phi = p[1], nu = p[2]))) }
  if (grp == "B") return(list(gen = function() design_B(name)))
  if (grp == "C") return(list(gen = function() {
    g <- design_C(name); Sd <- diag_locs(name, g$S, g$te)
    g$diag <- nrow(g$S) + seq_len(nrow(Sd)); g$S <- rbind(g$S, Sd)
    if (!is.null(g$stitch)) g$stitch$left <- c(g$stitch$left, Sd[, 1] < 0.5)
    g }))
  if (grp == "G") { p <- EXCH[[parts[2]]]
    return(list(gen = function() split_random(unif_locs()),
                cov_par = list(phi = p[1], nu = p[2]),
                tdf = as.numeric(sub("t", "", parts[3])))) }
  stop("unknown design: ", scen)
}

make_rep <- function(scen, b) {
  sp <- design_spec(scen)
  if (isTRUE(sp$fixed)) fixed_design()
  set.seed(7000 + 131 * b + sum(utf8ToInt(scen)))
  if (isTRUE(sp$fixed) && b == 1) set.seed(11)
  g <- sp$gen(); NN <- nrow(g$S)
  X <- matrix(runif(NN * 5), NN, 5)
  r <- if (!is.null(g$stitch)) {
    st <- g$stitch; z <- numeric(NN)
    for (k in 1:2) { ix <- if (k == 1) st$left else !st$left
      if (sum(ix) > 1) z[ix] <- sim_resid(g$S[ix, , drop = FALSE], phi = st$phi[k],
                                          nu = st$nu, nugget = FALSE) }
    z + sqrt(st$tau2) * rnorm(NN)
  } else do.call(sim_resid, c(list(g$S), sp$cov_par, list(tdf = if (is.null(sp$tdf)) NA else sp$tdf)))
  gain <- if (startsWith(scen, "B_")) do.call(available_gain,
            c(list(g$S, g$fit, g$cal, g$te), sp$cov_par)) else NA_real_
  list(y = friedman(X) + r, x = X, s = g$S, fit = g$fit, cal = g$cal, test = g$te,
       diag_idx = g$diag, gain = gain)
}
