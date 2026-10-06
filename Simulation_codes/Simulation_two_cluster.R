## The two-cluster example (Figures 1 and 2). One hundred observations in two clusters,
## the target inside the right one, an exponential covariance with phi = 6 and nugget
## 0.1, known mean and covariance, and 20,000 draws of the field. The three designs
## differ only in which observations are used for fitting and which for calibration.
## The locations are fixed: 50 uniform on each of the discs of radius 0.16 centred at
## (0.25, 0.5) and (0.75, 0.5), with the target at (0.75, 0.5) in the last row.
PHI <- 6; TAU2 <- 0.1; ALPHA <- 0.10; NREP <- 20000
H <- c(0.05, 0.20, 0.40, 1, 2)                       # LSCP bandwidths
S <- as.matrix(read.csv("Simulation_codes/two_cluster_locations.csv")); n <- 100
SITE <- S[n + 1, ]
set.seed(11)
C <- exp(-PHI * as.matrix(dist(S))) + TAU2 * diag(n + 1)
r_all <- matrix(rnorm(NREP * (n + 1)), NREP) %*% chol(C)      # one set of draws
r0 <- r_all[, n + 1]

krige <- function(idx, tgt) {
  a <- solve(C[idx, idx], C[idx, tgt, drop = FALSE])
  list(m = r_all[, idx] %*% a, v = C[1, 1] - colSums(C[idx, tgt, drop = FALSE] * a))
}
kth <- function(A, k) apply(A, 1, function(z) sort.int(z, partial = k)[k])

## the CVR benchmark: coefficient of variation of the k-th order statistic of m
## independent standard half-normal scores, by numerical integration
cv0 <- function(m) {
  k <- ceiling((1 - ALPHA) * (m + 1))
  dens <- function(x) { F <- 2 * pnorm(x) - 1
    k * choose(m, k) * F^(k - 1) * (1 - F)^(m - k) * 2 * dnorm(x) }
  m1 <- integrate(function(x) x * dens(x), 0, 10)$value
  m2 <- integrate(function(x) x^2 * dens(x), 0, 10)$value
  sqrt(m2 - m1^2) / m1
}
## coverage, mean width and CVR (population standard deviation) of a set of widths
entry <- function(method, design, covered, w, m) {
  fin <- all(is.finite(w))
  data.frame(design = design, method = method, coverage = mean(covered),
             width = if (fin) mean(w) else Inf,
             cvr = if (fin) sqrt(mean((w - mean(w))^2)) / mean(w) / cv0(m) else NA)
}

designs <- list(nearby = 51:100, mixed = c(1:25, 51:75), distant = 1:50)
table_rows <- list(); ranks <- list()
for (dn in names(designs)) {
  cal <- designs[[dn]]; fit <- setdiff(1:n, cal); m <- length(cal)
  k <- ceiling((1 - ALPHA) * (m + 1))
  ## naive split conformal
  q <- kth(abs(r_all[, cal]), k)
  table_rows[[length(table_rows) + 1]] <- entry("Naive", dn, abs(r0) <= q, 2 * q, m)
  ## GSCP: calibration and target kriged from the fitting block
  kc <- krige(fit, cal); k0 <- krige(fit, n + 1)
  sc <- abs(r_all[, cal] - kc$m) / matrix(sqrt(kc$v), NREP, m, byrow = TRUE)
  z0 <- abs(r0 - k0$m[, 1]) / sqrt(k0$v)
  qg <- kth(sc, k)
  table_rows[[length(table_rows) + 1]] <- entry("GSCP", dn, z0 <= qg, 2 * qg * sqrt(k0$v), m)
  ## LSCP: the plausibility w0 + sum of the weights of the scores at or above the target
  ## score must reach floor((m + 1) alpha)/(m + 1), the cut-off of the scp package
  d <- sqrt(colSums((t(S[cal, ]) - SITE)^2)); tM <- floor((m + 1) * ALPHA) / (m + 1)
  for (h in H) {
    ker <- exp(-d^2 / (2 * h^2)); w0 <- 1 / (1 + sum(ker)); wt <- ker / (1 + sum(ker))
    ql <- if (w0 >= tM) rep(Inf, NREP) else apply(sc, 1, function(z) {
      o <- order(z); tail <- rev(cumsum(rev(wt[o])))
      z[o][max(which(tail >= tM - w0))] })
    lab <- sprintf("LSCP, h=%s", format(h))
    table_rows[[length(table_rows) + 1]] <- entry(lab, dn, z0 <= ql, 2 * ql * sqrt(k0$v), m)
  }
  ## whitened: fitting, then calibration, then target
  ord <- c(fit, cal, n + 1); Lw <- t(chol(C[ord, ord]))
  e <- t(forwardsolve(Lw, t(r_all[, ord])))
  ew <- abs(e[, length(fit) + seq_len(m)]); zw <- abs(e[, n + 1])
  qw <- kth(ew, k)
  table_rows[[length(table_rows) + 1]] <-
    entry("Proposed (Whitened)", dn, zw <= qw, 2 * qw * Lw[n + 1, n + 1], m)
  ## rank of the target score among the calibration scores (Figure 2)
  if (dn != "mixed") ranks[[dn]] <- data.frame(design = dn,
    gscp = rowSums(sc < z0) + 1, whitened = rowSums(ew < zw) + 1)
}
dir.create("Simulation_results", showWarnings = FALSE)
write.csv(do.call(rbind, table_rows), "Simulation_results/two_cluster.csv", row.names = FALSE)
write.csv(do.call(rbind, ranks), "Simulation_results/two_cluster_ranks.csv", row.names = FALSE)
