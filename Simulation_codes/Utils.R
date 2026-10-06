## Methods shared by the simulation and the PM2.5 application.
##
## Every method takes (y, X, S, fit, cal, test): the response, covariates, coordinates,
## and the indices of the fitting, calibration and test observations.
##   naive              split conformal on the raw residuals
##   gscp_rf_oob        GSCP of Mao, Martin and Reich (2024), split form
##   lscp_rf_oob        LSCP with bandwidth h = 0.40
##   lscp_rf_oob_adapt  LSCP with the bandwidth selected by the interval score
##   w_fc               the whitened interval (dense Cholesky factor)
##   w_fc_nngp          the whitened interval with an NNGP factor, K = 25
##   w_fc_coordord, w_fc_maximin, w_fc_randord
##                      the whitened interval under other calibration orderings
suppressMessages({library(ranger); library(BRISC); library(RANN); library(scp)})

ALPHA_RUN <- c(0.20, 0.10, 0.05)
NTREE     <- 500
ETA_LOCAL <- 0.40      # LSCP bandwidth
ETA_VAL_N <- 100L      # validation points for the bandwidth selection
NNGP_K    <- 25L       # neighbours of the NNGP factor
N_CORES   <- if (.Platform$OS.type == "windows") 1L else parallel::detectCores()

## ---- mean and working covariance -------------------------------------------
## The covariance is fitted to the out-of-bag residuals of the forest, since in-sample
## forest residuals are artificially small.
mean_frame <- function(X, idx) {
  d <- data.frame(X[idx, , drop = FALSE]); names(d) <- paste0("x", seq_len(ncol(X))); d
}
fit_mean <- function(X, y, fit) {
  d <- mean_frame(X, fit); d$y <- y[fit]
  rf <- ranger(y ~ ., data = d, num.trees = NTREE)
  list(model = rf, oob = rf$predictions)
}
predict_mean <- function(mfit, X, idx)
  as.numeric(predict(mfit$model, data = mean_frame(X, idx))$predictions)

fit_covariance <- function(S, fit, r_oob) {
  th <- BRISC_estimation(S[fit, , drop = FALSE], r_oob,
                         matrix(1, length(fit), 1), verbose = FALSE)$Theta
  list(sigma2 = as.numeric(th["sigma.sq"]), phi = as.numeric(th["phi"]),
       tau2 = as.numeric(th["tau.sq"]),
       c00 = as.numeric(th["sigma.sq"]) + as.numeric(th["tau.sq"]))
}
cov_mat <- function(A, th)
  th$sigma2 * exp(-th$phi * as.matrix(dist(A))) + th$tau2 * diag(nrow(A))

krige <- function(S_cond, R, Sig_inv_r, S_new, th) {
  m <- numeric(nrow(S_new)); v <- numeric(nrow(S_new))
  for (i in seq_len(nrow(S_new))) {
    d <- sqrt(rowSums((S_cond - matrix(S_new[i, ], nrow(S_cond), 2, byrow = TRUE))^2))
    k <- th$sigma2 * exp(-th$phi * d); u <- backsolve(R, k, transpose = TRUE)
    m[i] <- sum(k * Sig_inv_r); v[i] <- th$c00 - sum(u^2)
  }
  list(m = m, v = pmax(v, .Machine$double.eps))
}

## the conformal interval from calibration scores, a centre and a scale
conformal_interval <- function(scores, centre, scale) {
  A <- sort(abs(scores)); m <- length(A)
  do.call(rbind, lapply(ALPHA_RUN, function(al) {
    k <- ceiling((1 - al) * (m + 1)); q <- if (k > m) Inf else A[k]
    data.frame(alpha = al, point = seq_along(centre),
               lower = centre - q * scale, upper = centre + q * scale,
               q = q, sd = scale) }))
}

## ---- naive split conformal --------------------------------------------------
method_naive <- function(y, X, S, fit, cal, test) {
  mfit <- fit_mean(X, y, fit)
  list(scores = y[cal] - predict_mean(mfit, X, cal),
       centre = predict_mean(mfit, X, test), scale = rep(1, length(test)))
}

## ---- the whitened interval --------------------------------------------------
## Orderings of the calibration block; each depends on the locations only.
order_cal <- function(Sc, how = "none") {
  if (how == "none")   return(seq_len(nrow(Sc)))
  if (how == "coord")  return(order(rowSums(Sc)))
  if (how == "random") return(sample.int(nrow(Sc)))
  D <- as.matrix(dist(Sc)); n <- nrow(Sc)          # maximin
  o <- integer(n); o[1] <- which.min(rowSums((Sc - rep(colMeans(Sc), each = n))^2))
  d <- D[, o[1]]
  for (k in 2:n) { d[o[seq_len(k - 1)]] <- -Inf
                   o[k] <- which.max(d); d <- pmin(d, D[, o[k]]) }
  o
}

method_w_fc <- function(y, X, S, fit, cal, test, ord = "none") {
  mfit  <- fit_mean(X, y, fit)
  r_oob <- y[fit] - mfit$oob
  th    <- fit_covariance(S, fit, r_oob)
  r_cal <- y[cal] - predict_mean(mfit, X, cal)
  nf <- length(fit); nc <- length(cal)
  pc <- order_cal(S[cal, , drop = FALSE], ord)
  cal <- cal[pc]; r_cal <- r_cal[pc]
  S_fc <- rbind(S[fit, , drop = FALSE], S[cal, , drop = FALSE])
  r_fc <- c(r_oob, r_cal)
  R <- chol(cov_mat(S_fc, th))
  e <- as.numeric(backsolve(R, r_fc, transpose = TRUE))[(nf + 1):(nf + nc)]
  Sir <- backsolve(R, backsolve(R, r_fc, transpose = TRUE))
  kt  <- krige(S_fc, R, Sir, S[test, , drop = FALSE], th)
  ## the diagnostic observations, if any, each transformed as a target given F and C
  dg <- list()
  if (!is.null(attr(y, "diag_idx"))) {
    di <- attr(y, "diag_idx")
    kd <- krige(S_fc, R, Sir, S[di, , drop = FALSE], th)
    zd <- (y[di] - predict_mean(mfit, X, di) - kd$m) / sqrt(kd$v)
    dg$kappa_diag <- sqrt(mean(zd^2))
  }
  list(scores = e, centre = predict_mean(mfit, X, test) + kt$m, scale = sqrt(kt$v),
       diag = dg)
}

## The NNGP factor keeps the fitting-then-calibration order, so each residual is
## conditioned on its K nearest predecessors.
method_w_fc_nngp <- function(y, X, S, fit, cal, test) {
  mfit  <- fit_mean(X, y, fit)
  r_oob <- y[fit] - mfit$oob
  th    <- fit_covariance(S, fit, r_oob)
  r_cal <- y[cal] - predict_mean(mfit, X, cal)
  nf <- length(fit); nc <- length(cal)
  S_fc <- rbind(S[fit, , drop = FALSE], S[cal, , drop = FALSE])
  r_fc <- c(r_oob, r_cal)
  e <- as.numeric(BRISC_decorrelation(
         S_fc, sim = matrix(r_fc, ncol = 1), sigma.sq = th$sigma2, tau.sq = th$tau2,
         phi = th$phi, n.neighbors = NNGP_K, cov.model = "exponential",
         stabilization = FALSE, verbose = FALSE)$output.data)[(nf + 1):(nf + nc)]
  St <- S[test, , drop = FALSE]
  nn <- nn2(S_fc, query = St, k = min(NNGP_K, nrow(S_fc)))$nn.idx
  mu <- numeric(nrow(St)); vv <- numeric(nrow(St))
  for (i in seq_len(nrow(St))) {
    nb  <- nn[i, ]
    Cnn <- cov_mat(S_fc[nb, , drop = FALSE], th)
    d0  <- sqrt(rowSums((S_fc[nb, , drop = FALSE] -
                         matrix(St[i, ], length(nb), 2, byrow = TRUE))^2))
    c0  <- th$sigma2 * exp(-th$phi * d0)
    a   <- tryCatch(solve(Cnn, c0), error = function(e) qr.solve(Cnn, c0))
    mu[i] <- sum(a * r_fc[nb]); vv[i] <- max(th$c00 - sum(c0 * a), 1e-12)
  }
  list(scores = e, centre = predict_mean(mfit, X, test) + mu, scale = sqrt(vv))
}

## ---- GSCP and LSCP -----------------------------------------------------------
## Bandwidth selection by the interval score, as in their select_eta(), with the
## validation points scored against the other calibration points and the grid
## extended to half the largest distance between calibration points.
eta_grid <- function(Sc, n = 40L) {
  md <- min(dist(Sc))
  seq(3 * md, max(20 * md, 0.5 * max(dist(Sc))), length.out = n)
}
interval_score <- function(lower, upper, y, alpha)
  (upper - lower) + (2 / alpha) * (lower - y) * (y < lower) +
    (2 / alpha) * (y - upper) * (y > upper)

select_eta_split <- function(delta, Sc, nu_c, al, n_val = ETA_VAL_N) {
  mm <- length(delta)
  if (mm < 20) return(NA_real_)
  grid  <- eta_grid(Sc)
  e_abs <- delta * sqrt(nu_c)
  D     <- as.matrix(dist(Sc))
  jv    <- if (mm <= n_val) seq_len(mm) else sort(sample.int(mm, n_val))
  sc <- vapply(grid, function(eta) {
    s <- vapply(jv, function(j) {
      oth <- seq_len(mm)[-j]
      ee  <- exp(-D[oth, j]^2 / (2 * eta^2))
      w   <- ee / (1 + sum(ee)); w0 <- 1 / (1 + sum(ee))
      target <- floor((length(oth) + 1) * al) / (length(oth) + 1) - w0
      if (target <= 0) return(Inf)
      dl <- delta[oth]; o <- order(dl)
      tl <- rev(cumsum(rev(w[o])))
      k  <- suppressWarnings(max(which(tl >= target)))
      q  <- if (!is.finite(k)) 0 else dl[o][k]
      h  <- q * sqrt(nu_c[j])
      interval_score(-h, h, e_abs[j], al)
    }, numeric(1))
    mean(s)
  }, numeric(1))
  if (all(!is.finite(sc))) return(NA_real_)
  grid[which.min(sc)]
}

## eta = NA gives GSCP, a number gives LSCP at that bandwidth, "auto" selects it.
method_scp <- function(y, X, S, fit, cal, test, eta = NA) {
  m <- ranger(y ~ ., data = cbind(mean_frame(X, fit), y = y[fit]), num.trees = NTREE)
  p_in <- predict(m, data = mean_frame(X, fit))$predictions
  mu_fit <- m$predictions; mu_fit[is.na(mu_fit)] <- p_in[is.na(mu_fit)]
  mu_cal <- predict(m, data = mean_frame(X, cal))$predictions
  mu_test <- predict(m, data = mean_frame(X, test))$predictions
  r_fit <- y[fit] - mu_fit; r_cal <- y[cal] - mu_cal
  Sf <- as.matrix(S[fit, , drop = FALSE]); Sc <- as.matrix(S[cal, , drop = FALSE])
  St <- as.matrix(S[test, , drop = FALSE]); mm <- length(cal)

  ## the covariance routine of the scp package, with the nugget floored
  th <- get_theta(Sf, r_fit, dists = seq(0.01, 0.2, 0.01))
  if (is.null(th) || any(!is.finite(th))) stop("get_theta failed on the fitting block")
  if (th[1] < 1e-3 * th[2]) th[1] <- 1e-3 * th[2]
  c00 <- th[1] + th[2]

  ## each calibration residual and each target kriged from the fitting block
  Ca <- scp:::mat_cov(as.matrix(dist(rbind(Sf, Sc, St))), th)
  iF <- seq_along(fit); iC <- length(fit) + seq_along(cal)
  iT <- length(fit) + mm + seq_along(test)
  L  <- chol(Ca[iF, iF, drop = FALSE])
  av <- backsolve(L, backsolve(L, r_fit, transpose = TRUE))
  kr <- function(ix) { K <- Ca[iF, ix, drop = FALSE]
    U <- backsolve(L, K, transpose = TRUE)
    list(m = as.numeric(crossprod(K, av)), nu = pmax(c00 - colSums(U^2), 1e-12)) }
  kc <- kr(iC); kt <- kr(iT)
  delta <- abs(r_cal - kc$m) / sqrt(kc$nu)

  ## GSCP uses the order statistic. LSCP keeps the largest score whose weighted upper
  ## tail reaches floor((m+1) alpha)/(m+1) - w0, the cut-off of the scp package.
  thr <- function(al, i, eta) {
    if (!is.finite(eta)) { k <- ceiling((1 - al) * (mm + 1))
                           return(if (k > mm) Inf else sort(delta)[k]) }
    d <- sqrt(rowSums((Sc - matrix(St[i, ], mm, 2, byrow = TRUE))^2))
    e <- exp(-d^2 / (2 * eta^2)); w <- e / (1 + sum(e)); w0 <- 1 / (1 + sum(e))
    target <- floor((mm + 1) * al) / (mm + 1) - w0
    if (target <= 0) return(Inf)
    o <- order(delta); tl <- rev(cumsum(rev(w[o])))
    k <- suppressWarnings(max(which(tl >= target)))
    if (!is.finite(k)) 0 else delta[o][k] }

  eta_al <- if (identical(eta, "auto"))
      vapply(ALPHA_RUN, function(al) select_eta_split(delta, Sc, kc$nu, al), numeric(1))
    else rep(eta, length(ALPHA_RUN))
  iv <- do.call(rbind, lapply(seq_along(ALPHA_RUN), function(j) {
    al <- ALPHA_RUN[j]
    q <- vapply(seq_along(test), function(i) thr(al, i, eta_al[j]), numeric(1))
    ctr <- mu_test + kt$m; half <- q * sqrt(kt$nu)
    data.frame(alpha = al, point = seq_along(test), lower = ctr - half, upper = ctr + half) }))
  list(intervals = iv, diag = list(eta_10 = eta_al[ALPHA_RUN == 0.10]))
}

METHODS <- list(
  naive             = method_naive,
  gscp_rf_oob       = function(...) method_scp(...),
  lscp_rf_oob       = function(...) method_scp(..., eta = ETA_LOCAL),
  lscp_rf_oob_adapt = function(...) method_scp(..., eta = "auto"),
  w_fc              = method_w_fc,
  w_fc_nngp         = method_w_fc_nngp,
  w_fc_coordord     = function(...) method_w_fc(..., ord = "coord"),
  w_fc_maximin      = function(...) method_w_fc(..., ord = "maximin"),
  w_fc_randord      = function(...) method_w_fc(..., ord = "random"))

## row-bind data frames whose columns differ, filling the missing ones with NA
rbind_fill <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  nm <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) { for (k in setdiff(nm, names(d))) d[[k]] <- NA
                                           d[nm] }))
}

## ---- run one method -----------------------------------------------------------
## Returns one row per test point at nominal 0.90, with coverage, width, the conformal
## threshold, and the distance to the nearest fitting observation.
run_method <- function(method, y, X, S, fit, cal, test) {
  res <- METHODS[[method]](y, X, S, fit, cal, test)
  iv  <- if (!is.null(res$intervals)) res$intervals
         else conformal_interval(res$scores, res$centre, res$scale)
  iv  <- iv[abs(iv$alpha - 0.10) < 1e-8, ]
  yt  <- y[test][iv$point]
  Sf  <- S[fit, , drop = FALSE]
  dist_fit <- vapply(test, function(i)
    min(sqrt(rowSums((Sf - matrix(S[i, ], nrow(Sf), 2, byrow = TRUE))^2))), numeric(1))
  out <- data.frame(method = method, point = iv$point,
                    covered = as.integer(yt >= iv$lower & yt <= iv$upper),
                    width = iv$upper - iv$lower,
                    q = if (is.null(iv$q)) NA_real_ else iv$q,
                    dist_fit = dist_fit[iv$point])
  for (nm in names(res$diag)) out[[nm]] <- res$diag[[nm]]
  out
}

## ---- summaries over replications ------------------------------------------------
## One row per replication: coverage over the test points (an unbounded interval
## counts as covering), mean width over the bounded intervals, number unbounded.
per_rep <- function(pt) {
  keys <- intersect(c("design", "method", "rep"), names(pt))
  out <- aggregate(pt[c("covered")], pt[keys], mean)
  out$width <- aggregate(pt$width, pt[keys], function(w) mean(w[is.finite(w)]))$x
  out$n_unbounded <- aggregate(pt$width, pt[keys], function(w) sum(!is.finite(w)))$x
  out$n_test <- aggregate(pt$width, pt[keys], length)$x
  for (nm in intersect(c("kappa_diag", "q", "eta_10", "gain"), names(pt)))
    out[[nm]] <- aggregate(pt[[nm]], pt[keys], function(z) z[1])$x
  out[order(out$design, out$method, out$rep), ]
}

## Coverage and width averaged over the replications with at least one bounded
## interval; reps counts all replications, n those that enter the averages.
summarise_reps <- function(pr) {
  ok <- pr[is.finite(pr$width), ]
  out <- aggregate(cbind(coverage = ok$covered, width = ok$width),
                   ok[c("design", "method")], mean)
  out$n <- aggregate(ok$width, ok[c("design", "method")], length)$x
  nu <- aggregate(pr$n_unbounded, pr[c("design", "method")], sum)
  names(nu)[3] <- "n_unbounded"
  rp <- aggregate(pr$rep, pr[c("design", "method")], length); names(rp)[3] <- "reps"
  out <- merge(merge(rp, out, all.x = TRUE), nu)
  out[order(out$design, out$method), ]
}

## The CVR: coefficient of variation of the per-replication mean width on a fixed
## design, divided by that of naive split conformal with m independent standard normal
## scores. Replications with an unbounded interval are left out.
cv0 <- function(m, alpha, B = 1e5, chunk = 5000) {
  k <- ceiling((1 - alpha) * (m + 1)); s1 <- 0; s2 <- 0; nn <- 0
  for (start in seq(1, B, by = chunk)) {
    b <- min(chunk, B - start + 1)
    w <- 2 * apply(matrix(abs(rnorm(b * m)), b, m), 1, function(z) sort.int(z, partial = k)[k])
    s1 <- s1 + sum(w); s2 <- s2 + sum(w^2); nn <- nn + b
  }
  mu <- s1 / nn
  sqrt(max(s2 / nn - mu^2, 0) * nn / (nn - 1)) / mu
}
summarise_cvr <- function(pr, m = 166) {
  set.seed(20260824); invisible(cv0(m, 0.05)); c0 <- cv0(m, 0.10)
  pr <- pr[pr$n_unbounded == 0, ]
  out <- aggregate(pr$width, pr[c("design", "method")], function(w) sd(w) / mean(w))
  names(out)[3] <- "cv"; out$cvr <- out$cv / c0
  out
}

## Run the given methods on replications 1..nrep of each design. Each replication is
## generated by make_rep() of Designs.R and seeded by design and replication, so every
## method sees the same data and random stream.
run_designs <- function(designs, methods, nrep) {
  rbind_fill(lapply(designs, function(scen) {
    cat(sprintf("  %s\n", scen))
    pts <- parallel::mclapply(seq_len(nrep), function(b) {
      d <- make_rep(scen, b)
      y <- d$y; if (!is.null(d$diag_idx)) attr(y, "diag_idx") <- d$diag_idx
      rbind_fill(lapply(methods, function(mt) {
        set.seed(900000 + 131 * b + sum(utf8ToInt(scen)))
        p <- run_method(mt, y, d$x, d$s, d$fit, d$cal, d$test)
        p$rep <- b; p$gain <- d$gain; p })) }, mc.cores = N_CORES)
    pts <- rbind_fill(pts); pts$design <- sub("^[A-Z]_", "", scen); pts }))
}
