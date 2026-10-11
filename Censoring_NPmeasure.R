# ============================================================
# SIMULATION: INDEPENDENCE TESTING UNDER RIGHT CENSORING
# ============================================================
# Implemented:
#   Naive Spearman, Kendall, distance covariance, HSIC
#   IPCW Spearman, Kendall, distance covariance, HSIC
#
# Not included:
#   IPCW spacing statistic (formula must first be finalized)
#
# Calibration:
#   Parametric null bootstrap using the known simulation DGP.
#   This is an oracle simulation benchmark, not a general
#   plug-in test for an unknown real-data distribution.
# ============================================================

if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
library(survival)

set.seed(2026)

# ------------------------------------------------------------
# 1. Complete-data generator
# ------------------------------------------------------------

generate_error <- function(n, error) {
  switch(
    error,
    normal  = rnorm(n),
    t5      = rt(n, df = 5) / sqrt(5 / 3),
    laplace = rexp(n) *
      sample(c(-1, 1), n, replace = TRUE) / sqrt(2),
    stop("Unknown error distribution: ", error)
  )
}

generate_data <- function(n, model, theta, error) {
  X <- if (model == "sinusoidal") {
    runif(n, -1, 1)
  } else {
    rnorm(n)
  }
  
  eps <- generate_error(n, error)
  
  signal <- switch(
    model,
    linear         = theta * X,
    cubic          = theta * X^3,
    quadratic      = theta * X^2,
    sinusoidal     = theta * sin(2 * pi * X),
    heteroscedastic = (sqrt(1 + theta * X^2) - 1) * eps,
    stop("Unknown model: ", model)
  )
  
  # The heteroscedastic model requires a variance effect,
  # while the other models use an additive error.
  Y <- if (model == "heteroscedastic") {
    1000 + sqrt(1 + theta * X^2) * eps
  } else {
    1000 + signal + eps
  }
  
  if (any(Y <= 0)) {
    stop("Nonpositive event times generated.")
  }
  
  data.frame(X = X, Y = Y)
}

# ------------------------------------------------------------
# 2. Calibrate censoring intercept on a separate pilot sample
# ------------------------------------------------------------

calibrate_censor_intercept <- function(
    model, theta, error, target_censoring,
    censoring_type, gamma1 = 0.5, sigma_c = 0.25,
    n_cal = 5000) {
  
  pilot <- generate_data(n_cal, model, theta, error)
  eta <- rnorm(n_cal, sd = sigma_c)
  
  x_effect <- if (censoring_type == "conditional") {
    gamma1 * pilot$X
  } else {
    rep(0, n_cal)
  }
  
  # C = exp(gamma0 + gamma1*X + eta)
  rate <- function(gamma0) {
    C <- exp(gamma0 + x_effect + eta)
    mean(C < pilot$Y)
  }
  
  lower <- min(log(pilot$Y) - x_effect - eta) - 10
  upper <- max(log(pilot$Y) - x_effect - eta) + 10
  
  uniroot(
    function(g) rate(g) - target_censoring,
    interval = c(lower, upper)
  )$root
}

generate_censored_data <- function(
    dat, gamma0, censoring_type,
    gamma1 = 0.5, sigma_c = 0.25) {
  
  x_effect <- if (censoring_type == "conditional") {
    gamma1 * dat$X
  } else {
    rep(0, nrow(dat))
  }
  
  C <- exp(gamma0 + x_effect + rnorm(nrow(dat), sd = sigma_c))
  T <- pmin(dat$Y, C)
  Delta <- as.integer(dat$Y <= C)
  
  data.frame(X = dat$X, T = T, Delta = Delta)
}

# ------------------------------------------------------------
# 3. IPCW estimation
# ------------------------------------------------------------

estimate_ipcw <- function(dat, censoring_type,
                          floor_g = 0.05) {
  n <- nrow(dat)
  
  if (sum(dat$Delta) == 0L) {
    return(list(
      ok = FALSE,
      reason = "No observed events",
      weights = rep(NA_real_, n),
      Ghat = rep(NA_real_, n),
      n_eff = NA_real_,
      max_weight = NA_real_
    ))
  }
  
  if (censoring_type == "independent") {
    # Treat censoring as the event in the censoring KM estimator.
    fit <- survfit(Surv(T, 1 - Delta) ~ 1, data = dat)
    
    Ghat <- rep(1, n)
    
    if (length(fit$time) > 0L && any(fit$n.event > 0)) {
      idx <- findInterval(dat$T, fit$time)
      use <- idx > 0L
      Ghat[use] <- fit$surv[idx[use]]
    }
    
  } else {
    # Conditional log-normal AFT model for censoring time.
    # An observed event time is a right-censoring bound for C.
    fit <- tryCatch(
      survreg(
        Surv(T, 1 - Delta) ~ X,
        data = dat,
        dist = "lognormal"
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit) || !is.finite(fit$scale) ||
        fit$scale <= 0) {
      return(list(
        ok = FALSE,
        reason = "Conditional censoring model failed",
        weights = rep(NA_real_, n),
        Ghat = rep(NA_real_, n),
        n_eff = NA_real_,
        max_weight = NA_real_
      ))
    }
    
    mu <- as.vector(predict(fit, newdata = dat, type = "lp"))
    Ghat <- plnorm(
      dat$T,
      meanlog = mu,
      sdlog = fit$scale,
      lower.tail = FALSE
    )
  }
  
  Ghat <- pmax(Ghat, floor_g)
  W_raw <- dat$Delta / Ghat
  total <- sum(W_raw)
  
  if (!is.finite(total) || total <= 0) {
    return(list(
      ok = FALSE,
      reason = "Weights cannot be normalized",
      weights = rep(NA_real_, n),
      Ghat = Ghat,
      n_eff = NA_real_,
      max_weight = NA_real_
    ))
  }
  
  W <- W_raw / total
  
  list(
    ok = TRUE,
    reason = NA_character_,
    weights = W,
    Ghat = Ghat,
    n_eff = 1 / sum(W^2),
    max_weight = max(W_raw)
  )
}

# ------------------------------------------------------------
# 4. Dependence statistics
# ------------------------------------------------------------

weighted_ecdf_at <- function(x, w) {
  vapply(x, function(z) sum(w[x <= z]), numeric(1))
}

ipcw_spearman <- function(x, y, w) {
  Fx <- weighted_ecdf_at(x, w)
  Fy <- weighted_ecdf_at(y, w)
  12 * sum(w * Fx * Fy) - 3
}

ipcw_kendall <- function(x, y, w) {
  n <- length(x)
  num <- den <- 0
  
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    wij <- w[i] * w[j]
    
    num <- num + sum(
      wij * sign((x[i] - x[j]) * (y[i] - y[j]))
    )
    den <- den + sum(wij)
  }
  
  if (den <= 0) return(NA_real_)
  num / den
}

weighted_dcov2 <- function(x, y, w) {
  w <- w / sum(w)
  
  A <- as.matrix(dist(x))
  B <- as.matrix(dist(y))
  
  a <- as.vector(A %*% w)
  b <- as.vector(B %*% w)
  
  Ac <- A - outer(a, rep(1, length(w))) -
    outer(rep(1, length(w)), a) + sum(w * a)
  
  Bc <- B - outer(b, rep(1, length(w))) -
    outer(rep(1, length(w)), b) + sum(w * b)
  
  sum(outer(w, w) * Ac * Bc)
}

gaussian_kernel <- function(x, bandwidth = NULL) {
  D2 <- as.matrix(dist(x))^2
  
  if (is.null(bandwidth)) {
    positive <- sqrt(D2[D2 > 0])
    bandwidth <- if (length(positive)) median(positive) else 1
  }
  
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    bandwidth <- 1
  }
  
  exp(-D2 / (2 * bandwidth^2))
}

weighted_hsic <- function(K, L, w) {
  w <- w / sum(w)
  
  center_kernel <- function(M, w) {
    r <- as.vector(M %*% w)
    M - outer(r, rep(1, length(w))) -
      outer(rep(1, length(w)), r) + sum(w * r)
  }
  
  Kc <- center_kernel(K, w)
  Lc <- center_kernel(L, w)
  
  sum(outer(w, w) * Kc * Lc)
}

calculate_statistics <- function(dat, w) {
  x <- dat$X
  y <- dat$T
  n <- nrow(dat)
  w0 <- rep(1 / n, n)
  
  c(
    naive_spearman = unname(cor(x, y, method = "spearman")),
    naive_kendall  = unname(cor(x, y, method = "kendall")),
    naive_dcov2    = weighted_dcov2(x, y, w0),
    naive_hsic     = weighted_hsic(
      gaussian_kernel(x), gaussian_kernel(y), w0
    ),
    ipcw_spearman  = ipcw_spearman(x, y, w),
    ipcw_kendall   = ipcw_kendall(x, y, w),
    ipcw_dcov2     = weighted_dcov2(x, y, w),
    ipcw_hsic      = weighted_hsic(
      gaussian_kernel(x), gaussian_kernel(y), w
    )
  )
}

# ------------------------------------------------------------
# 5. Parametric null bootstrap
# ------------------------------------------------------------
# Generates X and Y independently under theta = 0, using the
# scenario's censoring intercept and specified censoring model.
# The IPCW censoring model is refitted for each bootstrap sample.
# ------------------------------------------------------------

bootstrap_test <- function(
    dat, model, error, censoring_type, gamma0,
    B = 500, gamma1 = 0.5, sigma_c = 0.25,
    floor_g = 0.05, max_attempt_factor = 20) {
  
  ipcw_obs <- estimate_ipcw(dat, censoring_type, floor_g)
  
  if (!ipcw_obs$ok) {
    return(list(
      ok = FALSE,
      reason = ipcw_obs$reason,
      statistics = NULL,
      p_values = NULL,
      n_eff = NA_real_,
      max_weight = NA_real_,
      bootstrap_failure_rate = NA_real_
    ))
  }
  
  obs <- calculate_statistics(dat, ipcw_obs$weights)
  
  boot <- matrix(
    NA_real_, nrow = B, ncol = length(obs),
    dimnames = list(NULL, names(obs))
  )
  
  accepted <- 0L
  attempts <- 0L
  max_attempts <- max_attempt_factor * B
  
  while (accepted < B && attempts < max_attempts) {
    attempts <- attempts + 1L
    
    # The null is generated directly from the known simulation DGP.
    null_complete <- generate_data(
      n = nrow(dat),
      model = model,
      theta = 0,
      error = error
    )
    
    null_obs <- generate_censored_data(
      null_complete,
      gamma0 = gamma0,
      censoring_type = censoring_type,
      gamma1 = gamma1,
      sigma_c = sigma_c
    )
    
    ipcw_b <- estimate_ipcw(
      null_obs, censoring_type, floor_g
    )
    
    # No observed events or failed conditional fit:
    # reject this bootstrap draw and try another.
    if (!ipcw_b$ok) next
    
    stat_b <- tryCatch(
      calculate_statistics(null_obs, ipcw_b$weights),
      error = function(e) NULL
    )
    
    if (is.null(stat_b) || any(!is.finite(stat_b))) next
    
    accepted <- accepted + 1L
    boot[accepted, ] <- stat_b
  }
  
  if (accepted < B) {
    return(list(
      ok = FALSE,
      reason = "Could not obtain enough valid bootstrap draws",
      statistics = obs,
      p_values = NULL,
      n_eff = ipcw_obs$n_eff,
      max_weight = ipcw_obs$max_weight,
      bootstrap_failure_rate = 1 - accepted / attempts
    ))
  }
  
  signed_methods <- c(
    "naive_spearman", "naive_kendall",
    "ipcw_spearman", "ipcw_kendall"
  )
  
  pvals <- vapply(names(obs), function(m) {
    if (m %in% signed_methods) {
      exceed <- abs(boot[, m]) >= abs(obs[m])
    } else {
      exceed <- boot[, m] >= obs[m]
    }
    (1 + sum(exceed)) / (B + 1)
  }, numeric(1))
  
  list(
    ok = TRUE,
    reason = NA_character_,
    statistics = obs,
    p_values = pvals,
    n_eff = ipcw_obs$n_eff,
    max_weight = ipcw_obs$max_weight,
    bootstrap_failure_rate = 1 - B / attempts
  )
}

# ------------------------------------------------------------
# 6. One Monte Carlo replication
# ------------------------------------------------------------

simulate_one <- function(
    n, model, theta, error, censoring,
    censoring_type, gamma0,
    B = 500, gamma1 = 0.5, sigma_c = 0.25,
    floor_g = 0.05) {
  
  complete <- generate_data(n, model, theta, error)
  
  dat <- generate_censored_data(
    complete,
    gamma0 = gamma0,
    censoring_type = censoring_type,
    gamma1 = gamma1,
    sigma_c = sigma_c
  )
  
  test <- bootstrap_test(
    dat,
    model = model,
    error = error,
    censoring_type = censoring_type,
    gamma0 = gamma0,
    B = B,
    gamma1 = gamma1,
    sigma_c = sigma_c,
    floor_g = floor_g
  )
  
  if (!test$ok) {
    return(data.frame(
      ok = FALSE,
      reason = test$reason,
      n_eff = test$n_eff,
      max_weight = test$max_weight,
      achieved_censoring = mean(dat$Delta == 0),
      bootstrap_failure_rate = test$bootstrap_failure_rate
    ))
  }
  
  out <- c(
    as.list(test$statistics),
    as.list(setNames(
      test$p_values,
      paste0("p_", names(test$p_values))
    ))
  )
  
  data.frame(
    ok = TRUE,
    reason = NA_character_,
    as.data.frame(out, check.names = FALSE),
    n_eff = test$n_eff,
    max_weight = test$max_weight,
    achieved_censoring = mean(dat$Delta == 0),
    bootstrap_failure_rate = test$bootstrap_failure_rate,
    check.names = FALSE
  )
}

# ------------------------------------------------------------
# 7. Run simulation
# ------------------------------------------------------------

run_simulation <- function(
    design, N_MC = 2000, B = 500,
    n_cal = 5000, gamma1 = 0.5,
    sigma_c = 0.25, floor_g = 0.05,
    seed = 2026) {
  
  set.seed(seed)
  
  all_results <- vector("list", nrow(design) * N_MC)
  row_id <- 1L
  gamma_cache <- new.env(parent = emptyenv())
  
  for (s in seq_len(nrow(design))) {
    cfg <- design[s, ]
    
    key <- paste(
      cfg$model, cfg$theta, cfg$error,
      cfg$censoring, cfg$censoring_type,
      sep = "|"
    )
    
    if (!exists(key, envir = gamma_cache, inherits = FALSE)) {
      g0 <- calibrate_censor_intercept(
        model = cfg$model,
        theta = cfg$theta,
        error = cfg$error,
        target_censoring = cfg$censoring,
        censoring_type = cfg$censoring_type,
        gamma1 = gamma1,
        sigma_c = sigma_c,
        n_cal = n_cal
      )
      assign(key, g0, envir = gamma_cache)
    }
    
    gamma0 <- get(key, envir = gamma_cache)
    
    message(
      "Scenario ", s, "/", nrow(design),
      ": n=", cfg$n,
      ", model=", cfg$model,
      ", theta=", cfg$theta,
      ", censoring=", cfg$censoring,
      ", mechanism=", cfg$censoring_type
    )
    
    for (r in seq_len(N_MC)) {
      ans <- tryCatch(
        simulate_one(
          n = cfg$n,
          model = cfg$model,
          theta = cfg$theta,
          error = cfg$error,
          censoring = cfg$censoring,
          censoring_type = cfg$censoring_type,
          gamma0 = gamma0,
          B = B,
          gamma1 = gamma1,
          sigma_c = sigma_c,
          floor_g = floor_g
        ),
        error = function(e) {
          data.frame(
            ok = FALSE,
            reason = conditionMessage(e),
            n_eff = NA_real_,
            max_weight = NA_real_,
            achieved_censoring = NA_real_,
            bootstrap_failure_rate = NA_real_
          )
        }
      )
      
      meta <- data.frame(
        scenario = s,
        replication = r,
        n = cfg$n,
        model = cfg$model,
        theta = cfg$theta,
        error = cfg$error,
        censoring_target = cfg$censoring,
        censoring_type = cfg$censoring_type
      )
      
      all_results[[row_id]] <- cbind(meta, ans)
      row_id <- row_id + 1L
    }
  }
  
  do.call(rbind, all_results)
}

# ------------------------------------------------------------
# 8. Summarize empirical size and power
# ------------------------------------------------------------

summarize_simulation <- function(results, alpha = 0.05) {
  pcols <- grep("^p_", names(results), value = TRUE)
  
  if (!length(pcols)) {
    stop("No p-value columns found. Check simulation output.")
  }
  
  group_cols <- c(
    "n", "model", "theta", "error",
    "censoring_target", "censoring_type"
  )
  
  group_key <- interaction(
    results[, group_cols],
    drop = TRUE,
    lex.order = TRUE
  )
  
  groups <- split(results, group_key)
  output <- list()
  k <- 1L
  
  for (d in groups) {
    valid <- d$ok %in% TRUE
    
    for (pc in pcols) {
      p <- d[[pc]][valid]
      p <- p[is.finite(p)]
      
      rejection <- if (length(p)) mean(p < alpha) else NA_real_
      mcse <- if (length(p)) {
        sqrt(rejection * (1 - rejection) / length(p))
      } else {
        NA_real_
      }
      
      output[[k]] <- data.frame(
        n = d$n[1],
        model = d$model[1],
        theta = d$theta[1],
        error = d$error[1],
        censoring_target = d$censoring_target[1],
        censoring_type = d$censoring_type[1],
        method = sub("^p_", "", pc),
        rejection_rate = rejection,
        MCSE = mcse,
        attempted_replications = nrow(d),
        valid_replications = sum(valid),
        failure_rate = mean(!valid),
        mean_n_eff = mean(d$n_eff[valid], na.rm = TRUE),
        mean_max_weight = mean(d$max_weight[valid], na.rm = TRUE),
        achieved_censoring = mean(
          d$achieved_censoring[valid], na.rm = TRUE
        )
      )
      
      k <- k + 1L
    }
  }
  
  do.call(rbind, output)
}

# ------------------------------------------------------------
# 9. Pilot run: start small to verify that everything runs
# ------------------------------------------------------------

pilot_design <- expand.grid(
  n = c(50, 100),
  model = c("linear", "quadratic"),
  theta = c(0, 0.50),
  error = "normal",
  censoring = c(0.10, 0.30),
  censoring_type = c("independent", "conditional"),
  stringsAsFactors = FALSE
)

pilot_results <- run_simulation(
  design = pilot_design,
  N_MC = 10,
  B = 99,
  n_cal = 1000
)

pilot_summary <- summarize_simulation(pilot_results)

write.csv(
  pilot_results,
  "pilot_replication_results.csv",
  row.names = FALSE
)

write.csv(
  pilot_summary,
  "pilot_simulation_summary.csv",
  row.names = FALSE
)

print(head(pilot_summary))

# ------------------------------------------------------------
# 10. Full simulation: computationally intensive
# ------------------------------------------------------------

design <- expand.grid(
  n = c(50, 100, 200, 500),
  model = c("linear", "cubic", "quadratic",
            "sinusoidal", "heteroscedastic"),
  theta = c(0, 0.25, 0.50, 0.75, 1.00, 1.50),
  error = c("normal", "t5", "laplace"),
  censoring = c(0.10, 0.30, 0.50, 0.70),
  censoring_type = c("independent", "conditional"),
  stringsAsFactors = FALSE
)

full_results <- run_simulation(
  design = design,
  N_MC = 2000,
  B = 500,
  n_cal = 5000
)

full_summary <- summarize_simulation(full_results)

write.csv(
  full_results,
  "full_simulation_replication_results.csv",
  row.names = FALSE
)

write.csv(
  full_summary,
  "full_simulation_summary.csv",
  row.names = FALSE
)