# ==============================================================================
# PEDAGOGICAL R IMPLEMENTATION
#
# From distribution-free CDF bands to confidence regions for a
# location--scale family.
#
# The code assumes that inputs have the advertised form: n and B are
# positive integers, 0 < alpha < 1, probability bounds are ordered, and
# qfun is a vectorized quantile function.  This keeps the implementation
# focused on the statistical and geometric ideas.
# ==============================================================================


# ==============================================================================
# 1. DISTRIBUTION-FREE ORDER-STATISTIC BANDS
# ==============================================================================

# If T_1,...,T_B and a future statistic T_new are exchangeable, the kth
# calibration order statistic gives marginal coverage k/(B+1).  Choosing
# k = ceiling((B+1)(1-alpha)) therefore avoids an arbitrary interpolated
# sample quantile.
mc_rank_critical_value <- function(statistics, alpha) {
  B <- length(statistics)
  k <- ceiling((B + 1) * (1 - alpha))
  if (k > B) Inf else sort(statistics, partial = k)[k]
}


simulate_critical_value <- function(
    statistic,
    n,
    alpha = 0.05,
    B = 100000,
    seed = 123,
    ...
) {
  set.seed(seed)
  statistics <- numeric(B)
  for (b in seq_len(B)) {
    statistics[b] <- statistic(runif(n), ...)
  }
  mc_rank_critical_value(statistics, alpha)
}


# ------------------------------------------------------------------------------
# Kolmogorov--Smirnov
# ------------------------------------------------------------------------------

ks_statistic <- function(u) {
  u <- sort(u)
  n <- length(u)
  i <- seq_len(n)
  max(i / n - u, u - (i - 1) / n)
}


ks_bounds_mc <- function(
    n,
    alpha = 0.05,
    B = 100000,
    seed = 123,
    critical = NULL
) {
  if (is.null(critical)) {
    critical <- simulate_critical_value(
      ks_statistic, n, alpha, B, seed
    )
  }

  i <- seq_len(n)
  bounds <- data.frame(
    i = i,
    lower = pmax(0, i / n - critical),
    upper = pmin(1, (i - 1) / n + critical)
  )

  list(critical_value = critical, bounds = bounds)
}


# ------------------------------------------------------------------------------
# Dvoretzky--Kiefer--Wolfowitz
# ------------------------------------------------------------------------------

dkw_bounds <- function(n, alpha = 0.05) {
  epsilon <- sqrt(log(2 / alpha) / (2 * n))
  i <- seq_len(n)

  data.frame(
    i = i,
    lower = pmax(0, i / n - epsilon),
    upper = pmin(1, (i - 1) / n + epsilon)
  )
}


# ------------------------------------------------------------------------------
# Berk--Jones/Owen
# ------------------------------------------------------------------------------

# Bernoulli Kullback--Leibler divergence, with 0 log(0/y) = 0.
kl_bern <- function(p, q) {
  xlogy <- function(x, y) {
    value <- x * log(x / y)
    value[x == 0] <- 0
    value
  }

  xlogy(p, q) + xlogy(1 - p, 1 - q)
}


bj_statistic <- function(u) {
  u <- sort(u)
  n <- length(u)
  i <- seq_len(n)

  max(
    n * kl_bern((i - 1) / n, u),
    n * kl_bern(i / n, u)
  )
}


# Solve K(p,q) = critical/n on both sides of p.
bj_interval_at_p <- function(p, critical, n, tol = 1e-10) {
  if (is.infinite(critical)) return(c(lower = 0, upper = 1))

  target <- critical / n
  f <- function(q) kl_bern(p, q) - target

  lower <- if (p == 0) {
    0
  } else {
    uniroot(f, c(.Machine$double.eps, p), tol = tol)$root
  }

  upper <- if (p == 1) {
    1
  } else {
    uniroot(f, c(p, 1 - .Machine$double.eps), tol = tol)$root
  }

  c(lower = lower, upper = upper)
}


bj_bounds <- function(
    n,
    alpha = 0.05,
    B = 100000,
    seed = 123,
    critical = NULL
) {
  if (is.null(critical)) {
    critical <- simulate_critical_value(
      bj_statistic, n, alpha, B, seed
    )
  }

  i <- seq_len(n)
  lower <- vapply(
    i / n,
    function(p) bj_interval_at_p(p, critical, n)["lower"],
    numeric(1)
  )
  upper <- vapply(
    (i - 1) / n,
    function(p) bj_interval_at_p(p, critical, n)["upper"],
    numeric(1)
  )

  list(
    critical_value = critical,
    bounds = data.frame(i = i, lower = lower, upper = upper)
  )
}


# ------------------------------------------------------------------------------
# Duembgen--Wellner, with s = 1 and nu = 1 by default
# ------------------------------------------------------------------------------

C_nu <- function(t, nu = 1) {
  value <- rep(Inf, length(t))
  interior <- t > 0 & t < 1

  C <- log(log(exp(1) / (4 * t[interior] * (1 - t[interior]))))
  value[interior] <- C + nu * log(1 + C^2)
  value
}


# min_{v between p and q} C_nu(v).  The correction decreases toward 1/2.
C_nu_biv <- function(p, q, nu = 1) {
  lo <- pmin(p, q)
  hi <- pmax(p, q)

  ifelse(
    hi < 0.5,
    C_nu(hi, nu),
    ifelse(lo > 0.5, C_nu(lo, nu), 0)
  )
}


dw_discrepancy <- function(p, q, n, nu = 1) {
  n * kl_bern(p, q) - C_nu_biv(p, q, nu)
}


dw_statistic <- function(u, nu = 1) {
  u <- sort(u)
  n <- length(u)
  i <- seq_len(n)

  max(
    dw_discrepancy((i - 1) / n, u, n, nu),
    dw_discrepancy(i / n, u, n, nu)
  )
}


# The upper endpoint is the largest q >= p satisfying the DW inequality.
dw_upper_at_p <- function(p, n, critical, nu = 1, tol = 1e-10) {
  if (p == 1 || is.infinite(critical)) return(1)

  f <- function(q) dw_discrepancy(p, q, n, nu) - critical
  q_max <- 1 - .Machine$double.eps

  if (f(q_max) <= 0) 1 else uniroot(f, c(p, q_max), tol = tol)$root
}


dw_bounds <- function(
    n,
    alpha = 0.05,
    B = 100000,
    nu = 1,
    seed = 123,
    critical = NULL
) {
  if (is.null(critical)) {
    critical <- simulate_critical_value(
      dw_statistic, n, alpha, B, seed, nu = nu
    )
  }

  i <- seq_len(n)
  lower <- 1 - vapply(
    1 - i / n,
    function(p) dw_upper_at_p(p, n, critical, nu),
    numeric(1)
  )
  upper <- vapply(
    (i - 1) / n,
    function(p) dw_upper_at_p(p, n, critical, nu),
    numeric(1)
  )

  list(
    critical_value = critical,
    bounds = data.frame(i = i, lower = lower, upper = upper)
  )
}


# Calibrate once for a given (n, alpha), then reuse these bounds for every
# location--scale family in the outer simulation.
calibrate_bands <- function(
    n,
    alpha = 0.05,
    B = 100000,
    seed = 123
) {
  KS <- ks_bounds_mc(n, alpha, B, seed + 1)
  BJ <- bj_bounds(n, alpha, B, seed + 2)
  DW <- dw_bounds(n, alpha, B, seed = seed + 3)

  list(
    n = n,
    alpha = alpha,
    B = B,
    bounds = list(
      KS = KS$bounds,
      DKW = dkw_bounds(n, alpha),
      Berk_Jones = BJ$bounds,
      Duembgen_Wellner = DW$bounds
    ),
    critical_values = c(
      KS = KS$critical_value,
      Berk_Jones = BJ$critical_value,
      Duembgen_Wellner = DW$critical_value
    )
  )
}


# ==============================================================================
# 2. EXACT LOCATION--SCALE GEOMETRY
# ==============================================================================

# Upper envelope of lines c_j + m_j*s.  Slopes must be nondecreasing.
# The stack stores each surviving line and the point where it becomes active.
upper_envelope <- function(intercept, slope) {
  n <- length(slope)

  # Among parallel lines, only the highest intercept can be active.
  unique_c <- numeric(n)
  unique_m <- numeric(n)
  n_unique <- 0L

  for (j in seq_len(n)) {
    if (n_unique > 0L && slope[j] == unique_m[n_unique]) {
      unique_c[n_unique] <- max(unique_c[n_unique], intercept[j])
    } else {
      n_unique <- n_unique + 1L
      unique_c[n_unique] <- intercept[j]
      unique_m[n_unique] <- slope[j]
    }
  }

  unique_c <- unique_c[seq_len(n_unique)]
  unique_m <- unique_m[seq_len(n_unique)]

  hull_c <- numeric(n_unique)
  hull_m <- numeric(n_unique)
  start <- numeric(n_unique)
  size <- 0L

  for (j in seq_len(n_unique)) {
    cj <- unique_c[j]
    mj <- unique_m[j]

    while (size > 1L) {
      crossing <- (hull_c[size] - cj) / (mj - hull_m[size])
      if (crossing > start[size]) break
      size <- size - 1L
    }

    size <- size + 1L
    hull_c[size] <- cj
    hull_m[size] <- mj
    start[size] <- if (size == 1L) {
      -Inf
    } else {
      (hull_c[size - 1L] - cj) / (mj - hull_m[size - 1L])
    }
  }

  keep <- seq_len(size)
  list(
    intercept = hull_c[keep],
    slope = hull_m[keep],
    start = start[keep]
  )
}


eval_envelope <- function(envelope, s) {
  active <- findInterval(s, envelope$start)
  active <- pmax.int(1L, pmin.int(active, length(envelope$start)))
  envelope$intercept[active] + envelope$slope[active] * s
}


# Merge two already sorted breakpoint lists in linear time.
merge_sorted_unique <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  answer <- numeric(length(x) + length(y))
  i <- j <- 1L
  k <- 0L

  while (i <= length(x) || j <= length(y)) {
    if (j > length(y) || (i <= length(x) && x[i] <= y[j])) {
      value <- x[i]
      i <- i + 1L
    } else {
      value <- y[j]
      j <- j + 1L
    }

    if (k == 0L || value != answer[k]) {
      k <- k + 1L
      answer[k] <- value
    }
  }

  if (k == 0L) numeric(0) else answer[seq_len(k)]
}


# Evaluation is also linear when both the envelope starts and s are sorted.
eval_envelope_sorted <- function(envelope, s) {
  answer <- numeric(length(s))
  active <- 1L

  for (j in seq_along(s)) {
    while (active < length(envelope$start) &&
           envelope$start[active + 1L] <= s[j]) {
      active <- active + 1L
    }
    answer[j] <- envelope$intercept[active] + envelope$slope[active] * s[j]
  }

  answer
}


# H(s) = U(s)-L(s) is concave and piecewise affine.  Its nonnegative set
# is therefore one interval.  Scan the merged envelope breakpoints and
# interpolate its first and last zero crossings.
sigma_interval_from_envelopes <- function(env_L, env_V, tol = 1e-12) {
  starts_L <- env_L$start[is.finite(env_L$start) & env_L$start > 0]
  starts_V <- env_V$start[is.finite(env_V$start) & env_V$start > 0]
  knots <- c(0, merge_sorted_unique(starts_L, starts_V))

  L <- eval_envelope_sorted(env_L, knots)
  V <- eval_envelope_sorted(env_V, knots)
  gap <- -V - L
  gap_tol <- tol * pmax(1, abs(L), abs(V))
  gap[abs(gap) <= gap_tol] <- 0

  tail_slope <- -tail(env_V$slope, 1) - tail(env_L$slope, 1)
  slope_tol <- tol * max(
    1, abs(tail(env_V$slope, 1)), abs(tail(env_L$slope, 1))
  )

  feasible <- which(gap >= 0)

  if (length(feasible) == 0L) {
    if (tail_slope > slope_tol) {
      lower <- tail(knots, 1) - tail(gap, 1) / tail_slope
      return(c(lower = max(0, lower), upper = Inf))
    }
    return(c(lower = NA_real_, upper = NA_real_))
  }

  first <- feasible[1]
  last <- tail(feasible, 1)

  lower <- if (first == 1L) {
    0
  } else {
    x1 <- knots[first - 1L]
    x2 <- knots[first]
    y1 <- gap[first - 1L]
    y2 <- gap[first]
    x1 - y1 * (x2 - x1) / (y2 - y1)
  }

  upper <- if (last < length(knots)) {
    x1 <- knots[last]
    x2 <- knots[last + 1L]
    y1 <- gap[last]
    y2 <- gap[last + 1L]
    x1 - y1 * (x2 - x1) / (y2 - y1)
  } else if (tail_slope < -slope_tol) {
    tail(knots, 1) - tail(gap, 1) / tail_slope
  } else {
    Inf
  }

  c(lower = max(0, lower), upper = upper)
}


minimum_of_envelope <- function(envelope, tol = 1e-12) {
  knots <- c(0, envelope$start[is.finite(envelope$start) & envelope$start > 0])
  values <- eval_envelope_sorted(envelope, knots)
  slope_tol <- tol * max(1, abs(tail(envelope$slope, 1)))

  if (tail(envelope$slope, 1) < -slope_tol) -Inf else min(values)
}


empty_region <- function(a, b) {
  list(
    status = "empty",
    feasible = FALSE,
    bounded = TRUE,
    sigma_interval = c(lower = NA_real_, upper = NA_real_),
    sigma_length = 0,
    mu_interval = c(lower = NA_real_, upper = NA_real_),
    mu_length = 0,
    area = 0,
    a = a,
    b = b,
    boundary = NULL,
    polygon = NULL
  )
}


# Work in centered, dimensionless coordinates.  This is the actual geometric
# construction; the public wrapper below restores the original units.
location_scale_confidence_normalized <- function(
    x,
    lower,
    upper,
    qfun = qnorm,
    tol = 1e-12
) {
  n <- length(x)
  a <- qfun(lower)
  b <- qfun(upper)

  # These transformed endpoints would require an infinite parameter value.
  if (any(a == Inf) || any(b == -Inf)) return(empty_region(a, b))

  finite_a <- is.finite(a)
  finite_b <- is.finite(b)

  env_L <- if (any(finite_b)) {
    index <- rev(which(finite_b))
    upper_envelope(x[index], -b[index])
  } else NULL

  env_V <- if (any(finite_a)) {
    index <- which(finite_a)
    upper_envelope(-x[index], a[index])
  } else NULL

  # If one side has no finite constraints, the set has infinite area.
  if (is.null(env_L) || is.null(env_V)) {
    mu_lower <- if (is.null(env_L)) -Inf else minimum_of_envelope(env_L, tol)
    mu_upper <- if (is.null(env_V)) Inf else -minimum_of_envelope(env_V, tol)

    return(list(
      status = "unbounded",
      feasible = TRUE,
      bounded = FALSE,
      sigma_interval = c(lower = 0, upper = Inf),
      sigma_length = Inf,
      mu_interval = c(lower = mu_lower, upper = mu_upper),
      mu_length = Inf,
      area = Inf,
      a = a,
      b = b,
      boundary = NULL,
      polygon = NULL
    ))
  }

  sigma_interval <- sigma_interval_from_envelopes(env_L, env_V, tol)
  if (anyNA(sigma_interval) || sigma_interval["upper"] <= 0) {
    return(empty_region(a, b))
  }

  sigma_lower <- sigma_interval["lower"]
  sigma_upper <- sigma_interval["upper"]

  # Unbounded scale interval: retain the finite breakpoints and examine the
  # final affine pieces to decide which projections and measures diverge.
  if (is.infinite(sigma_upper)) {
    starts_L <- env_L$start[is.finite(env_L$start) & env_L$start > sigma_lower]
    starts_V <- env_V$start[is.finite(env_V$start) & env_V$start > sigma_lower]
    knots <- c(sigma_lower, merge_sorted_unique(starts_L, starts_V))

    L <- eval_envelope_sorted(env_L, knots)
    U <- -eval_envelope_sorted(env_V, knots)
    gap <- U - L
    gap_tol <- tol * pmax(1, abs(L), abs(U))
    gap[abs(gap) <= gap_tol] <- 0

    gap_slope <- -tail(env_V$slope, 1) - tail(env_L$slope, 1)
    slope_tol <- tol * max(
      1, abs(tail(env_V$slope, 1)), abs(tail(env_L$slope, 1))
    )
    prefix_area <- if (length(knots) == 1L) 0 else {
      sum(0.5 * (gap[-length(gap)] + gap[-1]) * diff(knots))
    }
    area <- if (gap_slope > slope_tol || tail(gap, 1) > tail(gap_tol, 1)) {
      Inf
    } else prefix_area

    lower_slope <- tail(env_L$slope, 1)
    upper_slope <- -tail(env_V$slope, 1)
    mu_lower <- if (lower_slope < -slope_tol) -Inf else min(L)
    mu_upper <- if (upper_slope > slope_tol) Inf else max(U)

    return(list(
      status = "unbounded",
      feasible = TRUE,
      bounded = FALSE,
      sigma_interval = sigma_interval,
      sigma_length = Inf,
      mu_interval = c(lower = mu_lower, upper = mu_upper),
      mu_length = if (is.finite(mu_lower) && is.finite(mu_upper)) {
        max(0, mu_upper - mu_lower)
      } else Inf,
      area = area,
      a = a,
      b = b,
      boundary = data.frame(
        sigma = knots, lower_mu = L, upper_mu = U, width = gap
      ),
      polygon = NULL
    ))
  }

  # Bounded region: every boundary is affine between consecutive knots, so
  # trapezoidal integration is exact rather than an approximation.
  starts_L <- env_L$start[
    is.finite(env_L$start) & env_L$start > sigma_lower &
      env_L$start < sigma_upper
  ]
  starts_V <- env_V$start[
    is.finite(env_V$start) & env_V$start > sigma_lower &
      env_V$start < sigma_upper
  ]
  knots <- c(
    sigma_lower,
    merge_sorted_unique(starts_L, starts_V),
    sigma_upper
  )
  knots <- unique(knots)

  L <- eval_envelope_sorted(env_L, knots)
  U <- -eval_envelope_sorted(env_V, knots)
  gap <- U - L
  gap[gap < 0 & abs(gap) <= tol * pmax(1, abs(L), abs(U))] <- 0

  area <- if (length(knots) == 1L) 0 else {
    sum(0.5 * (gap[-length(gap)] + gap[-1]) * diff(knots))
  }

  boundary <- data.frame(
    sigma = knots,
    lower_mu = L,
    upper_mu = U,
    width = gap
  )
  polygon <- rbind(
    data.frame(sigma = knots, mu = L, boundary = "lower"),
    data.frame(sigma = rev(knots), mu = rev(U), boundary = "upper")
  )

  list(
    status = "bounded",
    feasible = TRUE,
    bounded = TRUE,
    sigma_interval = sigma_interval,
    sigma_length = sigma_upper - sigma_lower,
    mu_interval = c(lower = min(L), upper = max(U)),
    mu_length = max(U) - min(L),
    area = area,
    a = a,
    b = b,
    boundary = boundary,
    polygon = polygon
  )
}


# Centering and scaling make numerical tolerances independent of measurement
# units.  The final lines simply restore the original location and scale.
location_scale_confidence <- function(
    x,
    lower,
    upper,
    qfun = qnorm,
    tol = 1e-12
) {
  x <- sort(x)
  center <- x[(length(x) + 1L) %/% 2L]
  data_scale <- max(abs(x - center))
  if (data_scale == 0) data_scale <- 1

  fit <- location_scale_confidence_normalized(
    (x - center) / data_scale,
    lower,
    upper,
    qfun,
    tol
  )

  fit$sigma_interval <- data_scale * fit$sigma_interval
  fit$sigma_length <- data_scale * fit$sigma_length
  fit$mu_interval <- center + data_scale * fit$mu_interval
  fit$mu_length <- data_scale * fit$mu_length
  fit$area <- if (fit$area == 0) 0 else data_scale^2 * fit$area

  if (!is.null(fit$boundary)) {
    fit$boundary$sigma <- data_scale * fit$boundary$sigma
    fit$boundary$lower_mu <- center + data_scale * fit$boundary$lower_mu
    fit$boundary$upper_mu <- center + data_scale * fit$boundary$upper_mu
    fit$boundary$width <- data_scale * fit$boundary$width
  }
  if (!is.null(fit$polygon)) {
    fit$polygon$sigma <- data_scale * fit$polygon$sigma
    fit$polygon$mu <- center + data_scale * fit$polygon$mu
  }

  fit
}


# ==============================================================================
# 3. MONTE CARLO COMPARISON
# ==============================================================================

compare_bands_mc <- function(
    n = 50,
    alpha = 0.05,
    M = 10000,
    band_calibration_B = 100000,
    seed = 123,
    ls.quantile.function = qnorm,
    bands = NULL,
    progress = interactive()
) {
  calibration <- if (is.null(bands)) {
    calibrate_bands(n, alpha, band_calibration_B, seed)
  } else if (!is.null(bands$bounds)) {
    bands
  } else NULL

  if (!is.null(calibration)) bands <- calibration$bounds
  methods <- names(bands)

  area <- mu_length <- sigma_length <- matrix(
    NA_real_, M, length(methods), dimnames = list(NULL, methods)
  )
  coverage <- matrix(FALSE, M, length(methods), dimnames = list(NULL, methods))
  status <- matrix("", M, length(methods), dimnames = list(NULL, methods))

  set.seed(seed)
  for (m in seq_len(M)) {
    # Sampling on the probability scale gives both the data and a direct
    # check of the simultaneous band event.
    u <- sort(runif(n))
    x <- ls.quantile.function(u)

    for (j in seq_along(methods)) {
      band <- bands[[j]]
      coverage[m, j] <- all(band$lower <= u & u <= band$upper)

      fit <- location_scale_confidence(
        x, band$lower, band$upper, ls.quantile.function
      )
      area[m, j] <- fit$area
      mu_length[m, j] <- fit$mu_length
      sigma_length[m, j] <- fit$sigma_length
      status[m, j] <- fit$status
    }

    if (progress && m %% 1000 == 0) {
      message("Completed ", m, " of ", M)
    }
  }

  standard_error <- function(z) {
    if (all(is.finite(z))) sd(z) / sqrt(length(z)) else NA_real_
  }

  empirical_coverage <- colMeans(coverage)
  summary <- data.frame(
    method = methods,
    expected_area = colMeans(area),
    se_area = apply(area, 2, standard_error),
    expected_mu_length = colMeans(mu_length),
    se_mu_length = apply(mu_length, 2, standard_error),
    expected_sigma_length = colMeans(sigma_length),
    se_sigma_length = apply(sigma_length, 2, standard_error),
    empirical_coverage = empirical_coverage,
    se_coverage = sqrt(empirical_coverage * (1 - empirical_coverage) / M),
    n_replications = M,
    n_empty = colSums(status == "empty"),
    n_bounded = colSums(status == "bounded"),
    n_unbounded = colSums(status == "unbounded"),
    row.names = NULL
  )

  list(
    summary = summary,
    area = area,
    mu_length = mu_length,
    sigma_length = sigma_length,
    coverage = coverage,
    region_status = status,
    bands = bands,
    band_calibration = calibration
  )
}


q_laplace <- function(p) {
  ifelse(p < 0.5, log(2 * p), -log(2) - log1p(-p))
}


# This function contains every experiment reported in the paper.  Each sample
# size is calibrated once; the resulting bands are reused across families.
run_reproduction_study <- function(
    alpha = 0.05,
    M = 10000,
    calibration_B = 100000,
    seed = 42,
    save_path = NULL,
    progress = interactive()
) {
  sizes <- c(20L, 50L, 100L, 200L)
  calibrations <- setNames(
    lapply(
      sizes,
      function(n) calibrate_bands(n, alpha, calibration_B, seed + n)
    ),
    sizes
  )

  run_family <- function(n, qfun) {
    compare_bands_mc(
      n = n,
      alpha = alpha,
      M = M,
      seed = seed,
      ls.quantile.function = qfun,
      bands = calibrations[[as.character(n)]],
      progress = progress
    )
  }

  results <- list(
    normal = setNames(lapply(sizes, run_family, qfun = qnorm), sizes),
    laplace = list(`50` = run_family(50L, q_laplace)),
    student_t2 = list(`50` = run_family(50L, function(p) qt(p, df = 2)))
  )

  if (!is.null(save_path)) saveRDS(results, save_path)
  results
}


# A compact example for the notebook.  Set draw = TRUE to plot the regions.
run_illustrative_example <- function(
    n = 50,
    alpha = 0.05,
    calibration_B = 100000,
    seed = 1,
    draw = interactive()
) {
  calibration <- calibrate_bands(n, alpha, calibration_B, seed)
  set.seed(seed)
  x <- rnorm(n, mean = 2, sd = 3)

  fits <- lapply(
    calibration$bounds,
    function(band) {
      location_scale_confidence(x, band$lower, band$upper, qnorm)
    }
  )

  summary <- data.frame(
    method = names(fits),
    status = vapply(fits, `[[`, character(1), "status"),
    area = vapply(fits, `[[`, numeric(1), "area"),
    mu_length = vapply(fits, `[[`, numeric(1), "mu_length"),
    sigma_length = vapply(fits, `[[`, numeric(1), "sigma_length")
  )

  if (draw) {
    bounded <- vapply(fits, function(fit) fit$status == "bounded", logical(1))
    old_par <- par(mfrow = grDevices::n2mfrow(sum(bounded)))
    on.exit(par(old_par))

    for (method in names(fits)[bounded]) {
      polygon_data <- fits[[method]]$polygon
      plot(
        polygon_data$sigma,
        polygon_data$mu,
        type = "n",
        xlab = expression(sigma),
        ylab = expression(mu),
        main = method
      )
      polygon(polygon_data$sigma, polygon_data$mu)
    }
  }

  list(sample = x, fits = fits, summary = summary)
}
