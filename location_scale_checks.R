# ============================================================
# VERIFICATION CHECKS FOR THE LOCATION--SCALE IMPLEMENTATION
#
# Run from the project directory with
#
#   Rscript location_scale_checks.R
#
# Sourcing location_scale_reference.R is intentionally safe: it
# defines functions but does not launch examples or simulations.
# ============================================================


load_reference_functions <- function(
    path = "location_scale_reference.R"
) {
  source(path, local = .GlobalEnv)
  invisible(TRUE)
}


relative_error <- function(observed, expected) {
  max(abs(observed - expected) / pmax(1, abs(expected)))
}


# ------------------------------------------------------------
# Quadratic reference formula
#
# This is retained only as an independent correctness oracle for
# moderate n. It is not used by the production implementation.
# ------------------------------------------------------------

sigma_interval_pairwise_reference <- function(
    x,
    a,
    b,
    tol = 1e-10
) {
  ii <- which(is.finite(b))
  jj <- which(is.finite(a))

  if (length(ii) == 0L || length(jj) == 0L) {
    return(list(
      feasible = TRUE,
      interval = c(lower = 0, upper = Inf)
    ))
  }

  D <- outer(b[ii], a[jj], function(bi, aj) aj - bi)
  Delta <- outer(x[ii], x[jj], function(xi, xj) xj - xi)

  if (any(abs(D) <= tol & Delta < -tol)) {
    return(list(
      feasible = FALSE,
      interval = c(lower = NA_real_, upper = NA_real_)
    ))
  }

  lower_candidates <- Delta[D < -tol] / D[D < -tol]
  upper_candidates <- Delta[D > tol] / D[D > tol]

  lower <- max(c(0, lower_candidates))
  upper <- if (length(upper_candidates) == 0L) {
    Inf
  } else {
    min(upper_candidates)
  }

  list(
    feasible = upper > 0 && lower <= upper + tol,
    interval = c(lower = lower, upper = upper)
  )
}


make_monotone_band <- function(n, half_width) {
  p <- seq_len(n) / (n + 1)
  list(
    lower = pmax(0, p - half_width),
    upper = pmin(1, p + half_width)
  )
}


# ------------------------------------------------------------
# 1. Envelope evaluation, including one segment and the final
#    active segment. These are direct regressions for the former
#    findInterval indexing defect.
# ------------------------------------------------------------

check_envelope_evaluation <- function(tolerance = 1e-14) {
  one_line <- upper_envelope(
    intercept = 2,
    slope = -3
  )
  one_error <- max(abs(
    eval_envelope(one_line, c(0, 1, 4)) - c(2, -1, -10)
  ))

  intercept <- c(0, 1, -2)
  slope <- c(-1, 0, 2)
  multi <- upper_envelope(intercept, slope)
  grid <- c(0, 0.5, 2, 10)
  direct <- vapply(
    grid,
    function(s) max(intercept + slope * s),
    numeric(1)
  )
  final_error <- max(abs(eval_envelope(multi, grid) - direct))

  nearly_parallel <- upper_envelope(
    intercept = c(0, -5e-15),
    slope = c(0, 1e-14)
  )
  parallel_grid <- c(0, 0.5, 1)
  parallel_direct <- vapply(
    parallel_grid,
    function(s) max(c(0, -5e-15) + c(0, 1e-14) * s),
    numeric(1)
  )
  parallel_error <- max(abs(
    eval_envelope(nearly_parallel, parallel_grid) - parallel_direct
  ))

  if (one_error > tolerance || final_error > tolerance ||
      parallel_error > tolerance || length(nearly_parallel$slope) != 2L) {
    stop("Envelope evaluation regression check failed.")
  }

  data.frame(
    check = c(
      "one-line envelope",
      "final envelope segment",
      "nearly parallel lines"
    ),
    maximum_absolute_error = c(one_error, final_error, parallel_error),
    passed = TRUE
  )
}


# ------------------------------------------------------------
# 2. Compare the linear scan with the pairwise formula.
# ------------------------------------------------------------

check_against_pairwise <- function(
    cases = 1000,
    seed = 4,
    tolerance = 1e-8
) {
  set.seed(seed)
  max_error <- 0

  for (case in seq_len(cases)) {
    n <- sample(3:80, 1L)
    x <- sort(rnorm(n))
    band <- make_monotone_band(
      n,
      half_width = runif(1L, 0.05, 0.45)
    )

    fit <- location_scale_confidence(
      x,
      band$lower,
      band$upper,
      qfun = qnorm
    )

    reference <- sigma_interval_pairwise_reference(
      x,
      qnorm(band$lower),
      qnorm(band$upper)
    )

    if (!identical(fit$feasible, reference$feasible)) {
      stop("Feasibility mismatch in randomized case ", case, ".")
    }

    if (fit$feasible) {
      error <- max(
        abs(fit$sigma_interval - reference$interval),
        na.rm = TRUE
      )

      if (is.finite(error)) {
        max_error <- max(max_error, error)
      }
      if (error > tolerance) {
        stop("Scale-interval mismatch in randomized case ", case, ".")
      }
    }
  }

  data.frame(
    check = "linear scan versus pairwise reference",
    cases = cases,
    maximum_absolute_error = max_error,
    passed = TRUE
  )
}


# ------------------------------------------------------------
# 3. Compare envelope values with direct maxima and verify area
#    against numerical integration.
# ------------------------------------------------------------

check_geometry <- function(seed = 19, tolerance = 1e-8) {
  set.seed(seed)
  n <- 120
  x <- sort(rnorm(n))
  band <- make_monotone_band(n, half_width = 0.18)

  fit <- location_scale_confidence(
    x,
    band$lower,
    band$upper,
    qfun = qnorm
  )

  stopifnot(
    fit$feasible,
    fit$bounded
  )

  grid <- seq(
    fit$sigma_interval["lower"],
    fit$sigma_interval["upper"],
    length.out = 501
  )

  finite_a <- is.finite(fit$a)
  finite_b <- is.finite(fit$b)
  direct_L <- vapply(
    grid,
    function(s) max(x[finite_b] - fit$b[finite_b] * s),
    numeric(1)
  )
  direct_U <- vapply(
    grid,
    function(s) min(x[finite_a] - fit$a[finite_a] * s),
    numeric(1)
  )

  envelope_L <- approx(
    fit$boundary$sigma,
    fit$boundary$lower_mu,
    xout = grid,
    ties = "ordered"
  )$y
  envelope_U <- approx(
    fit$boundary$sigma,
    fit$boundary$upper_mu,
    xout = grid,
    ties = "ordered"
  )$y
  boundary_error <- max(
    abs(direct_L - envelope_L),
    abs(direct_U - envelope_U)
  )

  integration_grid <- seq(
    fit$sigma_interval["lower"],
    fit$sigma_interval["upper"],
    length.out = 100001L
  )
  integration_lower <- approx(
    fit$boundary$sigma,
    fit$boundary$lower_mu,
    xout = integration_grid,
    ties = "ordered"
  )$y
  integration_upper <- approx(
    fit$boundary$sigma,
    fit$boundary$upper_mu,
    xout = integration_grid,
    ties = "ordered"
  )$y
  integration_width <- integration_upper - integration_lower
  numerical_area <- sum(
    0.5 *
      (integration_width[-length(integration_width)] +
         integration_width[-1L]) *
      diff(integration_grid)
  )
  area_error <- abs(fit$area - numerical_area)

  if (boundary_error > tolerance || area_error > tolerance) {
    stop("Envelope or area verification failed.")
  }

  data.frame(
    check = c(
      "direct boundaries versus envelopes",
      "exact versus numerical area"
    ),
    maximum_absolute_error = c(boundary_error, area_error),
    passed = TRUE
  )
}


# ------------------------------------------------------------
# 4. Check affine equivariance over extreme unit changes.
# ------------------------------------------------------------

check_affine_equivariance <- function(
    seed = 29,
    tolerance = 1e-8
) {
  set.seed(seed)
  n <- 80
  x <- sort(rnorm(n))
  band <- make_monotone_band(n, half_width = 0.17)
  baseline <- location_scale_confidence(
    x, band$lower, band$upper, qnorm
  )
  stopifnot(baseline$status == "bounded")

  scales <- c(1e-12, 1e-9, 1, 1e6)
  errors <- numeric(length(scales))
  for (k in seq_along(scales)) {
    multiplier <- scales[k]
    transformed <- location_scale_confidence(
      multiplier * x,
      band$lower,
      band$upper,
      qnorm
    )
    errors[k] <- max(
      relative_error(
        transformed$sigma_interval,
        multiplier * baseline$sigma_interval
      ),
      relative_error(
        transformed$mu_interval,
        multiplier * baseline$mu_interval
      ),
      relative_error(transformed$area, multiplier^2 * baseline$area)
    )
  }

  shift <- 1e6
  translated <- location_scale_confidence(
    x + shift,
    band$lower,
    band$upper,
    qnorm
  )
  shift_error <- max(
    relative_error(translated$sigma_interval, baseline$sigma_interval),
    relative_error(translated$mu_interval, baseline$mu_interval + shift),
    relative_error(translated$area, baseline$area)
  )

  if (any(errors > tolerance) || shift_error > tolerance) {
    stop("Affine-equivariance check failed.")
  }

  data.frame(
    transformation = c(paste0("scale ", scales), "shift 1e6"),
    maximum_relative_error = c(errors, shift_error),
    passed = TRUE
  )
}


# ------------------------------------------------------------
# 5. Check explicit empty and unbounded semantics.
# ------------------------------------------------------------

check_edge_cases <- function() {
  empty <- location_scale_confidence(
    x = c(0, 1),
    lower = c(0.5, 0.5),
    upper = c(0.5, 0.5),
    qfun = qnorm
  )

  unbounded <- location_scale_confidence(
    x = c(-1, 1),
    lower = c(0, 0),
    upper = c(1, 1),
    qfun = qnorm
  )

  u_sorted <- c(0.2, 0.8)
  unbounded_covers <- all(
    c(0, 0) <= u_sorted & u_sorted <= c(1, 1)
  )

  stopifnot(
    identical(empty$status, "empty"),
    empty$area == 0,
    empty$mu_length == 0,
    empty$sigma_length == 0,
    identical(unbounded$status, "unbounded"),
    !unbounded$bounded,
    unbounded_covers
  )

  data.frame(
    check = c(
      "empty set contributes zero size",
      "unbounded set may still cover"
    ),
    passed = TRUE
  )
}


# ------------------------------------------------------------
# 6. Inversion event equals the originating statistic event.
# ------------------------------------------------------------

check_band_inversion <- function(
    n = 18,
    cases = 500,
    seed = 37
) {
  set.seed(seed)
  critical <- list(KS = 0.25, BJ = 4, DW = 4)

  i <- seq_len(n)
  bands <- list(
    KS = data.frame(
      lower = pmax(0, i / n - critical$KS),
      upper = pmin(1, (i - 1) / n + critical$KS)
    ),
    BJ = bj_bounds(n = n, critical = critical$BJ)$bounds,
    DW = dw_bounds(n = n, critical = critical$DW)$bounds
  )

  mismatches <- setNames(integer(3), names(bands))
  for (case in seq_len(cases)) {
    u <- sort(runif(n))
    D <- max(i / n - u, u - (i - 1) / n)
    events <- c(
      KS = D <= critical$KS,
      BJ = bj_statistic(u) <= critical$BJ,
      DW = dw_statistic(u) <= critical$DW
    )
    inverted <- vapply(
      bands,
      function(band) all(band$lower <= u & u <= band$upper),
      logical(1)
    )
    mismatches <- mismatches + as.integer(events != inverted)
  }

  if (any(mismatches != 0L)) {
    stop("A band inversion does not reproduce its statistic event.")
  }

  data.frame(
    method = names(mismatches),
    cases = cases,
    mismatches = unname(mismatches),
    passed = TRUE
  )
}


# ------------------------------------------------------------
# 7. Rank calibration, reproducibility, and calibration reuse.
# ------------------------------------------------------------

check_calibration_and_reuse <- function() {
  critical <- mc_rank_critical_value(1:199, alpha = 0.1)
  stopifnot(critical == 180)

  conservative <- mc_rank_critical_value(1:10, alpha = 0.05)
  stopifnot(is.infinite(conservative))

  calibration_1 <- calibrate_bands(
    n = 12,
    alpha = 0.1,
    B = 199,
    seed = 9
  )
  calibration_2 <- calibrate_bands(
    n = 12,
    alpha = 0.1,
    B = 199,
    seed = 9
  )
  stopifnot(identical(
    calibration_1$critical_values,
    calibration_2$critical_values
  ))

  result <- compare_bands_mc(
    n = 12,
    alpha = 0.1,
    M = 30,
    band_calibration_B = 199,
    seed = 10,
    bands = calibration_1,
    progress = FALSE
  )
  stopifnot(
    identical(result$band_calibration$critical_values,
              calibration_1$critical_values),
    !anyNA(result$summary),
    all(result$summary$n_replications == 30),
    "se_coverage" %in% names(result$summary)
  )

  data.frame(
    check = c(
      "exchangeable-rank calibration",
      "seeded calibration reproducibility",
      "precomputed calibration reuse"
    ),
    passed = TRUE
  )
}


# ------------------------------------------------------------
# 8. Empirical scaling diagnostic.
#
# Timings are machine-dependent. The fitted log-log slope is
# reported, and a generous threshold guards against an accidental
# return to a quadratic implementation.
# ------------------------------------------------------------

check_scaling <- function(
    sizes = c(2000L, 4000L, 8000L, 16000L),
    repetitions = 5L,
    seed = 11,
    maximum_slope = 1.6
) {
  set.seed(seed)
  elapsed <- numeric(length(sizes))

  for (k in seq_along(sizes)) {
    n <- sizes[k]
    x <- sort(rnorm(n))
    band <- make_monotone_band(n, half_width = 0.12)

    elapsed[k] <- unname(system.time(
      for (r in seq_len(repetitions)) {
        fit <- location_scale_confidence(
          x,
          band$lower,
          band$upper,
          qfun = qnorm
        )
        stopifnot(fit$feasible)
      }
    )["elapsed"]) / repetitions
  }

  elapsed_for_fit <- pmax(elapsed, .Machine$double.eps)
  scaling_slope <- unname(coef(
    lm(log(elapsed_for_fit) ~ log(sizes))
  )[2L])
  passed <- is.finite(scaling_slope) && scaling_slope < maximum_slope
  if (!passed) {
    stop("Empirical scaling is inconsistent with the linear-time implementation.")
  }

  data.frame(
    n = sizes,
    seconds_per_call = elapsed,
    fitted_log_log_slope = scaling_slope,
    passed = passed
  )
}


run_all_checks <- function() {
  load_reference_functions()

  envelope <- check_envelope_evaluation()
  exactness <- check_against_pairwise()
  geometry <- check_geometry()
  equivariance <- check_affine_equivariance()
  edge_cases <- check_edge_cases()
  inversion <- check_band_inversion()
  calibration <- check_calibration_and_reuse()
  scaling <- check_scaling()

  print(envelope, row.names = FALSE)
  print(exactness, row.names = FALSE)
  print(geometry, row.names = FALSE)
  print(equivariance, row.names = FALSE)
  print(edge_cases, row.names = FALSE)
  print(inversion, row.names = FALSE)
  print(calibration, row.names = FALSE)
  print(scaling, row.names = FALSE)

  invisible(list(
    envelope = envelope,
    exactness = exactness,
    geometry = geometry,
    equivariance = equivariance,
    edge_cases = edge_cases,
    inversion = inversion,
    calibration = calibration,
    scaling = scaling
  ))
}


if (sys.nframe() == 0L) {
  run_all_checks()
}
