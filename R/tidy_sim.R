#' Configure Multi-Threading for rxode2 Simulation
#'
#' @param n_cores Number of threads to use (`NULL` defaults to `min(detectCores(), 4L)`).
#' @return Integer number of threads applied.
#' @noRd
.configure_threads <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    det <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 2L)
    if (is.na(det) || det < 1L) det <- 2L
    n_cores <- as.integer(min(det, 4L))
  } else {
    n_cores <- max(1L, as.integer(n_cores))
  }

  tryCatch(
    {
      if (exists("setRxThreads", where = asNamespace("rxode2"), mode = "function")) {
        rxode2::setRxThreads(n_cores)
      }
    },
    error = function(e) NULL
  )
  n_cores
}

#' Extract Variable Column Case-Insensitively from a Data Frame
#'
#' @param df A `data.frame`.
#' @param var_name Target variable name.
#' @return Numeric vector or `NULL` if not found.
#' @noRd
.extract_var_ci <- function(df, var_name) {
  idx <- match(tolower(var_name), tolower(names(df)))
  if (is.na(idx)) return(NULL)
  as.numeric(df[[idx]])
}

#' Build a Tidy Simulation Object from Simulated Trajectories
#'
#' @param sim_df `data.frame` containing `id`, `time`, and one or more simulated
#'   state/output columns (e.g., `Cc`, `Cp`, `effect`).
#' @param vars Optional character vector of output column names to extract. If
#'   `NULL`, automatically selects numeric columns other than `id`, `time`, `evid`.
#' @param obs_df Optional `data.frame` of observations (`time`, `variable`/`cmt`, `dv`, `is_tdm`).
#' @param dose_df Optional `data.frame` of dosing events (`time`, `amt`, `dur`).
#' @param probs Numeric vector of quantiles for uncertainty ribbons (default `c(0.05, 0.50, 0.95)`).
#' @return An object of class `"poso_sim"`.
#' @export
as_poso_sim <- function(sim_df,
                        vars = NULL,
                        obs_df = NULL,
                        dose_df = NULL,
                        probs = c(0.05, 0.50, 0.95)) {
  stopifnot(is.data.frame(sim_df))
  orig_names <- names(sim_df)
  names_low <- tolower(orig_names)
  names(sim_df) <- names_low

  if (!"time" %in% names(sim_df)) {
    stop("`sim_df` must contain a `time` column.", call. = FALSE)
  }
  if (!"id" %in% names(sim_df)) {
    sim_df$id <- 1L
  }

  n_rows <- nrow(sim_df)

  if (is.null(vars)) {
    exclude_cols <- c("id", "time", "evid", "amt", "dv", "dur", "rate", "ii", "addl", "cmt", "is_tdm")
    candidate_cols <- setdiff(names(sim_df), exclude_cols)
    matched_idx <- which(vapply(sim_df[candidate_cols], is.numeric, logical(1)))
    vars_low <- candidate_cols[matched_idx]
    vars_label <- orig_names[match(vars_low, names_low)]
  } else {
    vars_low <- tolower(vars)
    missing_vars <- vars[!vars_low %in% names(sim_df)]
    if (length(missing_vars) > 0L) {
      stop(
        sprintf(
          "Requested `vars` not found in simulation output: %s. Available columns: %s",
          paste(missing_vars, collapse = ", "),
          paste(orig_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    vars_label <- vars
  }

  stopifnot(length(vars_low) >= 1L)

  long_list <- lapply(seq_along(vars_low), function(i) {
    v_low <- vars_low[i]
    v_lbl <- vars_label[i]
    vals <- as.numeric(sim_df[[v_low]])
    if (length(vals) != n_rows) {
      stop(
        sprintf("Column `%s` has length %d, expected %d.", v_lbl, length(vals), n_rows),
        call. = FALSE
      )
    }
    data.frame(
      id = sim_df$id,
      time = sim_df$time,
      variable = v_lbl,
      value = vals,
      stringsAsFactors = FALSE
    )
  })
  long_df <- do.call(rbind, long_list)

  probs <- sort(unique(probs))
  grp_keys <- paste(long_df$variable, long_df$time, sep = "||")
  split_idx <- split(seq_len(nrow(long_df)), grp_keys)

  summary_rows <- lapply(split_idx, function(idx) {
    sub_df <- long_df[idx, , drop = FALSE]
    qs <- stats::quantile(sub_df$value, probs = probs, na.rm = TRUE, names = FALSE)
    out <- data.frame(
      variable = sub_df$variable[1L],
      time = sub_df$time[1L],
      mean = mean(sub_df$value, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (k in seq_along(probs)) {
      q_label <- paste0("q", sprintf("%02d", round(probs[k] * 100)))
      out[[q_label]] <- qs[k]
    }
    out
  })
  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(summary_df$variable, summary_df$time), , drop = FALSE]
  rownames(summary_df) <- NULL

  res <- list(
    trajectories = long_df,
    summary = summary_df,
    obs = obs_df,
    doses = dose_df,
    probs = probs,
    vars = vars_label
  )
  class(res) <- "poso_sim"
  res
}

#' Run and Tidy a posologyr Population Simulation with Multi-Threading by Default
#'
#' @param model An `rxode2` model function or object.
#' @param events Event table (`poso_ev_block` or `data.frame`).
#' @param n_sim Number of simulated subjects/replicates (default `200L`).
#' @param vars Optional character vector of variables to summarize.
#' @param seed Optional integer random seed.
#' @param n_cores Number of threads for `rxode2` (default `NULL` -> up to 4 cores).
#' @return A `"poso_sim"` object.
#' @export
poso_sim_tidy <- function(model,
                          events,
                          n_sim = 200L,
                          vars = NULL,
                          seed = NULL,
                          n_cores = NULL) {
  .configure_threads(n_cores = n_cores)
  if (!is.null(seed)) set.seed(as.integer(seed))
  ev_df <- as_poso_data(events, keep_meta = TRUE)

  obs_df <- NULL
  if ("DV" %in% names(ev_df)) {
    has_dv <- !is.na(ev_df$DV) & ev_df$EVID == 0L
    if (any(has_dv)) {
      obs_df <- data.frame(
        id = ev_df$ID[has_dv],
        time = ev_df$TIME[has_dv],
        dv = ev_df$DV[has_dv],
        variable = if (is.null(vars)) NA_character_ else vars[1L],
        is_tdm = if ("IS_TDM" %in% names(ev_df)) as.logical(ev_df$IS_TDM[has_dv]) else TRUE,
        stringsAsFactors = FALSE
      )
    }
  }

  dose_df <- NULL
  has_dose <- ev_df$EVID == 1L & !is.na(ev_df$AMT) & ev_df$AMT > 0
  if (any(has_dose)) {
    dose_df <- data.frame(
      id = ev_df$ID[has_dose],
      time = ev_df$TIME[has_dose],
      amt = ev_df$AMT[has_dose],
      dur = if ("DUR" %in% names(ev_df)) ev_df$DUR[has_dose] else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  single_ev <- ev_df[ev_df$ID == ev_df$ID[1L], , drop = FALSE]
  pop_ev <- do.call(rbind, lapply(seq_len(as.integer(n_sim)), function(i) {
    tmp <- single_ev
    tmp$ID <- i
    tmp
  }))
  pop_ev_clean <- as_poso_data(pop_ev, keep_meta = FALSE)

  sim_out <- tryCatch(
    {
      s <- posologyr::poso_simu_pop(dat = pop_ev_clean, prior_model = model, n_pop = as.integer(n_sim))
      if (is.data.frame(s)) {
        s
      } else if (is.list(s)) {
        df_cand <- NULL
        for (nm in c("simu", "sim", "out", "model", "data")) {
          if (nm %in% names(s) && (is.data.frame(s[[nm]]) || inherits(s[[nm]], "rxSolve"))) {
            df_cand <- as.data.frame(s[[nm]])
            break
          }
        }
        if (is.null(df_cand)) {
          for (el in s) {
            if (is.data.frame(el) || inherits(el, "rxSolve")) {
              df_cand <- as.data.frame(el)
              break
            }
          }
        }
        if (is.null(df_cand)) stop("No data.frame found in poso_simu_pop output.", call. = FALSE)
        df_cand
      } else {
        as.data.frame(s)
      }
    },
    error = function(e) {
      rx_obj <- if (is.function(model)) rxode2::rxode2(model) else model
      as.data.frame(rxode2::rxSolve(rx_obj, pop_ev_clean))
    }
  )

  if (!is.null(obs_df) && all(is.na(obs_df$variable))) {
    sim_low <- tolower(names(sim_out))
    cand_idx <- match(c("cc", "cp", "c", "conc"), sim_low)
    cand_idx <- cand_idx[!is.na(cand_idx)]
    if (length(cand_idx) > 0L) obs_df$variable <- names(sim_out)[cand_idx[1L]]
  }

  as_poso_sim(
    sim_df = sim_out,
    vars = vars,
    obs_df = obs_df,
    dose_df = dose_df
  )
}

#' Safely Align Two Trajectory Data Frames by Time (Handles Unequal Row Counts)
#'
#' @param df_ipred `data.frame` with `time` and simulated variable column.
#' @param df_pred Optional `data.frame` with `time` and simulated variable column.
#' @param var_name Target variable name (default `"Cc"`).
#' @return A `data.frame` with columns `time`, `ipred`, and optional `pred`.
#' @noRd
.align_ipred_pred <- function(df_ipred, df_pred = NULL, var_name = "Cc") {
  t_ip <- as.numeric(df_ipred$time)
  y_ip <- .extract_var_ci(df_ipred, var_name)
  if (is.null(y_ip)) {
    stop(sprintf("Variable `%s` not found in IPRED solution.", var_name), call. = FALSE)
  }

  out <- data.frame(time = t_ip, ipred = y_ip, stringsAsFactors = FALSE)
  out <- out[order(out$time), , drop = FALSE]
  out <- out[!duplicated(out$time), , drop = FALSE]

  if (!is.null(df_pred)) {
    t_pr <- as.numeric(df_pred$time)
    y_pr <- .extract_var_ci(df_pred, var_name)
    if (!is.null(y_pr)) {
      ord_pr <- order(t_pr)
      t_pr <- t_pr[ord_pr]
      y_pr <- y_pr[ord_pr]
      keep_pr <- !duplicated(t_pr)
      t_pr <- t_pr[keep_pr]
      y_pr <- y_pr[keep_pr]

      if (length(t_pr) == nrow(out) && isTRUE(all.equal(t_pr, out$time))) {
        out$pred <- y_pr
      } else {
        out$pred <- stats::approx(x = t_pr, y = y_pr, xout = out$time, rule = 2)$y
      }
    }
  }
  rownames(out) <- NULL
  out
}

#' Construct a Tidy MAP Trajectory Object Safely from Solution Data Frames
#'
#' Automatically aligns `time` across IPRED and PRED solutions even when row
#' counts differ due to event time points.
#'
#' @param ipred_df `data.frame` from `rxSolve` using individual MAP parameters.
#' @param pred_df Optional `data.frame` from `rxSolve` using typical population parameters (`ETA = 0`).
#' @param var_name Variable column to extract (default `"Cc"`).
#' @param obs_df Optional `data.frame` of observations (`time`, `dv`).
#' @return An object of class `"poso_map_tidy"`.
#' @export
as_poso_map_tidy <- function(ipred_df, pred_df = NULL, var_name = "Cc", obs_df = NULL) {
  aligned <- .align_ipred_pred(df_ipred = ipred_df, df_pred = pred_df, var_name = var_name)
  res <- list(
    indiv_df = aligned,
    obs_df = obs_df,
    var_name = var_name
  )
  class(res) <- "poso_map_tidy"
  res
}

#' One-Step Individual MAP & Typical PRED Trajectory Builder
#'
#' Solves both individual MAP (`IPRED`) and typical population (`ETA = 0`, `PRED`)
#' trajectories from a `posologyr::poso_estim_map()` result without manual
#' boilerplate.
#'
#' @param fit_map Result object from `posologyr::poso_estim_map()`.
#' @param model `rxode2` model function or object.
#' @param events Event table (`poso_ev_block` or `data.frame`).
#' @param var_name Output variable name (default `"Cc"`).
#' @param n_cores Optional number of threads (`NULL` -> auto up to 4).
#' @return A `"poso_map_tidy"` object.
#' @export
poso_map_tidy <- function(fit_map,
                          model,
                          events,
                          var_name = "Cc",
                          n_cores = NULL) {
  .configure_threads(n_cores = n_cores)
  ev_raw <- as_poso_data(events, keep_meta = TRUE)

  obs_df <- NULL
  if ("DV" %in% names(ev_raw)) {
    has_dv <- !is.na(ev_raw$DV) & ev_raw$EVID == 0L
    if (any(has_dv)) {
      obs_df <- data.frame(
        time = ev_raw$TIME[has_dv],
        dv = ev_raw$DV[has_dv],
        stringsAsFactors = FALSE
      )
    }
  }

  ev_clean <- as_poso_data(events, keep_meta = FALSE)
  rx_mod <- if (is.function(model)) rxode2::rxode2(model) else model

  ipred_params <- fit_map$indiv_param[1L, , drop = FALSE]
  ipred_sol <- as.data.frame(rxode2::rxSolve(rx_mod, params = ipred_params, events = ev_clean))

  pred_params <- ipred_params
  eta_cols <- grep("^ETA_", names(pred_params), value = TRUE)
  if (length(eta_cols) > 0L) {
    pred_params[eta_cols] <- 0
  }
  pred_sol <- as.data.frame(rxode2::rxSolve(rx_mod, params = pred_params, events = ev_clean))

  as_poso_map_tidy(
    ipred_df = ipred_sol,
    pred_df = pred_sol,
    var_name = var_name,
    obs_df = obs_df
  )
}

#' Summarize a Tidy Simulation Object
#'
#' @param object A `"poso_sim"` object.
#' @param ... Unused.
#' @export
summary.poso_sim <- function(object, ...) {
  object$summary
}
