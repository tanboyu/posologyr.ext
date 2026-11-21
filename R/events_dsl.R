#' Construct a posologyr Event Block
#'
#' Create a clean dosing or observation event block with optional repeated doses
#' (`ii`, `addl`) and infusion duration (`dur`) or rate (`rate`).
#'
#' @param time Numeric vector of event times (default `0`).
#' @param amt Numeric dose amount (default `0` for observation rows).
#' @param dur Optional infusion duration in time units.
#' @param rate Optional infusion rate in amount/time units.
#' @param ii Inter-dose interval when `addl > 0`.
#' @param addl Number of additional doses after each event in `time` (integer >= 0).
#' @param cmt Compartment index or name (default `1`).
#' @param evid Event ID (`1` = dose, `0` = observation).
#' @param dv Optional observed value (`NA_real_` by default).
#' @param is_tdm Logical indicating whether rows represent actual TDM observations.
#' @param id Subject ID (default `1L`).
#' @param covariates Named list or single-row `data.frame` of covariates to attach.
#' @return An object of class `"poso_ev_block"` (a `data.frame`).
#' @export
poso_ev <- function(time = 0,
                    amt = 0,
                    dur = NA_real_,
                    rate = NA_real_,
                    ii = NA_real_,
                    addl = 0L,
                    cmt = 1L,
                    evid = if (any(amt > 0)) 1L else 0L,
                    dv = NA_real_,
                    is_tdm = FALSE,
                    id = 1L,
                    covariates = NULL) {
  stopifnot(is.numeric(time), length(time) >= 1L)
  addl <- as.integer(addl)
  evid <- as.integer(evid)

  df <- data.frame(
    ID = as.integer(id),
    TIME = as.numeric(time),
    DV = as.numeric(dv),
    AMT = as.numeric(amt),
    DUR = as.numeric(dur),
    RATE = as.numeric(rate),
    II = as.numeric(ii),
    ADDL = addl,
    CMT = cmt,
    EVID = evid,
    IS_TDM = as.logical(is_tdm),
    stringsAsFactors = FALSE
  )

  if (!is.null(covariates)) {
    cov_df <- as.data.frame(covariates, stringsAsFactors = FALSE)
    if (nrow(cov_df) == 1L && nrow(df) > 1L) {
      cov_df <- cov_df[rep(1L, nrow(df)), , drop = FALSE]
    } else if (nrow(cov_df) != nrow(df)) {
      stop("`covariates` must have 1 row or the same number of rows as events.", call. = FALSE)
    }
    for (nm in names(cov_df)) {
      df[[nm]] <- cov_df[[nm]]
    }
  }

  class(df) <- unique(c("poso_ev_block", "data.frame"))
  df
}

#' Expand ADDL/II Dosing Rows into Explicit Time Rows
#'
#' @param df A `data.frame` containing `TIME`, `II`, `ADDL`, `EVID`.
#' @return A `data.frame` with `ADDL` expanded into explicit rows and `ADDL = 0L`.
#' @noRd
.expand_addl_rows <- function(df) {
  if (!all(c("TIME", "II", "ADDL") %in% names(df))) return(df)
  out_list <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    row_i <- df[i, , drop = FALSE]
    n_add <- if (!is.na(row_i$ADDL) && row_i$ADDL > 0L) as.integer(row_i$ADDL) else 0L
    if (n_add > 0L && !is.na(row_i$II) && row_i$II > 0) {
      times_all <- row_i$TIME + seq.int(0L, n_add) * row_i$II
      expanded <- row_i[rep(1L, length(times_all)), , drop = FALSE]
      expanded$TIME <- times_all
      expanded$ADDL <- 0L
      out_list[[i]] <- expanded
    } else {
      row_i$ADDL <- 0L
      out_list[[i]] <- row_i
    }
  }
  do.call(rbind, out_list)
}

#' Propagate Covariates Across Time-Sorted Rows
#'
#' @param df Time-sorted `data.frame`.
#' @param cov_cols Character vector of covariate column names.
#' @param method `"locf"` (last observation carried forward) or `"nocb"`.
#' @return Updated `data.frame`.
#' @noRd
.fill_covariates <- function(df, cov_cols, method = c("locf", "nocb")) {
  method <- match.arg(method)
  if (length(cov_cols) == 0L || nrow(df) <= 1L) return(df)
  for (col in cov_cols) {
    v <- df[[col]]
    if (method == "locf") {
      for (i in seq.int(2L, length(v))) {
        if (is.na(v[i]) && !is.na(v[i - 1L])) v[i] <- v[i - 1L]
      }
      non_na <- which(!is.na(v))
      if (length(non_na) > 0L && non_na[1L] > 1L) {
        v[seq_len(non_na[1L] - 1L)] <- v[non_na[1L]]
      }
    } else {
      for (i in rev(seq_len(length(v) - 1L))) {
        if (is.na(v[i]) && !is.na(v[i + 1L])) v[i] <- v[i + 1L]
      }
    }
    df[[col]] <- v
  }
  df
}

#' Sequence Two Event Blocks in Time (`mrgsolve`-style `%then%`)
#'
#' Shifts `e2` forward by the effective end time of `e1` and combines them.
#'
#' @param e1 First event block or `data.frame`.
#' @param e2 Second event block or `data.frame` starting at relative time `0`.
#' @return Combined `"poso_ev_block"` sorted by `ID` and `TIME`.
#' @export
`%then_poso%` <- function(e1, e2) {
  d1 <- .expand_addl_rows(as.data.frame(e1, stringsAsFactors = FALSE))
  d2 <- .expand_addl_rows(as.data.frame(e2, stringsAsFactors = FALSE))

  end_t1 <- max(d1$TIME, na.rm = TRUE)
  if ("II" %in% names(d1) && any(!is.na(d1$II))) {
    last_dose_idx <- which.max(d1$TIME)
    if (!is.na(d1$II[last_dose_idx]) && d1$II[last_dose_idx] > 0) {
      end_t1 <- max(end_t1, d1$TIME[last_dose_idx] + d1$II[last_dose_idx])
    }
  }

  d2$TIME <- d2$TIME + end_t1
  common_cols <- union(names(d1), names(d2))
  for (nm in setdiff(common_cols, names(d1))) d1[[nm]] <- NA
  for (nm in setdiff(common_cols, names(d2))) d2[[nm]] <- NA

  combined <- rbind(d1[, common_cols, drop = FALSE], d2[, common_cols, drop = FALSE])
  combined <- combined[order(combined$ID, combined$TIME, -combined$EVID), , drop = FALSE]
  rownames(combined) <- NULL
  class(combined) <- unique(c("poso_ev_block", "data.frame"))
  combined
}

#' Add Observation Time Points or TDM Measurements to an Event Table
#'
#' Quickly append dense prediction grid times (`evid = 0, dv = NA`) or actual
#' TDM sampling points (`evid = 0, dv = ...`) to an existing event table or
#' `data.frame`.
#'
#' @param x Event block or `data.frame`.
#' @param time Numeric vector of observation times.
#' @param dv Optional numeric vector of observed concentrations/responses.
#' @param cmt Compartment index or name (default `1L`).
#' @param is_tdm Logical indicating whether these points are real TDM observations.
#' @param cov_fill Covariate fill strategy (`"locf"` or `"nocb"`).
#' @return A `"poso_ev_block"` object sorted by `ID` and `TIME`.
#' @export
poso_add_obs <- function(x,
                         time,
                         dv = NA_real_,
                         cmt = 1L,
                         is_tdm = !all(is.na(dv)),
                         cov_fill = c("locf", "nocb")) {
  cov_fill <- match.arg(cov_fill)
  base_df <- as.data.frame(x, stringsAsFactors = FALSE)
  id_val <- if ("ID" %in% names(base_df)) base_df$ID[1L] else 1L

  if (length(dv) == 1L && length(time) > 1L) {
    dv_vec <- rep(as.numeric(dv), length(time))
  } else if (length(dv) == length(time)) {
    dv_vec <- as.numeric(dv)
  } else {
    stop("`dv` must have length 1 or match the length of `time`.", call. = FALSE)
  }

  obs_df <- data.frame(
    ID = as.integer(id_val),
    TIME = as.numeric(time),
    DV = dv_vec,
    AMT = 0,
    DUR = NA_real_,
    RATE = NA_real_,
    II = NA_real_,
    ADDL = 0L,
    CMT = cmt,
    EVID = 0L,
    IS_TDM = as.logical(is_tdm),
    stringsAsFactors = FALSE
  )

  all_cols <- union(names(base_df), names(obs_df))
  for (nm in setdiff(all_cols, names(base_df))) base_df[[nm]] <- NA
  for (nm in setdiff(all_cols, names(obs_df))) obs_df[[nm]] <- NA

  out <- rbind(base_df[, all_cols, drop = FALSE], obs_df[, all_cols, drop = FALSE])
  out <- out[order(out$ID, out$TIME, -out$EVID), , drop = FALSE]

  core_cols <- c("ID", "TIME", "DV", "AMT", "DUR", "RATE", "II", "ADDL", "CMT", "EVID", "IS_TDM")
  cov_cols <- setdiff(names(out), core_cols)
  out <- .fill_covariates(out, cov_cols = cov_cols, method = cov_fill)

  rownames(out) <- NULL
  class(out) <- unique(c("poso_ev_block", "data.frame"))
  out
}

#' Convert Event Table to posologyr / rxode2 Compatible Data Frame
#'
#' Expands repeated doses (`ADDL`/`II`), normalizes column types, and drops
#' helper metadata columns if desired.
#'
#' @param x Event block or `data.frame`.
#' @param keep_meta Keep `IS_TDM` column (default `TRUE`).
#' @return A plain `data.frame` ready for `posologyr` / `rxode2`.
#' @export
as_poso_data <- function(x, keep_meta = TRUE) {
  df <- .expand_addl_rows(as.data.frame(x, stringsAsFactors = FALSE))
  df <- df[order(df$ID, df$TIME, -df$EVID), , drop = FALSE]
  if (!keep_meta && "IS_TDM" %in% names(df)) {
    df$IS_TDM <- NULL
  }
  for (col in c("RATE", "II")) {
    if (col %in% names(df) && all(is.na(df[[col]]))) {
      df[[col]] <- NULL
    }
  }
  rownames(df) <- NULL
  df
}
