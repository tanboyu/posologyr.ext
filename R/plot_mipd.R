#' Clinical MIPD ggplot2 Theme
#'
#' @param base_size Base font size (default `12`).
#' @return A `ggplot2` theme object.
#' @export
theme_poso_mipd <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#EAECEF", linewidth = 0.4),
      strip.background = ggplot2::element_rect(fill = "#F4F6F9", color = "#D0D7DE"),
      strip.text = ggplot2::element_text(face = "bold", color = "#1F2328"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      legend.position = "top",
      legend.box.margin = ggplot2::margin(0, 0, -4, 0)
    )
}

#' Plot a Tidy Simulation / Posterior Object (`PKPDsim` & `mrgsolve` Style)
#'
#' @param x A `"poso_sim"` object.
#' @param vars Optional character vector of variables to plot.
#' @param show_ribbon Logical; plot outer uncertainty band (default `TRUE`).
#' @param show_indiv Logical; plot individual trajectories (default `FALSE`).
#' @param show_obs Logical; overlay observed TDM/data points (default `TRUE`).
#' @param show_doses Logical; mark dosing times along the bottom axis (default `TRUE`).
#' @param target_range Optional numeric vector of length 2 (`c(lower, upper)`) to draw
#'   a shaded therapeutic window.
#' @param log_y Logical; use log10 scale on y-axis (default `FALSE`).
#' @param ncol Number of facet columns when plotting multiple variables.
#' @param title Optional plot title.
#' @param xlab X-axis label (default `"Time"`).
#' @param ylab Y-axis label (default `"Concentration / Response"`).
#' @param ... Additional arguments passed to `theme_poso_mipd()`.
#' @return A `ggplot` object.
#' @export
plot.poso_sim <- function(x,
                          vars = NULL,
                          show_ribbon = TRUE,
                          show_indiv = FALSE,
                          show_obs = TRUE,
                          show_doses = TRUE,
                          target_range = NULL,
                          log_y = FALSE,
                          ncol = NULL,
                          title = NULL,
                          xlab = "Time",
                          ylab = "Concentration / Response",
                          ...) {
  stopifnot(inherits(x, "poso_sim"))
  target_vars <- if (is.null(vars)) x$vars else intersect(vars, x$vars)
  if (length(target_vars) == 0L) {
    stop("None of the requested `vars` found in simulation object.", call. = FALSE)
  }

  sum_df <- x$summary[x$summary$variable %in% target_vars, , drop = FALSE]
  traj_df <- x$trajectories[x$trajectories$variable %in% target_vars, , drop = FALSE]

  p <- ggplot2::ggplot(sum_df, ggplot2::aes(x = .data$time))

  if (!is.null(target_range) && length(target_range) == 2L) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = min(target_range),
      ymax = max(target_range),
      fill = "#2DA44E",
      alpha = 0.14
    ) +
      ggplot2::geom_hline(
        yintercept = target_range,
        linetype = "dashed",
        color = "#1A7F37",
        linewidth = 0.5
      )
  }

  if (isTRUE(show_indiv)) {
    p <- p + ggplot2::geom_line(
      data = traj_df,
      ggplot2::aes(y = .data$value, group = interaction(.data$id, .data$variable)),
      color = "#6E7781",
      alpha = 0.25,
      linewidth = 0.35
    )
  }

  q_cols <- grep("^q[0-9]{2}$", names(sum_df), value = TRUE)
  if (isTRUE(show_ribbon) && length(q_cols) >= 2L) {
    low_col <- q_cols[1L]
    high_col <- q_cols[length(q_cols)]
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[low_col]], ymax = .data[[high_col]]),
      fill = "#0969DA",
      alpha = 0.18
    )
  }

  med_col <- if ("q50" %in% names(sum_df)) "q50" else "mean"
  p <- p + ggplot2::geom_line(
    ggplot2::aes(y = .data[[med_col]]),
    color = "#0969DA",
    linewidth = 1.0
  )

  if (isTRUE(show_obs) && !is.null(x$obs) && nrow(x$obs) > 0L) {
    obs_sub <- x$obs
    if (!all(is.na(obs_sub$variable))) {
      obs_sub <- obs_sub[obs_sub$variable %in% target_vars | is.na(obs_sub$variable), , drop = FALSE]
      obs_sub$variable[is.na(obs_sub$variable)] <- target_vars[1L]
    } else {
      obs_sub$variable <- target_vars[1L]
    }
    if (nrow(obs_sub) > 0L) {
      p <- p + ggplot2::geom_point(
        data = obs_sub,
        ggplot2::aes(x = .data$time, y = .data$dv),
        color = "#CF222E",
        fill = "#CF222E",
        size = 2.6,
        shape = 21
      )
    }
  }

  if (isTRUE(show_doses) && !is.null(x$doses) && nrow(x$doses) > 0L) {
    p <- p + ggplot2::geom_rug(
      data = x$doses,
      ggplot2::aes(x = .data$time),
      sides = "b",
      color = "#8250DF",
      linewidth = 0.8,
      inherit.aes = FALSE
    )
  }

  if (length(target_vars) > 1L) {
    p <- p + ggplot2::facet_wrap(~variable, scales = "free_y", ncol = ncol)
  }

  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    theme_poso_mipd(...)
}

#' Plot Individual MAP Trajectory Overlaid with Observed TDM Data
#'
#' @param x A `"poso_map_tidy"` object containing `indiv_df` (`time`, `ipred`, optional `pred`)
#'   and optional `obs_df` (`time`, `dv`).
#' @param show_pred Logical; overlay population prediction curve if available (default `TRUE`).
#' @param pred_label Legend label for the second curve (default `"Population (PRED)"`).
#' @param ipred_label Legend label for the primary curve (default `"Individual MAP (IPRED)"`).
#' @param target_range Optional numeric vector of length 2 (`c(lower, upper)`).
#' @param log_y Logical; use log10 y-axis (default `FALSE`).
#' @param title Optional plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param ... Passed to `theme_poso_mipd()`.
#' @return A `ggplot` object.
#' @export
plot.poso_map_tidy <- function(x,
                               show_pred = TRUE,
                               pred_label = "Population (PRED)",
                               ipred_label = "Individual MAP (IPRED)",
                               target_range = NULL,
                               log_y = FALSE,
                               title = "Individual MAP Profile",
                               xlab = "Time",
                               ylab = "Concentration",
                               ...) {
  stopifnot(inherits(x, "poso_map_tidy"))
  df <- x$indiv_df

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time))

  if (!is.null(target_range) && length(target_range) == 2L) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = min(target_range),
      ymax = max(target_range),
      fill = "#2DA44E",
      alpha = 0.14
    ) +
      ggplot2::geom_hline(
        yintercept = target_range,
        linetype = "dashed",
        color = "#1A7F37",
        linewidth = 0.5
      )
  }

  if (isTRUE(show_pred) && "pred" %in% names(df)) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(y = .data$pred, linetype = pred_label),
      color = "#57606A",
      linewidth = 0.8
    )
  }

  p <- p + ggplot2::geom_line(
    ggplot2::aes(y = .data$ipred, linetype = ipred_label),
    color = "#0969DA",
    linewidth = 1.05
  )

  if (!is.null(x$obs_df) && nrow(x$obs_df) > 0L) {
    p <- p + ggplot2::geom_point(
      data = x$obs_df,
      ggplot2::aes(x = .data$time, y = .data$dv),
      color = "#CF222E",
      size = 2.8
    )
  }

  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }

  lt_values <- stats::setNames(c("solid", "dashed"), c(ipred_label, pred_label))

  p +
    ggplot2::scale_linetype_manual(
      name = NULL,
      values = lt_values
    ) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    theme_poso_mipd(...)
}

#' Build a Customizable Multi-Plot MIPD Bundle
#'
#' Generates a named list of `ggplot` objects (e.g. `prior`, `map`, `forecast`)
#' so users can either print all plots by default or select specific plot types.
#'
#' @param prior_sim Optional `"poso_sim"` object for prior simulation.
#' @param map_tidy Optional `"poso_map_tidy"` object for individual MAP fit.
#' @param forecast_tidy Optional `"poso_map_tidy"` object comparing adjusted vs unadjusted regimen.
#' @param include Character vector of plot names to include (`c("prior", "map", "forecast")` by default).
#' @param target_range Optional numeric vector of length 2 (`c(lower, upper)`).
#' @param log_y Logical; use log10 scale on y-axis (default `FALSE`).
#' @return An object of class `"poso_plot_bundle"`.
#' @export
poso_plot_bundle <- function(prior_sim = NULL,
                             map_tidy = NULL,
                             forecast_tidy = NULL,
                             include = c("prior", "map", "forecast"),
                             target_range = NULL,
                             log_y = FALSE) {
  include <- match.arg(include, choices = c("prior", "map", "forecast"), several.ok = TRUE)
  plots <- list()

  if ("prior" %in% include && !is.null(prior_sim)) {
    plots$prior <- plot(
      prior_sim,
      show_ribbon = TRUE,
      show_obs = TRUE,
      show_doses = TRUE,
      target_range = target_range,
      log_y = log_y,
      title = "1. Prior Population Distribution vs. Observed TDM"
    )
  }

  if ("map" %in% include && !is.null(map_tidy)) {
    plots$map <- plot(
      map_tidy,
      show_pred = TRUE,
      target_range = target_range,
      log_y = log_y,
      title = "2. Individual MAP Fit (IPRED) vs. Population (PRED)"
    )
  }

  if ("forecast" %in% include && !is.null(forecast_tidy)) {
    plots$forecast <- plot(
      forecast_tidy,
      show_pred = TRUE,
      pred_label = "Unadjusted Regimen",
      ipred_label = "Adjusted Regimen",
      target_range = target_range,
      log_y = log_y,
      title = "3. Regimen Forecast: Adjusted vs. Unadjusted"
    )
  }

  structure(list(plots = plots), class = "poso_plot_bundle")
}

#' Plot a poso_plot_bundle Object (Combined if patchwork Available, Else Printed/Returned)
#'
#' @param x A `"poso_plot_bundle"` object.
#' @param which Optional character vector of plot keys to draw.
#' @param ncol Number of columns when combining multiple plots (default `1L`).
#' @param ... Unused.
#' @return A combined `patchwork` object if installed, otherwise the list of `ggplot` objects.
#' @export
plot.poso_plot_bundle <- function(x, which = NULL, ncol = 1L, ...) {
  stopifnot(inherits(x, "poso_plot_bundle"))
  pl <- x$plots
  if (!is.null(which)) {
    pl <- pl[intersect(names(pl), which)]
  }
  if (length(pl) == 0L) {
    stop("No matching plots found in bundle.", call. = FALSE)
  }

  if (requireNamespace("patchwork", quietly = TRUE) && length(pl) > 1L) {
    combined <- patchwork::wrap_plots(pl, ncol = ncol)
    print(combined)
    return(invisible(combined))
  }

  for (nm in names(pl)) {
    print(pl[[nm]])
  }
  invisible(pl)
}
