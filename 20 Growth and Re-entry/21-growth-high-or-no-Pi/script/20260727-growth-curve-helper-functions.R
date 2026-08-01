# Helper functions for the 2026-07-27 growth-curve analysis

# Supplementary display helpers. Data import, data shaping, and the
# calculations used to obtain growth statistics remain in the main Rmd.

plot_growth_curves <- function(inoculum) {
  ggplot2::ggplot(
    dplyr::filter(analysis_data, Inoculate == inoculum),
    ggplot2::aes(
      x = Elapsed, y = Measurements, color = Genotype,
      group = interaction(Date, Well)
    )
  ) +
    ggplot2::geom_line(alpha = 0.25, linewidth = 0.45) +
    ggplot2::geom_line(
      data = dplyr::filter(growth_means, Inoculate == inoculum),
      ggplot2::aes(
        x = Elapsed, y = mean_OD, color = Genotype, group = Genotype
      ),
      linewidth = 1.25
    ) +
    ggplot2::facet_wrap(
      ~Treatment, nrow = 1,
      labeller = ggplot2::labeller(Treatment = treatment_labels)
    ) +
    ggplot2::scale_color_manual(values = species_colors, drop = FALSE) +
    ggplot2::labs(
      title = paste("Inoculum:", inoculum, "phase"),
      x = "Elapsed time (hours)",
      y = "Background-corrected OD600",
      color = "Species"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(face = 3)
    )
}

mean_sd <- function(x) {
  tibble::tibble(
    y = mean(x, na.rm = TRUE),
    ymin = mean(x, na.rm = TRUE) - stats::sd(x, na.rm = TRUE),
    ymax = mean(x, na.rm = TRUE) + stats::sd(x, na.rm = TRUE)
  )
}

plot_metric_summary <- function(data, inoculum, title) {
  ggplot2::ggplot(
    dplyr::filter(data, Inoculate == inoculum),
    ggplot2::aes(x = Genotype, y = Value, color = Genotype)
  ) +
    ggplot2::geom_jitter(
      width = 0.10, height = 0, alpha = 0.65, size = 1.8
    ) +
    ggplot2::stat_summary(
      fun.data = mean_sd, geom = "errorbar", width = 0.18
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "point", shape = 21,
      fill = "white", size = 2.8, stroke = 0.9
    ) +
    ggplot2::facet_grid(
      Metric ~ Treatment, scales = "free_y",
      labeller = ggplot2::labeller(Treatment = treatment_labels)
    ) +
    ggplot2::scale_color_manual(values = species_colors, drop = FALSE) +
    ggplot2::labs(
      title = paste(title, "-", inoculum, "inoculum"),
      x = NULL, y = NULL, color = "Species"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "none"
    )
}

plot_growthcurver_fits <- function(inoculum) {
  data <- dplyr::filter(growthcurver_fitted, Inoculate == inoculum)

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Elapsed, color = Genotype,
      group = interaction(Date, Well)
    )
  ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Measurements), alpha = 0.35, linewidth = 0.4
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = fitted_OD), linetype = 2, linewidth = 0.75
    ) +
    ggplot2::facet_grid(
      Treatment ~ Genotype,
      labeller = ggplot2::labeller(Treatment = treatment_labels)
    ) +
    ggplot2::scale_color_manual(values = species_colors, drop = FALSE) +
    ggplot2::labs(
      title = paste("growthcurver fits:", inoculum, "inoculum"),
      subtitle = "Observed curves are solid; logistic fits are dashed",
      x = "Elapsed time (hours)",
      y = "Background-corrected OD600"
    ) +
    ggplot2::theme(legend.position = "none")
}

plot_gcplyr_sensitivity <- function(inoculum) {
  ggplot2::ggplot(
    dplyr::filter(gcplyr_sensitivity, Inoculate == inoculum),
    ggplot2::aes(
      x = rate_13, y = rate_loess, color = Genotype, shape = Treatment
    )
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey70") +
    ggplot2::geom_point(size = 2.3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = species_colors, drop = FALSE) +
    ggplot2::labs(
      title = paste("gcplyr rate sensitivity:", inoculum, "inoculum"),
      x = "13-point maximum growth rate",
      y = "Drop-two + LOESS + 9-point maximum growth rate",
      color = "Species", shape = "Treatment"
    ) +
    ggplot2::theme(legend.position = "bottom")
}

plot_method_comparison <- function(inoculum) {
  ggplot2::ggplot(
    dplyr::filter(method_comparison_long, Inoculate == inoculum),
    ggplot2::aes(
      x = gcplyr, y = growthcurver, color = Genotype, shape = Treatment
    )
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey70") +
    ggplot2::geom_point(size = 2.2, alpha = 0.8) +
    ggplot2::facet_wrap(~Metric, scales = "free", nrow = 1) +
    ggplot2::scale_color_manual(values = species_colors, drop = FALSE) +
    ggplot2::labs(
      title = paste("Method comparison:", inoculum, "inoculum"),
      x = "gcplyr estimate", y = "growthcurver estimate",
      color = "Species", shape = "Treatment"
    ) +
    ggplot2::theme(legend.position = "bottom")
}

#' Diagnose gcplyr lag-time or maximum-growth-rate estimates
#'
#' Recalculates the per-capita derivative for selected curves, identifies its
#' maximum, and plots either the lag-time tangent construction or the maximum
#' per-capita growth rate and its corresponding minimum doubling time.
#'
#' @param df A data frame containing Date, Well, Elapsed, and Measurements.
#'   Measurements must already be background-corrected.
#' @param date One or more dates coercible with as.Date().
#' @param wells One or more plate-well names, such as c("B2", "C2").
#' @param diagnostic Either "lag_time" or "max_growth_rate".
#' @param baseline Either "minimum" for the global minimum observed density or
#'   "first_minimum" for the first local minimum detected by gcplyr. Used only
#'   for the lag-time diagnostic.
#' @param window_width_n Number of observations used by gcplyr::calc_deriv().
#'
#' @return A ggplot object. The per-curve estimates used in the plot are
#'   available with attr(plot, "diagnostic_summary").
plot_gcplyr_diagnostic <- function(
    df,
    date,
    wells,
    diagnostic = c("lag_time", "max_growth_rate"),
    baseline = c("minimum", "first_minimum"),
    window_width_n = 5
) {
  diagnostic <- match.arg(diagnostic)
  baseline <- match.arg(baseline)
  date <- as.Date(date)
  wells <- as.character(wells)

  required_cols <- c("Date", "Well", "Elapsed", "Measurements")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "df is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  diagnostic_data <- df |>
    dplyr::filter(
      .data$Date %in% date,
      as.character(.data$Well) %in% wells
    ) |>
    dplyr::group_by(.data$Date, .data$Well) |>
    dplyr::arrange(.data$Elapsed, .by_group = TRUE) |>
    dplyr::mutate(
      diagnostic_percap = gcplyr::calc_deriv(
        x = .data$Elapsed,
        y = .data$Measurements,
        percapita = TRUE,
        blank = 0,
        window_width_n = window_width_n,
        trans_y = "log",
        warn_ungrouped = FALSE
      )
    ) |>
    dplyr::ungroup()

  if (nrow(diagnostic_data) == 0) {
    stop("No curves matched the requested date and wells.")
  }

  if (any(diagnostic_data$Measurements <= 0, na.rm = TRUE)) {
    stop(
      "The selected curves contain corrected OD values at or below zero. ",
      "The log-scale lag diagnostic is undefined for those values."
    )
  }

  diagnostic_summary <- diagnostic_data |>
    dplyr::group_by(.data$Date, .data$Well) |>
    dplyr::summarize(
      baseline_dens = if (baseline == "minimum") {
        gcplyr::min_gc(.data$Measurements)
      } else {
        gcplyr::first_minima(
          x = .data$Elapsed,
          y = .data$Measurements,
          return = "y"
        )
      },
      max_percap = gcplyr::max_gc(.data$diagnostic_percap),
      max_percap_index = gcplyr::which_max_gc(.data$diagnostic_percap),
      max_percap_time = .data$Elapsed[.data$max_percap_index],
      max_percap_dens = .data$Measurements[.data$max_percap_index],
      min_doubling_time = gcplyr::doubling_time(y = .data$max_percap),
      lag_time = suppressWarnings(
        gcplyr::lag_time(
          x = .data$Elapsed,
          y = .data$Measurements,
          deriv = .data$diagnostic_percap,
          blank = 0,
          y0 = .data$baseline_dens
        )
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      tangent_intercept =
        log(.data$max_percap_dens) -
        .data$max_percap * .data$max_percap_time,
      lag_label = sprintf(
        "lag = %.2f h\nmax rate = %.2f h\u207B\u00B9",
        .data$lag_time,
        .data$max_percap
      ),
      max_growth_label = sprintf(
        "max rate = %.2f h\u207B\u00B9\nat %.2f h\nmin doubling time = %.2f h",
        .data$max_percap,
        .data$max_percap_time,
        .data$min_doubling_time
      )
    )

  diagnostic_theme <- ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 14, face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  if (diagnostic == "lag_time") {
    p <- ggplot2::ggplot(
      diagnostic_data,
      ggplot2::aes(x = .data$Elapsed, y = log(.data$Measurements))
    ) +
      ggplot2::geom_point(size = 0.9, alpha = 0.7) +
      ggplot2::geom_line(linewidth = 0.45, color = "grey40") +
      ggplot2::geom_abline(
        data = diagnostic_summary,
        ggplot2::aes(
          slope = .data$max_percap,
          intercept = .data$tangent_intercept
        ),
        inherit.aes = FALSE,
        color = "#D55E00",
        linewidth = 0.8
      ) +
      ggplot2::geom_point(
        data = diagnostic_summary,
        ggplot2::aes(
          x = .data$max_percap_time,
          y = log(.data$max_percap_dens)
        ),
        inherit.aes = FALSE,
        color = "#D55E00",
        size = 2.2
      ) +
      ggplot2::geom_vline(
        data = diagnostic_summary,
        ggplot2::aes(xintercept = .data$lag_time),
        color = "#0072B2",
        linetype = 2,
        linewidth = 0.7
      ) +
      ggplot2::geom_hline(
        data = diagnostic_summary,
        ggplot2::aes(yintercept = log(.data$baseline_dens)),
        color = "#0072B2",
        linetype = 3,
        linewidth = 0.7
      ) +
      ggplot2::geom_label(
        data = diagnostic_summary,
        ggplot2::aes(
          x = Inf,
          y = Inf,
          label = .data$lag_label
        ),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.1,
        size = 3,
        linewidth = 0,
        fill = scales::alpha("white", 0.8)
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$Date, .data$Well),
        scales = "free_y",
        labeller = ggplot2::label_both
      ) +
      ggplot2::labs(
        title = "gcplyr lag-time diagnostic",
        subtitle = paste(
          "Baseline:",
          ifelse(
            baseline == "minimum",
            "minimum observed OD",
            "first local minimum"
          )
        ),
        x = "Elapsed time (hours)",
        y = "log(background-corrected OD600)",
        caption = paste(
          "Orange: maximum-growth tangent and originating point.",
          "Blue: estimated lag time and baseline density."
        )
      ) +
      diagnostic_theme
  } else {
    p <- ggplot2::ggplot(
      diagnostic_data,
      ggplot2::aes(x = .data$Elapsed, y = .data$diagnostic_percap)
    ) +
      ggplot2::geom_line(linewidth = 0.55, color = "grey40") +
      ggplot2::geom_point(
        data = diagnostic_summary,
        ggplot2::aes(
          x = .data$max_percap_time,
          y = .data$max_percap
        ),
        inherit.aes = FALSE,
        color = "#D55E00",
        size = 2.5
      ) +
      ggplot2::geom_label(
        data = diagnostic_summary,
        ggplot2::aes(
          x = Inf,
          y = Inf,
          label = .data$max_growth_label
        ),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.1,
        size = 3,
        linewidth = 0,
        fill = scales::alpha("white", 0.8)
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$Date, .data$Well),
        scales = "free_y",
        labeller = ggplot2::label_both
      ) +
      ggplot2::coord_cartesian(ylim = c(-1, NA)) +
      ggplot2::labs(
        title = "gcplyr maximum-growth-rate diagnostic",
        subtitle = paste(
          "Per-capita derivative calculated with a",
          window_width_n,
          "point window"
        ),
        x = "Elapsed time (hours)",
        y = "Per-capita growth rate (h\u207B\u00B9)",
        caption = paste(
          "The orange point marks the maximum per-capita growth rate.",
          "Minimum doubling time is log(2) divided by that rate."
        )
      ) +
      diagnostic_theme
  }

  attr(p, "diagnostic_summary") <- diagnostic_summary
  if (diagnostic == "lag_time") {
    attr(p, "lag_summary") <- diagnostic_summary
  } else {
    attr(p, "max_growth_summary") <- diagnostic_summary
  }
  p
}

#' Plot the gcplyr lag-time diagnostic
#'
#' Backward-compatible convenience wrapper for plot_gcplyr_diagnostic().
plot_gcplyr_lag_diagnostic <- function(
    df,
    date,
    wells,
    baseline = c("minimum", "first_minimum"),
    window_width_n = 5
) {
  plot_gcplyr_diagnostic(
    df = df,
    date = date,
    wells = wells,
    diagnostic = "lag_time",
    baseline = baseline,
    window_width_n = window_width_n
  )
}

#' Plot the gcplyr maximum-growth-rate diagnostic
plot_gcplyr_max_growth_diagnostic <- function(
    df,
    date,
    wells,
    window_width_n = 5
) {
  plot_gcplyr_diagnostic(
    df = df,
    date = date,
    wells = wells,
    diagnostic = "max_growth_rate",
    window_width_n = window_width_n
  )
}
