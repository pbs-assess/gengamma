plot_aic <- function(dat) {
  bp +
    geom_tile(
    data = dat,
    mapping = aes(
      x = 1, y = fspecies2,
      width = Inf, height = 1, fill = factor(odd_species2)
    )
  ) +
    geom_point(
      data = dat, #ll_diff |> filter(all_ok) |> filter(region == "GOA"), # models that converged
      aes(x = aic_w, y = fspecies2, colour = family, shape = family),
      stroke = 0.75, size = 1.75, position = ggstance::position_dodgev(height = 0.8)
    )  +
    scale_colour_manual(values = family_colours) +
    scale_shape_manual(values = family_shapes) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = scales::label_percent(suffix = ""),
      limits = c(0, 1)) +
      # limits = c(0, 1.27)) +
    labs(x = parse(text = "AIC~weight~(`%`)"), y = "Species",
      shape = "Family",
      colour = "Family") +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.placement = "outside",
          strip.text.y = element_text(face = "bold"))
}

plot_cv <- function(dat, centre_min = TRUE) {
  scaler <- 1
  xlabels <- c("> 400", "250", "0")
  if (!centre_min) {
    scaler <- -1
    xlabels <- c("> -400", "-250", "0")
  }
  bp +
    geom_tile(
      data = dat,
      mapping = aes(
        x = scaler, y = fspecies2,
        width = Inf, height = 1, fill = factor(odd_species2)
    )) +
    geom_point(
      data = dat,
      aes(x = scaler * pmin(diff_ll, 500), y = fspecies2, colour = family, shape = family),
      stroke = 0.75, size = 1.75, position = ggstance::position_dodgev(height = 0.8)
    )  +
    scale_colour_manual(values = family_colours) +
    scale_shape_manual(values = family_shapes) +
    scale_x_continuous(breaks = scaler * c(400, 200, 0),
      labels = c("> -400", "-250", "0"),
      limits = sort(scaler * c(500, 0))) +
    labs(x = "abs(sumLL diff)", y = "Species",
      shape = "Family",
      colour = "Family") +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.text.y = element_blank())
}

plot_q <- function(dat) {
  bp +
    geom_tile(
      data = dat,
      mapping = aes(
        x = 1, y = fspecies2,
        width = Inf, height = 1, fill = factor(odd_species2)
    )) +
    geom_text(data = filter(dat, family == "delta-gengamma"),
    aes(x = 1, y = fspecies2, label = round(est_q, digits = 1)),
      hjust = 1, size = 3, colour = "grey50") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()#,
          # axis.title.x = element_text()
          ) +
    scale_x_continuous(limits =c(0.3, 1.25)) +
    coord_cartesian(expand = FALSE) +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.text.y = element_blank())
}

plot_n <- function(dat) {
  bp +
    geom_tile(
      data = dat,
      mapping = aes(
        x = 1, y = fspecies2,
        width = Inf, height = 1, fill = factor(odd_species2)
    )) +
    geom_text(data = filter(dat, family == "delta-gengamma"),
      # aes(x = 1, y = fspecies2, label = paste0("(", count, ")")),
      aes(x = 1, y = fspecies2, label = format(round(rate, 2), nsmall = 2)),
      hjust = 1, size = 3, colour = "grey50") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "white"), #element_blank(),
          axis.ticks.x = element_line(colour = "white"), # element_blank()) +
          axis.title.x = element_blank()
          ) +
    scale_x_continuous(limits =c(0.3, 1.25)) +
    coord_cartesian(expand = FALSE) +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.text.y = element_blank())
}

no_x <- function() {
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
  )
}

invis_x_axis <- function() {
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = "white"),
    axis.ticks.x = element_line(colour = "white")
  )
}

get_q_rank <- function(dat, region) {
  # browser()
  dat |>
    filter(family == "delta-gengamma") |>
    select(region, species, est_q) |>
    arrange(est_q) |>
    mutate(q_rank = row_number()) |>
  left_join(dat |> select(-est_q)) |>
  mutate(fspecies = forcats::fct_rev(species)) |>
  mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(q_rank))) |>
  mutate(odd_species2 = ifelse(as.numeric(fspecies2) %% 2 == 0, 1, 0))
}

get_legend <- function(plot) {
  gtable <- ggplotGrob(plot)
  legend <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}
