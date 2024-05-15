# Analyze the distribution of mortality rates in a cohort of newborns

# Init ------------------------------------------------------------

library(tidyverse)

source('src/00-global.R')

path <- list(
  lifetab = 'out/lifetab.rds',
  hazards = 'out/hazards.rds',
  out = 'out',
  distribution = 'out/distribution.rds'
)

distribution <- list()
fig <- list()

# Load data -------------------------------------------------------

# neonatal and infant death counts and exposures by stratum
lifetab <- readRDS(path$lifetab)
# estimated neonatal hazard trajectories by stratum
hazards <- readRDS(path$hazards)

# Mean-mode convergence -------------------------------------------

# plot age specific mortality of reference group vs.
# age specific population mortality and add the age
# specific ratios
distribution$mean_mode_ratio <-
  lifetab$neonatal_j %>%
  # add mortality of modal stratum
  # Regular or high bw, Term or post-term, Apgar 9
  left_join(
    filter(lifetab$neonatal_jk, k == lifetab$neonatal_modal_k$k) %>%
      select(j, m_jr = m_jk)
  ) %>%
  mutate(
    mean_mode_ratio =
      m_j/m_jr
  )

fig$memora <-
  distribution$mean_mode_ratio %>%
  ggplot(aes(x = x_j)) +
  geom_segment(
    aes(
      xend = x_j,
      y = m_jr,
      yend = m_j
    ),
    size = 0.4
  ) +
  geom_point(
    aes(y = m_jr),
    size = 0.6
  ) +
  geom_point(
    aes(y = m_j),
    size = 0.6
  ) +
  geom_point(
    aes(
      y = exp(0.5*log(m_jr*m_j)),
    ),
    color = 'white',
    size = 5
  ) +
  geom_text(
    aes(
      y = exp(0.5*log(m_jr*m_j)),
      label = formatC(mean_mode_ratio, format = 'd')
    ),
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = c(0, 7, 14, 21, 28)
  ) +
  scale_y_continuous(
    trans = 'log10',
    breaks = figspec$bl$hazard1,
    labels = figspec$bl$hazard1
  ) +
  coord_cartesian(
    ylim = c(6*10^(-6), 10^(-2))
  ) +
  figspec$MyGGplotTheme(size = 7, ar = 0.5, axis = 'xy') +
  labs(
    x = 'Days since birth',
    y = 'Deaths per person-day'
  )

# Taylor's law ----------------------------------------------------

# Taylor's law
distribution$taylor <-
  lifetab$neonatal_jk %>%
  group_by(j, x_jk) %>%
  summarise(
    m_j = sum(m_jk*pi_jk),
    var_j = sum((m_jk-m_j)^2*pi_jk)
  ) %>%
  rename(x_j = x_jk) %>%
  ungroup()

distribution$taylor_regression <-
  lm(log(var_j) ~ log(m_j), distribution$taylor)

# plot-taylors-law
fig$taylorslaw <-
  hazards$cohort_survival_x %>%
  ggplot() +
  # loglog-linear weighted least squares fit
  geom_abline(
    slope = distribution$taylor_regression$coefficients[2],
    intercept = distribution$taylor_regression$coefficients[1],
    color = 'grey80', size = 0.5
  ) +
  annotate(
    'text',
    x = -8, y = -18,
    size = 3.5/2,
    label =
      paste0(
        "log*Var*bgroup('[', h[k](x), ']')==",
        formatC(
          distribution$taylor_regression$coefficients[1],
          digits = 1, format = 'f'
        ),
        formatC(
          distribution$taylor_regression$coefficients[2],
          digits = 1, format = 'f', flag = '+'
        ),
        '%*%',
        "log*E*bgroup('[', h[k](x), ']')"
      ),
    parse = TRUE
  ) +
  # life-table death rates and variances
  geom_point(
    aes(x = log(m_j), y = log(var_j)),
    size = 0.2,
    data = distribution$taylor
  ) +
  # predicted from random effects NB regression
  geom_line(
    aes(
      x = log(h_x),
      y = log(disp_h_xk*h_x)
    ),
    size = 0.3
  ) +
  geom_text(
    aes(
      x = log(m_j)+0.1,
      y = log(var_j)-0.2,
      label = x_j
    ),
    size = 3.5/2,
    hjust = 0,
    data =
      distribution$taylor %>%
      filter(x_j %in% c(0:5, 10, 20))
  ) +
  scale_y_continuous(
    breaks = log(10^seq(-8, -3, 1)),
    labels =
      parse(
        text =
          c('10^-8', 'scriptstyle(-7)',
            'scriptstyle(-6)',
            'scriptstyle(-5)',
            'scriptstyle(-4)',
            '10^-3')
      )
  ) +
  scale_x_continuous(
    breaks =
      log(10^seq(-5, -2, 1)),
    labels =
      parse(
        text =
          c('10^-5', 'phantom()^-4', 'phantom()^-3', '10^-2')
      ),
    limits = c(log(10^(-5)), log(10^(-2)))
  ) +
  figspec$MyGGplotTheme(ar = 0.9, size = 7, axis = 'xy', grid = 'xy') +
  labs(
    x = 'Average mortality across strata',
    y = 'Variance of mortality across strata'
  )

fig$taylorslaw

# Hazard quantiles --------------------------------------------------------

# age specific quantiles of group specific hazards
distribution$hazard_quantiles_xk <-
  hazards$cohort_survival_xk %>%
  group_by(x) %>%
  arrange(x, h_xk) %>%
  mutate(
    # cumulative population proportion
    cum_p = c(0, head(cumsum(pi_kx), -1)),
    # cumulative hazard
    cum_h = h_xk
  ) %>%
  ungroup() %>%
  select(
    x, k,
    birthweight_c4,
    gestation_at_delivery_c4,
    apgar5,
    cum_p,
    cum_h
  )

# age specific quantiles of group specific mortality rates
distribution$deathrate_quantiles_jk <-
  lifetab$neonatal_jk %>%
  group_by(x_jk) %>%
  arrange(x_jk, m_jk) %>%
  mutate(
    cum_p = c(0, head(cumsum(pi_jk), -1)),
    cum_h = m_jk
  ) %>%
  ungroup() %>%
  select(
    x_jk, k,
    birthweight_c4,
    gestation_at_delivery_c4,
    apgar5,
    cum_p,
    cum_h
  )

# which quantile is the average hazard
distribution$quantile_of_mean_hazard_j <-
  left_join(
    distribution$hazard_quantiles_xk,
    hazards$cohort_survival_x
  ) %>%
  group_by(x) %>%
  arrange(x, cum_h) %>%
  select(x, cum_h, h_x, cum_p) %>%
  filter(cum_h >= h_x) %>%
  dplyr::slice(1) %>%
  ungroup()

# Cumulative distribution of mortality rates
# across strata at days 0 and 27 following birth.
# The black lines correspond to the model predictions and the
# grey lines to the life-table estimates.
# Annotated is the quantile of average mortality at both ages
# as an indicator for the skewness of the distribution.
fig$quantiles <-
  distribution$hazard_quantiles_xk %>%
  filter(x %in% c(0, 27)) %>%
  ggplot() +
  geom_step(
    aes(x = cum_h, y = cum_p, group = x_jk),
    data =
      distribution$deathrate_quantiles_jk %>%
      filter(x_jk %in% c(0, 27)),
    size = 0.3,
    color = 'grey80'
  ) +
  geom_step(
    aes(
      x = cum_h, y = cum_p,
      size = as.factor(x),
      group = x
    ),
    show.legend = FALSE
  ) +
  annotate(
    'text',
    x = 0.1,
    y = 0.95,
    label = 'Day of birth',
    size = 3.5/2
  ) +
  annotate(
    'text',
    x = 0.001,
    y = 1.04,
    label = 'Day 28',
    fontface = 'bold',
    size = 3.5/2
  ) +
  geom_segment(
    aes(
      x = cum_h, xend = cum_h,
      y = 0, yend = cum_p
    ),
    linetype = 3,
    size = 0.2,
    data =
      distribution$quantile_of_mean_hazard_j %>%
      filter(x %in% c(0, 27))
  ) +
  geom_segment(
    aes(
      x = cum_h, xend = 10^(-6),
      y = cum_p, yend = cum_p
    ),
    linetype = 3,
    size = 0.2,
    data =
      distribution$quantile_of_mean_hazard_j %>%
      filter(x %in% c(0, 27))
  ) +
  geom_text(
    aes(
      x = 10^(-6)*1.2, y = cum_p*1.02,
      label = formatC(cum_p, digits = 2),
      hjust = 0
    ),
    size = 3.5/2,
    data =
      distribution$quantile_of_mean_hazard_j %>%
      filter(x %in% c(0, 27))
  ) +
  scale_size_discrete(range = c(0.3, 0.5)) +
  scale_x_log10(
    breaks =
      c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels =
      parse(
        text =
          c('10^-6',
            'phantom()^-5', 'phantom()^-4',
            'phantom()^-3', 'phantom()^-2', 'phantom()^-1',
            '10^0')
      ),
    limits =
      c(1e-6, 2),
    expand =
      c(0,0)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c('0', '.25', '.5', '.75', '1'),
    limits = c(0, 1.1),
    expand = c(0,0)
  ) +
  labs(
    x = 'Mortality/hazard rate',
    y = 'Cumulative population proportion'
  ) +
  figspec$MyGGplotTheme(ar = 0.9, size = 7, axis = 'xy', grid = '')

# Density of death rates over age ---------------------------------

distribution$deathrate_density <- list()

# calculate mortality rate densities on an age*mortality grid
FindIntervalStart <-
  function (x, breaks) {
    breaks[.bincode(x = x, breaks = breaks,
                    # [a, b)
                    right = FALSE, include.lowest = FALSE)]
  }
distribution$deathrate_density <-
  within(distribution$deathrate_density, {
    age_breaks = c(0, 7, 14, 21, 27)
    bins = 100
    y_breaks_start =
      10^(seq(log10(1e-6), log10(4), length.out = bins))
    y_breaks_width =
      10^(diff(log10(y_breaks_start)))[1]
    x_breaks_start =
      seq(0, 30, by = 1)
    x_breaks_width =
      diff(x_breaks_start)[1]
  })

# estimate density of death rates over time by binning death rates
# with probability of bin given by relative exposure contributed by
# groups in that bin
distribution$deathrate_density$neonatal_density_j <-
  lifetab$neonatal_jk %>%
  mutate(
    m_start =
      FindIntervalStart(m_jk, distribution$deathrate_density$y_breaks_start),
    x_start =
      FindIntervalStart(x_jk, distribution$deathrate_density$x_breaks_start)
  ) %>%
  group_by(m_start, x_start) %>%
  summarise(
    D_jm = sum(D_jk),
    O_jm = sum(O_jk)
  ) %>%
  ungroup() %>%
  complete(
    x_start,
    m_start,
    fill = list(O_jm = 0, D_jm = 0)
  ) %>%
  group_by(x_start) %>%
  mutate(
    p_jm = O_jm/sum(O_jm)
  ) %>%
  ungroup() %>%
  mutate(
    m_width =
      distribution$deathrate_density$y_breaks_width*m_start-m_start,
    x_width =
      distribution$deathrate_density$x_breaks_width,
    y  = x_start-p_jm*distribution$deathrate_density$bins
  )

distribution$deathrate_density$neonatal_density_j <-
  distribution$deathrate_density$neonatal_density_j %>%
  mutate(m_start = c(m_start[-1], NA)) %>%
  bind_rows(., distribution$deathrate_density$neonatal_density_j) %>%
  arrange(x_start, m_start)

fig$densities <-
  distribution$deathrate_density$neonatal_density_j %>%
  filter(x_start %in% distribution$deathrate_density$age_breaks) %>%
  ggplot() +
  geom_ribbon(
    aes(x = m_start, ymin = x_start, ymax = x_start-log1p(p_jm*300),
        group = x_start),
    fill = 'grey80', color = NA
  ) +
  geom_path(
    aes(x = m_j, y = x_j),
    size = 2, color = 'white',
    data = lifetab$neonatal_j
  ) +
  geom_path(
    aes(x = m_j, y = x_j),
    data = lifetab$neonatal_j
  ) +
  geom_point(
    aes(x = m_j, y = x_j),
    data = lifetab$neonatal_j
  ) +
  scale_x_log10() +
  scale_y_continuous(
    breaks = distribution$deathrate_density$age_breaks
  ) +
  labs(
    x = 'Deaths per person-hour',
    y = 'Days since birth'
  ) +
  coord_flip(ylim = c(-3, 27), clip = 'off') +
  figspec$MyGGplotTheme()

# Counterfactual hazards ------------------------------------------

fig$counterfactual <-
  hazards$cohort_survival_x %>%
  ggplot(aes(x = x)) +
  geom_line(
    aes(x = x, y = h_x),
    lwd = 0.5, col = 'black'
  ) +
  geom_line(
    aes(x = x, y = h_tilde_x),
    lwd = 0.3, col = 'black'
  ) +
  annotate(
    'text',
    x = 20, y = 2.7*1e-5,
    label = 'With selection',
    fontface = 'bold',
    size = 3
  ) +
  annotate(
    'text',
    x = 20, y = 4.8*1e-5,
    label = 'Without selection',
    size = 3
  ) +
  scale_x_continuous(
    breaks = figspec$bl$age
  ) +
  scale_y_log10(
    breaks = figspec$bl$hazard2,
    labels = parse(text = names(figspec$bl$hazard2)),
    limit = c(1e-5, 2.5e-3),
    expand = c(0,0)
  ) +
  figspec$MyGGplotTheme(ar = 0.8, size = 7) +
  labs(
    x = 'Days since birth',
    y = 'Deaths per person-day'
  )

# the ratio of counterfactual vs. estimated hazard at x = 27 is
hazards$cohort_survival_x %>%
  filter(x == 27) %>%
  transmute(h_tilde_x/h_x)

# Export ----------------------------------------------------------

saveRDS(distribution, file = path$distribution)

figspec$ExportFigure(
  fig$memora, path = path$out,
  filename = 'memora',
  device = 'pdf'
)

figspec$ExportFigure(
  fig$taylorslaw, path = path$out,
  filename = 'taylorslaw',
  device = 'pdf'
)

figspec$ExportFigure(
  fig$quantiles, path = path$out,
  filename = 'quantiles',
  device = 'pdf'
)

figspec$ExportFigure(
  fig$densities, path = path$out,
  filename = 'densities',
  device = 'pdf'
)

figspec$ExportFigure(
  fig$counterfactual, path = path$out,
  filename = 'counterfactual',
  device = 'pdf'
)
