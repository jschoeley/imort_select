# Estimate stratum specific hazard trajectories
#
# Jonas Schöley

# 1. Model the stratum and age-specific neonatal hazards (h_xk) with a
#    multilevel regression model informed by death counts and exposures
# 2. Derive stratum specific and population level survival statistics
#    from modeled hazards
# 3. Plot modeled hazards against life-table death rate estimates

# Init ------------------------------------------------------------

library(tidyverse)
library(lme4)

source('src/00-global.R')

path <- list(
  lifetab = 'out/lifetab.rds',
  hazards = 'out/hazards.rds',
  out = 'out'
)

hazards <- list()
fig <- list()

# Load data -------------------------------------------------------

# neonatal and infant death counts and exposures by stratum
lifetab <- readRDS(path$lifetab)

# Model hazards by stratum ----------------------------------------

# add degree 2 polynomial basis for log1p(age) to data
hazards$neonatal_regression_jk <-
  lifetab$neonatal_regression_jk %>%
  mutate(
    log1p_x_jk = log1p(x_jk),
    log1p_x_jk_pow2 = (log1p(x_jk))^2
  )

# fit multilevel poisson regression
hazards$fit <-
  glmer(
    D_jk ~
      1 + log1p_x_jk + log1p_x_jk_pow2 +
      (1 + log1p_x_jk + log1p_x_jk_pow2 |
         gestation_at_delivery_c4) +
      (1 + log1p_x_jk + log1p_x_jk_pow2 |
         gestation_at_delivery_c4:birthweight_c4) +
      (1 + log1p_x_jk + log1p_x_jk_pow2 |
         gestation_at_delivery_c4:birthweight_c4:apgar5) +
      offset(log(O_jk)),
    family = poisson(link = 'log'),
    data = hazards$neonatal_regression_jk,
    verbose = 2, nAGQ = 0,
    control =
      glmerControl(
        optimizer = 'bobyqa',
        optCtrl   = list(maxfun = 100000)
      )
  )

# Stratum parameters ----------------------------------------------

# Recover stratum specific hazard parameters.
# The model fit to each group k is
#
# log(E(D_jk)) = b0_k + b1_k*log(x_jk+1) + b2_k*log(x_jk+1)^2 + log(O_jk),
#
# with the b_k's being the sums of group specific random effects.
# We recover the b_k's by predicting from the model with certain
# variables set to 0.
hazards$parameters_k <-
  lifetab$neonatal_k %>%
  mutate(
    b0_k =
      predict(
        hazards$fit,
        newdata = mutate(., log1p_x_jk = 0, log1p_x_jk_pow2 = 0, O_jk = 1)
      ),
    b1_k =
      predict(
        hazards$fit,
        newdata =
          mutate(., log1p_x_jk = 1, log1p_x_jk_pow2 = 0, O_jk = 1)
      ) - b0_k,
    b2_k =
      predict(
        hazards$fit,
        newdata =
          mutate(., log1p_x_jk = 0, log1p_x_jk_pow2 = 1, O_jk = 1)
      ) - b0_k
  )

# Stratum survival ------------------------------------------------

# calculate the stratum specific hazard (h_xk), cumulative
# hazard (H_xk), survival (S_xk), and proportion on all survivors (pi_kx)
# for the first 4 weeks of life

# Hazard function: degree 2 polynomial over log(age + 1)
hazards$Hzrd <-
  function (x, b0, b1, b2) {
    exp(b0 + b1*log1p(x) + b2*log1p(x)^2)
  }

# Derivative of hazard
hazards$dHzrd <-
  function (x, b0, b1, b2) {
    (x+1)^(b1-1)*exp(b0+b2*log1p(x)^2)*(b1+2*b2*log1p(x))
  }

# Cumulative hazard function
hazards$CumHzrd <-
  function (x, b0, b1, b2) {
    b2 <- as.complex(b2)
    term1 <- sqrt(pi)*exp(b0-(b1+1)^2/(4*b2))
    term2 <- 2*sqrt(b2)
    out <-
      -term1 *
      (
        RcppFaddeeva::erfi( (b1+1) / term2 ) -
          RcppFaddeeva::erfi( (b1+2*b2*log1p(x)+1) / term2 )
      ) /
      term2
    return(Re(out))
  }

# derive stratum specific survival trajectories
# from parameter estimates
hazards$cohort_survival_xk <-
  hazards$parameters_k %>%
  crossing(x = c(seq(0, 5, 0.1), 6:27)) %>%
  mutate(
    # stratum specific hazard
    h_xk = hazards$Hzrd(x, b0_k, b1_k, b2_k),
    # derivative in stratum specific hazard
    dh_xk = hazards$dHzrd(x, b0_k, b1_k, b2_k),
    # cumulative stratum specific hazard
    H_xk = hazards$CumHzrd(x, b0_k, b1_k, b2_k),
    # stratum specific probability to survive up to x
    S_xk = exp(-H_xk),
  ) %>%
  group_by(k) %>%
  mutate(
    # relative decline of hazard from age 0 to age x
    D_xk = (h_xk-h_xk[x==0])/h_xk[x==0],
  ) %>%
  group_by(x) %>%
  mutate(
    # probability of stratum k at age x
    pi_kx = S_xk*pi_k/sum(S_xk*pi_k)
  ) %>%
  ungroup() %>%
  select(x, everything())

# Population survival ---------------------------------------------

# calculate population level statistics of interest

# marginal (population level) statistics
hazards$cohort_survival_x <-
  hazards$cohort_survival_xk %>%
  group_by(x) %>%
  summarise(

    # marginal survival
    S_x = sum(S_xk*pi_k),
    # marginal hazard
    h_x = sum(h_xk*pi_kx),
    # marginal log-hazard
    lh_x = log(h_x),

    # components of the Vaupel-Zhang equality
    # Note that the share of the marginal hazard decline
    # explained by population heterogeneity is the same no matter
    # if the relative or the absolute derivative of the hazard is
    # used:
    # h'(x) = E[h'(x|k)] - Var[h(x|k)]
    # log(h(x))' = h'(x)/h(x) =
    #   E[h'(x|k)]/h(x) - Var[h(x|k)]/h(x)
    # -Var[h(x|k)]/h(x) / h'(x)/h(x) = -Var[h(x|k)]/h'(x)

    # expectation of hazard derivative across strata
    avg_dh_xk = sum(dh_xk*pi_kx),
    # bdagger
    bdagger_x = avg_dh_xk/h_x,
    # variance of hazards across strata
    var_h_xk = sum(pi_kx*(h_xk-h_x)^2),
    # dispersion of hazards across strata
    disp_h_xk = var_h_xk/h_x,
    # derivative of marginal hazard
    dh_x = avg_dh_xk - var_h_xk,
    # relative derivative of marginal hazard
    b_x = bdagger_x - disp_h_xk,
    # compositional share on hazard derivative
    v_x = -disp_h_xk/b_x*100,

    # entropy of distribution of survivors across strata
    entropy_p_kx = -sum(log(pi_kx)*pi_kx),

    # counterfactual population hazard with population proportions
    # fixed at pk0
    h_tilde_x = sum(pi_k*h_xk),
  ) %>%
  ungroup() %>%
  mutate(
    # relative drop in counterfactual population hazard at age x
    # compared to age 0
    D_tilde_x =
      exp((h_tilde_x-h_tilde_x[x==0])/h_tilde_x[x==0]),
    # relative drop in population hazard at age x
    # compared to age 0
    D_x =
      exp((h_x-h_x[x==0])/h_x[x==0])
  ) %>%
  ungroup()

# Plot survival ---------------------------------------------------

# plot stratum specific survival trajectories
hazards$cohort_survival_xk %>%
  ggplot(aes(x = x, y = S_xk)) +
  geom_line(aes(color = apgar5)) +
  facet_grid(birthweight_c4 ~ gestation_at_delivery_c4) +
  scale_color_viridis_d(option = 'A') +
  scale_x_continuous(breaks = figspec$bl$age, trans = 'log1p') +
  scale_y_continuous(trans = 'log1p') +
  figspec$MyGGplotTheme(panel_border = T)

# Info on specific strata -----------------------------------------

# mortality at birth of the low mortality strata
hazards$cohort_survival_xk %>%
  filter(x == 0) %>%
  arrange(h_xk) %>%
  select(k, x, h_xk, birthweight_c4, gestation_at_delivery_c4, apgar5) %>%
  head()

# mortality at birth of the high mortality strata
hazards$cohort_survival_xk %>%
  filter(x == 0) %>%
  arrange(h_xk) %>%
  select(k, x, h_xk, birthweight_c4, gestation_at_delivery_c4, apgar5) %>%
  tail()

# neonatal hazard decline of the low mortality strata
hazards$cohort_survival_xk %>%
  filter(k == 226, x %in% c(0, 27)) %>%
  pull(h_xk) %>%
  {(.[2]-.[1])/.[1]}

# neonatal hazard decline of the high mortality strata
hazards$cohort_survival_xk %>%
  filter(k == 1, x %in% c(0, 27)) %>%
  pull(h_xk) %>%
  {(.[2]-.[1])/.[1]}

# Plot hazards ----------------------------------------------------

hazard_plot <- list()

# don't plot missing category
hazard_plot$cohort_survival_xk <-
  hazards$cohort_survival_xk %>%
  filter(
    birthweight_c4 != '(Missing)',
    gestation_at_delivery_c4 != '(Missing)',
    apgar5 != '(Missing)'
  )
hazard_plot$lt_neonatal_jk <-
  lifetab$neonatal_jk %>%
  filter(
    birthweight_c4 != '(Missing)',
    gestation_at_delivery_c4 != '(Missing)',
    apgar5 != '(Missing)',
    D_jk != 0
  )

# plot observed versus predicted
fig$hazards <-
  hazard_plot$cohort_survival_xk %>%
  ggplot() +
  # life-table stratum specific mortality estimates
  geom_point(
    aes(x = x_jk, y = m_jk, size = D_jk),
    size = 0.1,
    shape = 16,
    color = 'grey50',
    data = hazard_plot$lt_neonatal_jk,
    show.legend = FALSE
  ) +
  # predicted stratum specific hazard
  geom_line(
    aes(
      x = x, y = h_xk,
      group = apgar5, color = apgar5
    ),
    size = 0.2
  ) +
  # life-table population mortality
  geom_point(
    aes(x = x_j, y = m_j),
    data = lifetab$neonatal_j,
    size = 0.3,
    shape = 16
  ) +
  # predicted population hazard
  geom_line(
    aes(x = x, y = h_x),
    size = 0.3,
    data = hazards$cohort_survival_x
  ) +
  scale_color_manual(
    values = viridis::magma(14)[1:11]
  ) +
  scale_x_continuous(
    breaks = figspec$bl$age
  ) +
  scale_y_log10(
    breaks = 10^(seq(-6, 0, 1)),
    labels =
      parse(
        text =
          c('10^-6', 'scriptstyle(-5)', 'scriptstyle(-4)',
            'scriptstyle(-3)', 'scriptstyle(-2)',
            'scriptstyle(-1)', '10^0')
      )
  ) +
  facet_grid(
    gestation_at_delivery_c4 ~ birthweight_c4
  ) +
  labs(
    x = 'Days since birth',
    y = 'Deaths per person-day',
    color = 'Apgar'
  ) +
  figspec$MyGGplotTheme(
    axis = '', size = 7, ar = 0.8, panel_border = TRUE, grid = ''
  )

# Export ----------------------------------------------------------

saveRDS(hazards, path$hazards)

figspec$ExportFigure(
  fig$hazards, path = path$out,
  filename = 'hazards',
  device = 'pdf'
)
