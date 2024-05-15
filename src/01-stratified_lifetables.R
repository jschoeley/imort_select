# Prepare stratified neonatal life-tables and calculate summaries
#
# Jonas Schöley

# Init ------------------------------------------------------------

library(tidyverse)

cnst <- list(
  # birthweight and gestation labels
  recode_bweight =
    c('Extremely low' = 'Extremely low [0,1000)g',
      'Very low' = 'Very low [1000,1500)g',
      'Low' = 'Low [1500,2500)g',
      'Regular or high' = 'Regular [2500,4200)g',
      'Regular or high' = 'High 4200g+',
      '(Missing)' = '(Missing)'),
  recode_gestation =
    c('Extremely preterm' = 'Extremely preterm <28w',
      'Very preterm' = 'Very preterm [28,32)w',
      'Moderate to late preterm' = 'Moderate to late preterm [32,37)w',
      'Term or post-term' = 'Term 37w+',
      '(Missing)' = '(Missing)')
)

path <- list(
  usinfants = 'dat/usinfants.rds',
  lifetab = 'out/lifetab.rds'
)

dat <- list()
lifetab <- list()

# Functions -------------------------------------------------------

# Get Stratified Life-table From Survival Data
#
# Aggregate individual level survival data into a stratified
# life-tables with prespecified age-groups using the assumptions of
# piecewise-constant hazards.
#
# @df a data frame
# @x time at censoring or event
# @event event indicator, 0 if censoring, 1 if event
# @cuts vector of cuts for life-table age groups
# @... strata
GetLifeTable <- function (df, x, event, cuts, ...) {

  x_i = enquo(x); event_i = enquo(event); strata = enexprs(...)

  lt <-
    df %>%
    # aggregate individual level event and censoring times
    # into predefined age groups
    mutate(
      j =
        .bincode(
          x = !!x_i, breaks = cuts,
          # [a, b)
          right = FALSE, include.lowest = FALSE
        ),
      x_jk =
        cuts[j]
    ) %>%
    group_by(!!!strata, j, x_jk) %>%
    summarise(
      # total observed deaths in interval
      D_jk = sum(!!event_i),
      # total observed censorings in interval
      C_jk = sum(!(!!event_i)),
      # average time spent in interval for
      # those who leave during interval (by death or censoring)
      a_jk = mean(!!x_i - first(x_jk))
    ) %>%
    ungroup() %>%
    # add predefined age intervals to table where no deaths or
    # censorings occoured
    complete(
      nesting(x_jk = head(cuts, -1), j = 1:length(head(cuts, -1))),
      nesting(!!!strata),
      fill = list(D_jk = 0, C_jk = 0, a_jk = NA)
    ) %>%
    arrange(!!!strata, x_jk) %>%
    group_by(!!!strata) %>%
    mutate(k = cur_group_id()) %>%
    group_by(k) %>%
    # calculate life-table
    mutate(
      # width of the interval
      n_jk =
        c(diff(x_jk), last(cuts)-last(x_jk)),
      # observed population alive at start of interval
      N_jk =
        head(cumsum(c(sum(C_jk+D_jk), -(C_jk+D_jk))), -1),
      # total person-time of exposure over interval
      O_jk =
        (N_jk-D_jk-C_jk)*n_jk +
        ifelse(is.na(a_jk), 0, a_jk*(D_jk+C_jk)),
      # mortality rate over interval
      m_jk =
        ifelse(O_jk == 0, 0, D_jk/O_jk),
      # probability of surviving interval
      # assuming constant force of mortality
      p_jk =
        exp(-m_jk*n_jk),
      # probability of death during the interval given survival
      # to interval
      q_jk =
        1-p_jk,
      # probability of survival until interval start
      l_jk =
        cumprod(c(1, head(p_jk, -1))),
      # probability of death in interval
      d_jk =
        c(-diff(l_jk), last(l_jk)*last(q_jk))
    ) %>%
    # distribution of stratum-specific exposures by age
    group_by(j) %>%
    mutate(
      pi_jk = O_jk/sum(O_jk)
    ) %>%
    ungroup() %>%
    select(
      k, j, !!!strata, x_jk, n_jk, N_jk, a_jk, O_jk, D_jk, C_jk,
      m_jk, q_jk, l_jk, d_jk, pi_jk
    )

  return(lt)

}

DescribeSurvivalByStratum <- function (df, death, ...) {

  death_i = enquo(death); strata = enexprs(...)

  df %>%
    group_by(!!!strata) %>%
    summarise(
      births = n(),
      deaths = sum(!!death_i),
      survivors = births-deaths,
      prb_surv = 1 - deaths/births
    ) %>%
    mutate(
      prp_births = births/sum(births),
      prp_deaths = deaths/sum(deaths),
      diff_prp_births = survivors/sum(survivors) - prp_births
    ) %>%
    select(
      -survivors
    )

}

# Load data -------------------------------------------------------

# individual level data on infant deaths to birth cohorts US 1995-2012
dat$usinfants <- readRDS(path$usinfants)

# Survival data preparation ---------------------------------------

# basic data pooling, selection and recoding
dat$usinfants <- dat$usinfants %>%
  # subset to study period
  filter(date_of_delivery_y %in% 2008:2012) %>%
  # add variables for survival analysis
  mutate(
    # age at death in fractional days
    age_at_death_fd =
      case_when(
        age_at_death_c == 'Under 1 hour' ~ 0,
        age_at_death_c == '1-23 hours' ~ 1/24,
        TRUE ~ as.numeric(age_at_death_d)
      ),
    # width of observation interval
    age_at_death_fd_width =
      case_when(
        age_at_death_c == 'Under 1 hour' ~ 1/24,
        age_at_death_c == '1-23 hours' ~ 23/24,
        !(is.na(age_at_death_d)) ~ 1
      ),
    # infant death indicator
    infant_death =
      !is.na(age_at_death_d),
    # infant survival time in fractional days
    # observations right-censored at day 365
    infant_survtime_fd =
      if_else(infant_death, age_at_death_fd, 365),
    infant_survtime_fd_width =
      if_else(infant_death, age_at_death_fd_width, 0),
    # neonatal death indicator
    neonatal_death =
      !is.na(age_at_death_d) & age_at_death_d < 28,
    # neonatal survival time in fractional days
    # observations right-censored at day 28
    neonatal_survtime_fd =
      if_else(neonatal_death, age_at_death_fd, 28),
    neonatal_survtime_fd_width =
      if_else(neonatal_death, age_at_death_fd_width, 0),
    neonatal_survtime_d =
      if_else(neonatal_death, age_at_death_d, 28L),
    neonatal_survtime_d_width =
      if_else(neonatal_death, 1, 0)
  ) %>%
  # make NAs explicit, recode, relabel and reorder levels
  mutate(
    birthweight_c4 =
      fct_explicit_na(birthweight_c5) %>%
      # merge regular and high birthweight groups
      fct_recode(!!!cnst$recode_bweight) %>%
      fct_relevel(c('Extremely low', 'Very low', 'Low',
                    'Regular or high', '(Missing)')),
    gestation_at_delivery_c4 =
      fct_explicit_na(gestation_at_delivery_c4) %>%
      fct_recode(!!!cnst$recode_gestation) %>%
      fct_relevel(c('Extremely preterm', 'Very preterm',
                    'Moderate to late preterm',
                    'Term or post-term',
                    '(Missing)')),
    apgar5 =
      fct_explicit_na(as.factor(apgar5)) %>%
      fct_relevel(as.character(0:10))
  ) %>%
  select(
    # survival
    infant_death, infant_survtime_fd, infant_survtime_fd_width,
    neonatal_death, neonatal_survtime_fd, neonatal_survtime_fd_width,
    neonatal_survtime_d, neonatal_survtime_d_width,
    # strata
    birthweight_c4, gestation_at_delivery_c4, apgar5
  )

# Create stratified life-tables -----------------------------------

# create life-tables stratified by the intersections of
# birthweight,
# gestation at delivery,
# 5 minute Apgar score

# stratified neonatal life-tables
lifetab$neonatal_jk <-
  dat$usinfants %>%
  # infant life-tables with daily age-grouping
  # stratified by birthweight, age of gestation and
  # five minute apgar score
  GetLifeTable(
    x = neonatal_survtime_fd + 0.5*neonatal_survtime_fd_width,
    event = neonatal_death, cuts = c(0:29),
    birthweight_c4, gestation_at_delivery_c4, apgar5
  ) %>%
  filter(x_jk != 28)

# neonatal life-table w/o 0 exposure for log-linear regression
lifetab$neonatal_regression_jk <-
  lifetab$neonatal_jk %>% filter(O_jk != 0)

# Create population life-tables -----------------------------------

# population neonatal life-table
lifetab$neonatal_j <-
  lifetab$neonatal_jk %>%
  group_by(j, x_jk) %>%
  summarise(
    D_j = sum(D_jk), N_j = sum(N_jk),
    O_j = sum(O_jk), m_j = D_j/O_j
  ) %>%
  rename(x_j = x_jk) %>%
  ungroup()

# population infant life-table
lifetab$infant_j <-
  dat$usinfants %>%
  GetLifeTable(
    x = infant_survtime_fd + 0.5*infant_survtime_fd_width,
    event = infant_death, cuts = 0:366,
  ) %>%
  filter(x_jk != 366, !is.na(O_jk)
  ) %>%
  select(j, x_j = x_jk, D_j = D_jk, N_j = N_jk, O_j = O_jk, m_j = m_jk)

# Generate life-table summaries -----------------------------------

# information about the strata in the neonatal life-table
lifetab$neonatal_k <-
  lifetab$neonatal_jk %>%
  group_by(
    k, birthweight_c4, gestation_at_delivery_c4, apgar5
  ) %>%
  summarise(
    # share of stratum k on total population at birth
    pi_k = pi_jk[j == 1]
  ) %>%
  ungroup()

# most prevalent group
lifetab$neonatal_modal_k <-
  lifetab$neonatal_k %>% filter(pi_k == max(pi_k))

# total number of births, deaths and
# person-days under observation
lifetab$birth_death_exposure <-
  lifetab$neonatal_j %>%
  summarise(
    D = sum(D_j),
    N = first(N_j),
    O = sum(O_j)
  )

# summary table of neonatal population considered in this analysis
lifetab$summary <-
  dat$usinfants %>%
  {
    bind_rows(
      'Total' =
        DescribeSurvivalByStratum(., neonatal_death, NA) %>%
        rename(level = `NA`),
      '5 minute Apgar' =
        DescribeSurvivalByStratum(., neonatal_death, apgar5) %>%
        rename(level = apgar5),
      'Gestation at delivery' =
        DescribeSurvivalByStratum(., neonatal_death, gestation_at_delivery_c4) %>%
        rename(level = gestation_at_delivery_c4),
      'Birthweight' =
        DescribeSurvivalByStratum(., neonatal_death, birthweight_c4) %>%
        rename(level = birthweight_c4),
      .id = 'variable'
    )
  }

# Exports ---------------------------------------------------------

saveRDS(lifetab, path$lifetab)
