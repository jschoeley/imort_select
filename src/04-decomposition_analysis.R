# Decompose the changing distribution of death rates over age
# into selection and convergence components
#
# Jonas Schöley

# Init ------------------------------------------------------------

library(tidyverse)

source('src/00-global.R')

path <- list(
  lifetab = 'out/lifetab.rds',
  hazards = 'out/hazards.rds',
  out = 'out',
  decomp = 'out/decomp.rds'
)

decomp <- list()
fig <- list()

# Load data -------------------------------------------------------

# neonatal and infant death counts and exposures by stratum
lifetab <- readRDS(path$lifetab)
# estimated neonatal hazard trajectories by stratum
hazards <- readRDS(path$hazards)

# Vaupel-Zhang decomposition --------------------------------------

# we already performed the decomposition in a prior step
decomp$cohort_survival_x <- hazards$cohort_survival_x

# Vaupel-Zhang decomposition table
decomp$vaupel_zhang <-
  with(decomp$cohort_survival_x %>% filter(x %in% c(0, 1, 7, 14, 21, 27)),{
    tibble(
      age =
        x,
      h_x =
        formatC(h_x, format = 'e', digits = 1),
      dh_x =
        formatC(dh_x, format = 'e', digits = 1),
      direct =
        formatC(avg_dh_xk, format = 'e', digits = 1),
      compos =
        formatC(var_h_xk, format = 'e', digits = 1),
      pct =
        formatC(v_x, format = 'f', digits = 1),
    )
  })

# Kitagawa decomposition ------------------------------------------

# Kitagawa decomposition
#
# Decomposing differences in weighted averages of measures into
# differences in weights and differences in measures
#
# @g1_k measure in population 1 by group k
# @g2_k measure in population 2 by group k
# @p1_k weight in population 1 by group k
# @p2_k weight in population 2 by group k
DecomposeKitagawa <- function (g1_k, g2_k, p1_k, p2_k) {
  d_g <- 0.5*sum((p1_k+p2_k)*(g2_k-g1_k))
  d_p <- 0.5*sum((g1_k+g2_k)*(p2_k-p1_k))

  return(c(d_g = d_g, d_p = d_p))
}

{

  # decomposition of age-differences in
  # population average mortality
  decomp$kitagawa_deathrate <-
    lifetab$neonatal_jk %>%
    filter(x_jk %in% c(0, 1, 7, 14, 21, 27)) %>%
    group_by(k) %>%
    mutate(
      lead_m_jk = lead(m_jk),
      lead_pi_jk = lead(pi_jk),
      diff_m_jk = lead_m_jk - m_jk,
      diff_pi_jk = lead_pi_jk - pi_jk
    ) %>%
    group_by(j, x_jk) %>%
    do({
      # kitawaga decomposition of rate difference
      kitawaga <-
        DecomposeKitagawa(.$m_jk, .$lead_m_jk, .$pi_jk, .$lead_pi_jk)
      tibble(
        # average mortality at age j
        avg_m_jk =
          sum(.$m_jk*.$pi_jk),
        # rate component of mortality difference
        d_avg_m_jk_rate = kitawaga[1],
        # compositional component of mortality difference
        d_avg_m_jk_comp = kitawaga[2],
        # total mortality difference
        d_avg_m_jk = d_avg_m_jk_rate + d_avg_m_jk_comp
      )
    }) %>%
    ungroup() %>%
    mutate(
      # relative mortality difference
      rd_avg_m_jk = d_avg_m_jk / avg_m_jk,
      # rate component of relative mortality difference
      rd_avg_m_jk_rate = d_avg_m_jk_rate / avg_m_jk,
      # compositional component of relative mortality difference
      rd_avg_m_jk_comp = d_avg_m_jk_comp / avg_m_jk,
      # share of relative mortality decline due to selection
      RDS = rd_avg_m_jk_comp / rd_avg_m_jk
    ) %>%
    rename(
      x_j = x_jk
    )
  decomp$kitagawa_deathrate

  decomp$kitagawa_deathrate_tab <-
    with(decomp$kitagawa_deathrate, {
      tibble(
        # the construction `c(rbind(vector, NA))` injects
        # NAs in between the values of the vector
        age = c(rbind(x_j, NA)),
        rate = c(rbind(avg_m_jk, NA)),
        arrow = rep(c('', '$\\mapsto$'), length.out = 12),
        total_change = c(rbind(NA, d_avg_m_jk)),
        total_change_pct = c(rbind(NA, rd_avg_m_jk)),
        direct_change = c(rbind(NA, d_avg_m_jk_rate)),
        direct_change_pct = c(rbind(NA, rd_avg_m_jk_rate)),
        comp_change = c(rbind(NA, d_avg_m_jk_comp)),
        comp_change_pct = c(rbind(NA, rd_avg_m_jk_comp)),
        selection_effect = c(rbind(NA, RDS)),
      )
    }) %>%
    slice(-nrow(.)) %>%
    mutate_at(
      vars(rate, total_change, direct_change, comp_change),
      ~ formatC(., format = 'e', digits = 1) %>%
        ifelse(. == ' NA', '', .)
    ) %>%
    mutate_at(
      vars(
        total_change_pct,
        direct_change_pct,
        comp_change_pct,
        selection_effect
      ),
      ~ formatC(.*100, format = 'f', digits = 1) %>%
        ifelse(. == 'NA', '', .)
    ) %>%
    mutate_all(
      ~ ifelse(is.na(.), '', .) %>%
        ifelse(. == 'NA', '', .)
    )

}

{

  # decomposition of age-differences in
  # population variance
  decomp$kitagawa_variance <-
    lifetab$neonatal_jk %>%
    filter(x_jk %in% c(0, 1, 7, 14, 21, 27)) %>%
    # for each age add
    # - mean mortality
    # - squared distance of stratum specific mortality rate
    #   to mean mortality
    # - variance of mortality rates
    group_by(j) %>%
    mutate(
      avg_m_jk = sum(m_jk*pi_jk),
      dst_m_jk = (m_jk - avg_m_jk)^2,
      var_m_jk = sum(dst_m_jk*pi_jk)
    ) %>%
    # age differences of squared distances and relative exposures
    arrange(k, j) %>%
    group_by(k) %>%
    mutate(
      lead_dst_m_jk = lead(dst_m_jk),
      lead_pi_jk = lead(pi_jk),
      lead_var_m_jk = lead(var_m_jk),
      diff_dst_m_jk = lead_dst_m_jk - dst_m_jk,
      diff_pi_jk = lead_pi_jk - pi_jk,
      diff_var_m_jk = lead_var_m_jk - var_m_jk
    ) %>%
    group_by(j, x_jk) %>%
    do({
      # kitawaga decomposition of rate difference
      kitawaga <-
        DecomposeKitagawa(
          .$dst_m_jk, .$lead_dst_m_jk, .$pi_jk, .$lead_pi_jk
        )
      tibble(
        # variance of mortality rates
        var_m_jk = .$var_m_jk[1],
        # convergence component of variance difference
        d_var_m_jk_direct = kitawaga[1],
        # compositional component of variance difference
        d_var_m_jk_comp = kitawaga[2],
        # total variance difference
        d_var_m_jk = d_var_m_jk_direct + d_var_m_jk_comp,
        # observed total variance difference
        obs_d_var_m_jk = .$diff_var_m_jk[1]
      )
    }) %>%
    ungroup() %>%
    mutate(
      # relative variance difference
      rd_var_m_jk = d_var_m_jk / var_m_jk,
      # direct (convergence) component of
      # relative variance difference
      rd_var_m_jk_direct = d_var_m_jk_direct / var_m_jk,
      # compositional component of
      # relative variance difference
      rd_var_m_jk_comp = d_var_m_jk_comp / var_m_jk,
      # share of relative variance decline due to selection
      RDS = rd_var_m_jk_comp / rd_var_m_jk
    ) %>%
    rename(
      x_j = x_jk
    )
  decomp$kitagawa_variance

  decomp$kitagawa_variance_tab <-
    with(decomp$kitagawa_variance, {
      tibble(
        # the construction `c(rbind(vector, NA))` injects
        # NAs in between the values of the vector
        age = c(rbind(x_j, NA)),
        var = c(rbind(var_m_jk, NA)),
        arrow = rep(c('', '$\\mapsto$'), length.out = 12),
        total_change = c(rbind(NA, d_var_m_jk)),
        total_change_pct = c(rbind(NA, rd_var_m_jk)),
        direct_change = c(rbind(NA, d_var_m_jk_direct)),
        direct_change_pct = c(rbind(NA, rd_var_m_jk_direct)),
        comp_change = c(rbind(NA, d_var_m_jk_comp)),
        comp_change_pct = c(rbind(NA, rd_var_m_jk_comp)),
        selection_effect = c(rbind(NA, RDS)),
      )
    }) %>%
    slice(-nrow(.)) %>%
    mutate_at(
      vars(var, total_change, direct_change, comp_change),
      ~ formatC(., format = 'e', digits = 1) %>%
        ifelse(. == ' NA', '', .)
    ) %>%
    mutate_at(
      vars(
        total_change_pct,
        direct_change_pct,
        comp_change_pct,
        selection_effect
      ),
      ~ formatC(.*100, format = 'f', digits = 1) %>%
        ifelse(. == 'NA', '', .)
    ) %>%
    mutate_all(
      ~ ifelse(is.na(.), '', .) %>%
        ifelse(. == 'NA', '', .)
    )

}

{

  # decomposition of age-differences in
  # ratio among population mortality to modal mortality
  decomp$kitagawa_memora <-
    lifetab$neonatal_jk %>%
    left_join(lifetab$neonatal_j, by = 'j') %>%
    filter(x_jk %in% c(0, 1, 7, 14, 21, 27)) %>%
    # for each age add
    # - ratio of group specific mortality and apgar 9 mortality
    group_by(j) %>%
    mutate(
      r_jk = m_jk/m_jk[k == lifetab$neonatal_modal_k$k]
    ) %>%
    group_by(k) %>%
    mutate(
      lead_r_jk = lead(r_jk),
      lead_pi_jk = lead(pi_jk),
      diff_r_jk = lead_r_jk - r_jk,
      diff_pi_jk = lead_pi_jk - pi_jk
    ) %>%
    group_by(j, x_jk) %>%
    do({
      # kitawaga decomposition of rate difference
      kitawaga <-
        DecomposeKitagawa(
          .$r_jk, .$lead_r_jk, .$pi_jk, .$lead_pi_jk
        )
      tibble(
        # ratio of population hazard to modal hazard
        r_j = .$m_j[1]/.$m_jk[.$k==lifetab$neonatal_modal_k$k],
        # convergence component of memora difference
        d_r_j_direct = kitawaga[1],
        # compositional component of memora difference
        d_r_j_comp = kitawaga[2],
        # total memora difference
        d_r_j = d_r_j_direct + d_r_j_comp
      )
    }) %>%
    ungroup() %>%
    mutate(
      # relative memora difference
      rd_r_j = d_r_j / r_j,
      # direct (convergence) component of relative memora difference
      rd_r_j_direct = d_r_j_direct / r_j,
      # compositional component of relative memora difference
      rd_r_j_comp = d_r_j_comp / r_j,
      # share of relative memora decline due to selection
      RDS = rd_r_j_comp / rd_r_j
    ) %>%
    rename(
      x_j = x_jk
    )
  decomp$kitagawa_memora

  # format as a nice table
  decomp$kitagawa_memora_tab <-
    with(decomp$kitagawa_memora, {
      tibble(
        # the construction `c(rbind(vector, NA))` injects
        # NAs in between the values of the vector
        age = c(rbind(x_j, NA)),
        r = c(rbind(r_j, NA)),
        arrow = rep(c('', '$\\mapsto$'), length.out = 12),
        total_change = c(rbind(NA, d_r_j)),
        total_change_pct = c(rbind(NA, rd_r_j)),
        direct_change = c(rbind(NA, d_r_j_direct)),
        direct_change_pct = c(rbind(NA, rd_r_j_direct)),
        comp_change = c(rbind(NA, d_r_j_comp)),
        comp_change_pct = c(rbind(NA, rd_r_j_comp)),
        selection_effect = c(rbind(NA, RDS)),
      )
    }) %>%
    slice(-nrow(.)) %>%
    mutate_at(
      vars(r, total_change, direct_change, comp_change),
      ~ formatC(., format = 'd') %>%
        ifelse(. == ' NA', '', .)
    ) %>%
    mutate_at(
      vars(
        total_change_pct,
        direct_change_pct,
        comp_change_pct,
        selection_effect
      ),
      ~ formatC(.*100, format = 'f', digits = 1) %>%
        ifelse(. == 'NA', '', .)
    ) %>%
    mutate_all(
      ~ ifelse(is.na(.), '', .) %>%
        ifelse(. == 'NA', '', .)
    )

}

# Export ----------------------------------------------------------

saveRDS(decomp, file = path$decomp)
