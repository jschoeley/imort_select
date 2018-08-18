#'---
#'title: Life-table and decomposition analysis
#'author: Jonas Schöley
#'date: 2018-02-21
#'---

#'## Strata
#'
#' clinical variables:
#'   - birthweight
#'   - gestation at delivery
#'   - 5 minute apgar score
#'   - presence and severity of congenital anomalies
#'   - plurality
#'   - presence of birth injury
#'   - sex
#'
#' social strata:
#'   - education of mother
#'   - race and hispanic orig of mother
#'   - martial status of mother
#'   - residence status of mother
#'
#' maternal risk factors:
#'   - age of mother
#'   - alcohol or tobacco use during pregnancy

# Init --------------------------------------------------------------------

# install.packages(c('numDeriv', 'tidyverse', 'rlang', 'Hmisc', 'knitr', 'kableExtra', 'haven'))

library(tidyverse)
library(knitr)
library(kableExtra)

# set available memory to 100GB (important only on windows)
#memory.limit(size = 100000)

# Input -------------------------------------------------------------------

load('./priv/data/02-harmonized/2017-10-29-ideath.RData')

# we focus on those born in years 2005 and 2010 and throw away anything
# non essential. memory is precious. we add event and time variables for the
# construction of the survival curve.
ideath_sub <-
  ideath %>%
  # subset to study period
  filter(date_of_delivery_y %in% 2005:2010) %>%
  # add variables for survival analysis
  mutate(
    # age at death in (completed) hours
    # we integrate additional information
    # available for the day of birth
    age_at_death_h =
      ifelse(age_at_death_d > 0,  # if death not at first day
             age_at_death_d*24, # convert age in days to hours
             # otherwise check if death happened in first hour
             # or hour 1-23 and code accordingly
             ifelse(age_at_death_c == '1-23 hours', 1, 0)),
    # so now we have interval censored data on the age at death in hours,
    # to deal with it we add the width of the age interval of death [x, x+nx)
    age_at_death_h_width = 24,
    age_at_death_h_width = ifelse(age_at_death_h == 0, 1, age_at_death_h_width),
    age_at_death_h_width = ifelse(age_at_death_h == 1, 23, age_at_death_h_width),
    # neonatal death indicator, i.e. death during days [0, 28) of life
    death = age_at_death_d < 28, death = ifelse(is.na(death), FALSE, death),
    # by definition a neonatal death is only registered if it occours within
    # the first 28 days of life, therefore we right censor at 28 days or 672 hours
    survtime_d = ifelse(death == TRUE, age_at_death_d, 28),
    # age at death or censoring in (completed) hours
    survtime_h = ifelse(death == TRUE, age_at_death_h, 672),
    survtime_h_width = ifelse(death == TRUE, 24, 0),
    survtime_h_width = ifelse(survtime_h == 0, 1, survtime_h_width),
    survtime_h_width = ifelse(survtime_h == 1, 23, survtime_h_width)
  ) %>%
  select(
    # basic survival information
    death, survtime_h, survtime_h_width,
    # clinical variables
    birthweight_c3, gestation_at_delivery_c4, apgar5_c3,
    congenital_anomalies_c3, plurality_c2, birth_injury, sex,
    # social strata
    education_of_mother_c2, race_and_hispanic_orig_of_mother_c2,
    martial_status_of_mother,
    # maternal risk factors
    age_of_mother_c3,
    alcto_use_during_pregnancy
  )
rm(ideath)

# export a stata copy
# haven::write_dta(ideath_sub,
#                  path = './priv/data/02-harmonized/2017-10-29-ideath.dta',
#                  version = 14)

# Analysis functions ------------------------------------------------------

#' Calculate Life-tables From Individual Level Survival Times
#'
#' Individual level survival times may be interval or right censored.
#' The life-table can be arbitrarily abridged.
#'
#' @param df a data frame
#' @param x time until death or censoring
#' @param nx width of interval [x, x+nx), i.e. precision of measurement
#' @param death death (TRUE) or censored (FALSE)
#' @param cuts cutpoint for age intervals in abridged life-table [x1, ..., xn)
#' @param ... grouping variable
GetLifeTable <- function (df, x, nx, death, cuts, ...) {
  x_i = enquo(x); nx_i = enquo(nx); death_i = enquo(death); strata = quos(...)

  FindIntervalStart <- function (x, breaks) {
    breaks[.bincode(x = x, breaks = breaks, right = FALSE, include.lowest = FALSE)]
  }

  lt <-
    df %>%
    # aggregation into pre-defined age-groups
    mutate(x0 = FindIntervalStart(!!x_i, cuts)) %>%
    group_by(..., x0, add = FALSE) %>%
    summarise(
      nDx = sum(!!death_i),
      nCx = sum(!(!!death_i)),
      # average time spent in interval for
      # those who leave during interval (by death or censoring)
      nax = mean(!!x_i) - first(x0) + 0.5*mean(!!nx_i)
    ) %>%
    arrange(..., x0) %>%
    group_by(..., add = FALSE) %>%
    # calculation of life-table columns
    mutate(
      nx = c(diff(x0), last(cuts)-last(x0)),
      # assuming no late entry
      Nx = head(cumsum(c(sum(nDx, nCx), -(nCx+nDx))), -1),
      # distribution of censoring or death
      nfx = (nCx + nDx) / sum(nCx + nDx),
      # life-table
      nqx = nDx/Nx,
      lx = Nx/first(Nx),
      ndx = lx*nqx,
      nLx = lead(lx, n = 1, default = NA) * nx + (nfx*nax),
      nLx = ifelse(is.na(nLx), last(nax), nLx),
      nmx = ndx / nLx
    ) %>% ungroup() %>%
    mutate(id = group_indices(., ...)) %>%
    # fill 'gaps' in life-table
    complete(x0 = head(cuts, -1), distinct(., id, !!!strata),
             fill = list(nDx = 0, nCx = 0, ndx = 0, nmx = 0, nqx = 0, nax = 0)) %>%
    arrange(id, x0) %>%
    group_by(id, add = FALSE) %>%
    mutate(
      nx = c(diff(x0), last(cuts)-last(x0)),
      Nx = head(cumsum(c(sum(nDx, nCx), -(nCx+nDx))), -1),
      nEx = (Nx-nDx)*nx + nax*(nDx+nCx),
      nEx = ifelse(is.infinite(nEx), nax*(nDx+nCx), nEx),
      lx = Nx/first(Nx),
      nLx = ifelse(is.na(nLx), lx*nx, nLx),
      Tx = rev(cumsum(rev(nLx))), ex = Tx/lx
    ) %>% ungroup() %>%
    select(id, ..., x = x0, nx, Nx, nEx, nDx, nCx,
           lx, ndx, nqx, nax, nmx, nLx, Tx, ex)

  return(lt)

}

#' Fit Cubic-Splines to Life-table Columns
#'
#' @param df a data frame
#' @param x name of age column (unquoted)
#' @param lx name of survival column (unquoted)
#' @details We fit a spline to the log-transformed survival. In order to avoid
#' taking the log of 0 and to avoid fitting a Hyman filtered spline to log(y) =
#' 0 (the spline performs bad if y = 0), we use the value transformation $g(y) =
#' \log(y+2)$, with the inverse $g^{-1}(y) = exp(y)-2$. We transform the domain
#' as well using $q(x) = \log(x+1)$ with inverse $q^{-1}(x)=1/(x-1)$ in order to
#' capture the rapid decline in hazard of death right right after birth observed
#' in infant mortality.
FitLTSpline <- function (df, x, lx) {
  require(rlang)
  x = enquo(x); lx = enquo(lx)

  logloglxFun <- splinefun(x = log(unlist(df[,quo_name(x)])+1),
                           y = log(unlist(df[,quo_name(lx)])+2),
                           method = 'hyman')

  fnct <- function (x, type = 'lx') {
    lxFun <- function (x) {
      lx = exp(logloglxFun(log(x+1)))-2
      # in rare occurences the lx will be negative (bordering to 0)
      # we fix that
      ifelse(lx<0, 0, lx)
    }
    dxFun <- function (x) { -lxFun(x) * (1/(x+1)) * logloglxFun(log(x+1), deriv = 1) }
    hxFun <- function (x) {
      hx = dxFun(x) / lxFun(x)
      # set hazard to 0 if lx = 0
      ifelse(is.nan(hx), 0, hx)
    }
    ddxhxFun <- function (x) { numDeriv::grad(hxFun, x) }
    switch(type, lx = lxFun(x), dx = dxFun(x), hx = hxFun(x), ddxhx = ddxhxFun(x))
  }

  return(fnct)
}

#' @param df a data frame
#' @param Nx name of age-specific population size column (unquoted)
#' @param nDx name of age-specific death count column (unquoted)
#' @param nCx name of age-specific censoring count column (unquoted)
#' @param ... names of stratificatio variables (unquoted)
SummariseLT <- function(df, Nx, nDx, nCx, nEx, nx, nax, ...) {
  Nx = enquo(Nx); nDx = enquo(nDx); nCx = enquo(nCx); nEx = enquo(nEx); nx = enquo(nx)

  df %>%
    group_by(..., add = FALSE) %>%
    summarise(
      N = first(Nx),
      E = sum(nEx),
      D = sum(nDx),
      C = sum(nCx),
      probD = D / N,
      CMR = D / E
    ) %>%
    ungroup() %>%
    mutate(pN = N / sum(N),
           pD = D / sum(D))
}

#' Vaupel-Zhang Decomposition
#'
#' Decompose the change in mortality over x into a direct and a compositional
#' component.
#'
#' @param df a data frame
#' @param x time/age
#' @param mux_z mortality in group z
#' @param dmux_z derivative of mortality in group z
#' @param pxz joint probability of being in x and z
#' @param pz_x conditional probability of z given x
DecomposeVR <- function (df, x, mux_z, dmux_z, pxz, pz_x) {
  x = enquo(x); mux_z = enquo(mux_z); dmux_z = enquo(dmux_z)
  pxz = enquo(pxz); pz_x = enquo(pz_x)

  df %>%
    group_by(!!x, add = FALSE) %>%
    summarise(
      # P(X>x) = SUM_z [P(X>x & z)]
      lx = sum(!!pxz),
      mux = Hmisc::wtd.mean(!!mux_z, !!pz_x),
      direct_fx = Hmisc:: wtd.mean(!!dmux_z, !!pz_x),
      compos_fx = -Hmisc:: wtd.var(!!mux_z, !!pz_x, method = 'ML'),
      dmux = direct_fx + compos_fx,
      p_compos_fx = compos_fx / dmux
    ) %>% ungroup() %>%
    mutate(
      # d(x) = mu(x) * S(x)
      dx = mux*lx
    ) %>%
    select(x, lx,dx, mux, dmux, direct_fx, compos_fx, p_compos_fx)
}

IntegrateVZDecomposition <- function (df, x, dmux, direct, compos, xout) {
  x = enquo(x); dmux = enquo(dmux); direct = enquo(direct); compos = enquo(compos)

  length_xout <- length(xout)

  # linearly interpolate the decomposition results
  dmuxFun <- approxfun(x = df[[quo_name(x)]],
                       y = df[[quo_name(dmux)]],
                       method = 'linear')
  directfxFun <-approxfun(x = df[[quo_name(x)]],
                          y = df[[quo_name(direct)]],
                          method = 'linear')
  composfxFun <- approxfun(x = df[[quo_name(x)]],
                           y = df[[quo_name(compos)]],
                           method = 'linear')

  diff_mux <- vector('numeric', length = length_xout-1)
  diff_direct <- vector('numeric', length = length_xout-1)
  diff_compos <- vector('numeric', length = length_xout-1)

  # integrate results
  for (i in 1:(length_xout-1)) {
    diff_mux[i] <-
      integrate(dmuxFun,
                lower = xout[i], upper = xout[i+1]
      )$value
    diff_direct[i] <-
      integrate(directfxFun,
                lower = xout[i], upper = xout[i+1]
      )$value
    diff_compos[i] <-
      integrate(composfxFun,
                lower = xout[i], upper = xout[i+1]
      )$value
  }

  data_frame(
    x = head(xout, -1),
    nx = diff(xout),
    diff_mux,
    diff_direct,
    p_direct = diff_direct/diff_mux,
    diff_compos,
    p_compos = diff_compos/diff_mux
  )

}

# Workflow functions ------------------------------------------------------

AnalyzeThis <- function(df, x, nx, death, cuts, xout, xlab, title, ...) {

  x_i = enquo(x); nx_i = enquo(nx); death_i = enquo(death)

  # 1. stratified individual level survival times to abridged life-table
  lt <- GetLifeTable(df, x = !!x_i, nx = !!nx_i, death = !!death_i,
                     cuts = cuts, ...)

  # 2. summary statistics by strata
  lt_summary <- SummariseLT(lt, Nx = Nx, nDx = nDx, nCx = nCx, nEx = nEx,
                            nx = nx, nax = nax,
                            id, ...)

  # 3. interpolate life-table using splines
  lt_interpol <-
    left_join(lt, lt_summary) %>%
    group_by(id, ..., add = FALSE) %>%
    # interpolate life-tables using splines
    do({
      SP <- FitLTSpline(., x = x, lx = lx)
      data_frame(
        x = xout,
        lx_z = SP(xout, type = 'lx'),
        mux_z = SP(xout, type = 'hx'),
        dmux_z = SP(xout, type = 'ddxhx'),
        # probability to be in group z and to survive to x P(X>x & z)
        pxz = .$pN[1]*lx_z
      )
    }) %>%
    group_by(x, add = FALSE) %>%
    # probability to be in group z given survival to age x: P(z | X>x)
    # P(X>x & z) = P(X>x) * P(z | X>x)
    # P(z | X>x) = P(X>x & z) / P(X>x) = P(X>x & z) / sum_z[ P(X>x & z) ]
    mutate(pz_x = pxz / sum(pxz)) %>% ungroup()

  # 4. calculate percent point change in population share by z over study period
  diff_pz <-
    lt_interpol %>%
    filter(x %in% range(x)) %>%
    group_by(...) %>%
    summarise(
      pz1 = first(pz_x),
      pz2 = last(pz_x),
      diff_pz = diff(pz_x)
    )

  # 5. decompose
  lt_decomp <- DecomposeVR(lt_interpol,
                           x = x, mux_z = mux_z, dmux_z = dmux_z,
                           pxz = pxz, pz_x = pz_x)
  # integrate decomposition results over cuts
  lt_decompo_disc <-
    IntegrateVZDecomposition(lt_decomp, x = x, dmux = dmux, direct = direct_fx,
                             compos = compos_fx, xout = cuts[!is.infinite(cuts)])

  # 6. plot
  y_break = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 5e-01)
  x_break = unique(lt[['x']])

  plot_data <-
    lt_summary %>%
    mutate(cmr_rank = rank(CMR)) %>%
    right_join(., lt_interpol)

  hzrd_plot <-
    ggplot(plot_data) +
    annotation_logticks(sides = 'lr', scaled = TRUE, base = 10, colour = 'grey92') +
    geom_line(aes(x = x, y = mux_z, color = cmr_rank, group = id)) +
    geom_text(aes(x = x, y = mux_z, label = id, color = cmr_rank, group = id),
              data = filter(plot_data, x == 1)) +
    geom_line(aes(x = x, y = mux),
              color = 'red', size = 1,
              data = lt_decomp) +
    scale_y_continuous('Deaths per person-hour of exposure',
                       trans = 'log10', breaks = y_break,
                       limits = c(1e-07, 0.5),
                       expand = c(0,0)) +
    scale_x_continuous('Age',
                       trans = 'sqrt', breaks = x_break, labels = xlab,
                       expand = c(0,0)) +
    theme_minimal() +
    theme(aspect.ratio = 0.7,
          panel.grid.minor = element_blank(),
          legend.position = 'none') +
    ggtitle(title)

  list(lt_summary = lt_summary, lt = lt,
       lt_interpol = lt_interpol, diff_pz = diff_pz,
       lt_decomp = lt_decomp, lt_decompo_disc = lt_decompo_disc,
       hzrd_plot = hzrd_plot)

}

# Decomposition analysis --------------------------------------------------

# cut points for abridged infant life-table
intervals_h = c(0, 1, 24, 48, 72, 96, 120, 144, 168, 336, 504, 672, Inf)
intervals_h_fine = seq(0, 672, 0.1)

xlab = c('Birth', '1 hr', '1 day', '2d', '3d', '4d', '5d', '6d',
         '7d', '14d', '21d', '28d')

## margins

# overall
ilt_overall <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'Overall')

## all strata

# all strata strata
ilt_multi_strata <-
  AnalyzeThis(
    ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
    cuts = c(0, 24, 168, Inf), xout = c(0, 24, 168, 672),
    xlab = c('0', '24', '168'), title = 'Interaction',
    birthweight_c3, gestation_at_delivery_c4, apgar5_c3,
    congenital_anomalies_c3, plurality_c2
    #sex, birth_injury,
    #education_of_mother_c2, race_and_hispanic_orig_of_mother_c2,
    #martial_status_of_mother,
    #age_of_mother_c3, alcto_use_during_pregnancy
  )

ilt_multi_strata$lt %>%
  mutate_at(vars(birthweight_c, gestation_at_delivery_c, apgar5,
            congenital_anomalies, plurality),
            funs(as.factor(ifelse(is.na(.), 'unknown', as.character(.))))) %>%
  filter(x %in% c(0, 1)) %>%
  haven::write_dta(data = .,
                   path = './out/2017-09-01-ilt_multi_strata.dta',
                   version = 13)

## clinical strata

# birthweight
ilt_birthweight_c3 <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by birthweight',
              birthweight_c3); ilt_birthweight_c3
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_birthweight_c3.pdf'), width = 7, height = 5)

# gestation at delivery
ilt_gestation_at_delivery_c4 <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by gestation at delivery',
              gestation_at_delivery_c4); ilt_gestation_at_delivery_c4
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_gestation_at_delivery_c4.pdf'), width = 7, height = 5)

# apgar5
ilt_apgar5_c3 <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by APGAR 5 score',
              apgar5_c3); ilt_apgar5_c3
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_apgar5_c3.pdf'), width = 7, height = 5)

# congenital anomalies
ilt_congenital_anomalies_c3 <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by presence and severity of congenital conditions',
              congenital_anomalies_c3); ilt_congenital_anomalies_c3
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_congenital_anomalies_c3.pdf'), width = 7, height = 5)

# plurality
ilt_plurality <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by plurality',
              plurality); ilt_plurality
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_plurality.pdf'), width = 7, height = 5)

# birth injury
ilt_birth_injury <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by occurence of birth injuries',
              birth_injury); ilt_birth_injury
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_birth_injury.pdf'), width = 7, height = 5)

# sex
ilt_sex <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by sex',
              sex); ilt_sex
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_sex.pdf'), width = 7, height = 5)

## maternal risk factors

# mothers age
ilt_age_of_mother_c3 <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by age of mother',
              age_of_mother_c3); ilt_age_of_mother_c3
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_age_of_mother_c3.pdf'), width = 7, height = 5)

# tobacoo or alcohol use during pregnancy
ilt_tobalc_use <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by occurence of tobacco or alcohol consumption during pregnancy',
              alcto_use_during_pregnancy); ilt_tobalc_use
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_tobalc_use.pdf'), width = 7, height = 5)

## social strata

# mothers origin
ilt_mothers_origin <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by origin of mother',
              race_and_hispanic_orig_of_mother_c2); ilt_mothers_origin
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_mothers_origin.pdf'), width = 7, height = 5)

# mothers education
ilt_mothers_education <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by education of mother',
              education_of_mother_c2); ilt_mothers_education
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_mothers_education.pdf'), width = 7, height = 5)


# mothers martial status
ilt_mothers_martial_status <-
  AnalyzeThis(ideath_sub, x = survtime_h, nx = survtime_h_width, death = death,
              cuts = intervals_h, xout = intervals_h_fine,
              xlab = xlab, title = 'by martial status of mother',
              martial_status_of_mother); ilt_mothers_martial_status
ggsave(paste0('./out/fig/', Sys.Date(), '-ilt_martial_status_of_mother.pdf'), width = 7, height = 5)

# Frailty analysis --------------------------------------------------------

frail <-
  ilt_multi_strata$lt_summary %>%
  mutate(z = CMR/ilt_overall$lt_summary$CMR)

x = rep(frail$z, times = frail$N)
x = x[x > 0 & !is.na(x)]
fitdistrplus::fitdist(x, distr = "gamma")


ilt_multi_strata$lt_summary %>%
  ggplot() +
  geom_histogram(aes(x = CMR/ilt_overall$lt_summary$CMR,
                     y = ..density.., weight = pN),
                 bins = 50) +
  geom_vline(xintercept = 1) +
  scale_x_continuous(trans = 'log', breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)) +
  theme_minimal()

# Print summary table -----------------------------------------------------

bind_rows(
  lapply(list(
    birthweight_c = ilt_birthweight_c,
    gestation_at_delivery = ilt_gestation_at_delivery_c,
    apgar5 = ilt_apgar5,
    congenital_anomalies = ilt_congenital_anomalies,
    plurality = ilt_plurality,
    birth_injury = ilt_birth_injury,
    sex = ilt_sex,
    mothers_education = ilt_mothers_education,
    mothers_origin = ilt_mothers_origin,
    mothers_martial_status = ilt_mothers_martial_status,
    mothers_residence_status = ilt_mothers_residence_status,
    age_of_mother = ilt_age_of_mother,
    tobalc_use = ilt_tobalc_use
  ), '[[', 'lt_summary'),
  .id = 'stratum'
) %>% mutate_at(vars(-stratum, -id, -N, -E, -D, -C, -probD, -CMR, -pN, -pD),
                funs(ifelse(is.na(.), '', as.character(.)))) %>%
  unite(col = cat, -stratum, -id, -N, -E, -D, -C, -probD, -CMR, -pN, -pD, sep = '') %>%
  mutate(cat = ifelse(cat == '', NA, cat)) %>%
  mutate(
    pN = pN*100,
    pD = pD*100,
    probD_ratio = probD  / ilt_overall$lt_summary$probD*100
    ) %>%
  select(cat, N, pN, D, pD, probD, probD_ratio) %>%
  mutate_at(vars(N, D),
            formatC, format = 'f', digits = 0, big.mark = ',') %>%
  mutate_at(vars(pN, pD),
            formatC, format = 'f', digits = 1) %>%
  mutate_at(vars(probD),
            formatC, format = 'f', digits = 4) %>%
  mutate_at(vars(probD_ratio),
            formatC, format = 'f', digits = 2) %>%
  kable(align = c('r'), format = 'latex', booktabs = TRUE, longtable = TRUE,
        escape = FALSE, row.names = FALSE,
        col.names = c('stratum',
                      'N', 'pct.',
                      'D', 'pct.',
                      'by sub-group', 'compared to margin')
  ) %>%
  kable_styling(latex_options = 'repeat_header', font_size = 6) %>%
  add_header_above(c(' ' = 1, 'Births' = 2, 'Neonatal deaths' = 2,
                     'Probability of neonatal death' = 2)) %>%
  group_rows(group_label = 'Birthweight',
             start_row =  1, end_row = 6) %>%
  group_rows(group_label = 'Gestation at birth',
             start_row =  7, end_row = 14) %>%
  group_rows(group_label = '5 minute APGAR score',
             start_row =  15, end_row = 26) %>%
  group_rows(group_label = 'Presence and severity of congenital anomalies',
             start_row =  27, end_row = 30) %>%
  group_rows(group_label = 'Plurality',
             start_row =  31, end_row = 34) %>%
  group_rows(group_label = 'Presence of birth injury',
             start_row =  35, end_row = 37) %>%
  group_rows(group_label = 'Sex',
             start_row =  38, end_row = 39) %>%
  group_rows(group_label = 'Education of mother',
             start_row =  40, end_row = 44) %>%
  group_rows(group_label = 'Race and hispanic origin of mother',
             start_row =  45, end_row = 49) %>%
  group_rows(group_label = 'Martial status of mother',
             start_row =  50, end_row = 51) %>%
  group_rows(group_label = 'Residence status of mother',
             start_row =  52, end_row = 55) %>%
  group_rows(group_label = 'Age of mother',
             start_row =  56, end_row = 61) %>%
  group_rows(group_label = 'Alcohol or tobacco use during pregnancy',
             start_row =  62, end_row = 64)

# Print decomposition table -----------------------------------------------

table_ages = c(0, 1, 24, 168, 672)

bind_rows(
  birthweight_c = ilt_birthweight_c$lt_decomp,
  gestation_at_delivery = ilt_gestation_at_delivery_c$lt_decomp,
  apgar5 = ilt_apgar5$lt_decomp,
  congenital_anomalies = ilt_congenital_anomalies$lt_decomp,
  plurality = ilt_plurality$lt_decomp,
  birth_injury = ilt_birth_injury$lt_decomp,
  sex = ilt_sex$lt_decomp,
  mothers_education = ilt_mothers_education$lt_decomp,
  mothers_origin = ilt_mothers_origin$lt_decomp,
  mothers_martial_status = ilt_mothers_martial_status$lt_decomp,
  mothers_residence_status = ilt_mothers_residence_status$lt_decomp,
  age_of_mother = ilt_age_of_mother$lt_decomp,
  tobalc_use = ilt_tobalc_use$lt_decomp,
  #multi_strata = ilt_multi_strata$lt_decomp,
  .id = 'stratum'
) %>%
  filter(x %in% table_ages) %>%
  group_by(stratum) %>%
  mutate(p_direct_fx = direct_fx/dmux,
         rel_dmux = dmux / mux) %>% ungroup() %>%
  select(x, mux, dmux, rel_dmux, direct_fx, p_direct_fx, compos_fx, p_compos_fx) %>%
  mutate_at(vars(rel_dmux, p_direct_fx, p_compos_fx),
            formatC, format = 'f', digits = 2) %>%
  mutate_at(vars(mux, dmux, direct_fx, compos_fx),
            formatC, format = 'e', digits = 1) %>%
  kable(align = c('c'), format = 'latex', booktabs = TRUE, longtable = TRUE,
        escape = FALSE, row.names = FALSE,
        col.names = c('$x$', '$\\bar\\mu(x)$',
                      '$\\dot{\\bar{\\mu}}(x)$',
                      '$\\dot{\\bar{\\mu}}(x)/\\bar\\mu(x)$',
                      '$\\bar{\\dot{\\mu}}_z(x)$',
                      '$\\bar{\\dot{\\mu}}_z(x)/\\dot{\\bar{\\mu}}(x)$',
                      '$-\\sigma^2_{\\mu_z}(x)$',
                      '$-\\sigma^2_{\\mu_z}(x)/\\dot{\\bar{\\mu}}(x)$')
  ) %>%
  kable_styling(latex_options = 'repeat_header', font_size = 6) %>%
  add_header_above(c('Age' = 1, 'Population hazard' = 1,
                     'Absolute and relative derivative of population hazard' = 2,
                     'Absolute and relative direct effect' = 2,
                     'Absolute and relative selection effect' = 2)) %>%
  group_rows(group_label = 'Birthweight',
             start_row =  1, end_row = 5) %>%
  group_rows(group_label = 'Gestation at birth',
             start_row =  6, end_row = 10) %>%
  group_rows(group_label = '5 minute APGAR score',
             start_row =  11, end_row = 15) %>%
  group_rows(group_label = 'Presence and severity of congenital anomalies',
             start_row =  16, end_row = 20) %>%
  group_rows(group_label = 'Plurality',
             start_row =  21, end_row = 25) %>%
  group_rows(group_label = 'Presence of birth injury',
             start_row =  26, end_row = 30) %>%
  group_rows(group_label = 'Sex',
             start_row =  31, end_row = 35) %>%
  group_rows(group_label = 'Education of mother',
             start_row =  36, end_row = 40) %>%
  group_rows(group_label = 'Race and hispanic origin of mother',
             start_row =  41, end_row = 45) %>%
  group_rows(group_label = 'Martial status of mother',
             start_row =  46, end_row = 50) %>%
  group_rows(group_label = 'Residence status of mother',
             start_row =  51, end_row = 55) %>%
  group_rows(group_label = 'Age of mother',
             start_row =  56, end_row = 60) %>%
  group_rows(group_label = 'Alcohol or tobacco use during pregnancy',
             start_row =  61, end_row = 65)