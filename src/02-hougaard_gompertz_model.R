# Analyze Hougaard-Gompertz model fit to neonatal life-tables
#
# Jonas Schöley

# Init ------------------------------------------------------------

library(tidyverse)
library(brms)

path <- list(
  lifetab = 'out/lifetab.rds',
  hougaard = 'out/hougaard.rds',
  out = 'out'
)

hougaard <- list()

source('src/00-global.R')

# Load data -------------------------------------------------------

lifetab <- readRDS(path$lifetab)

# Functions -------------------------------------------------------

# We fit the a multiplicative frailty model with a frailty
# distribution invented by Hougaard (1986) and a Gompertz
# baseline hazard to the infant life-table for the total population
# and make predictions regarding the role of mortality selection.

# h0(x) = a*exp(b*x)
# H0(x) = a*(exp(b*x)-1)/b
# h(x|z) = z*h0(x)
# Z ~ Hougaard(x; alpha, delta, theta)
#
# restricting E[Z] = 1 gives parametrization
#
# Z ~ Hougaard(x; alpha, theta^(1-alpha), theta),
#
# resulting in population hazard
#
# h(x) = a*exp(b*x)*((theta + a*(exp(b*x)-1)/b)/theta)^(alpha-1)
#
# with logarithm
#
# log h(x) = log(a) + b*x +
#   (alpha-1)*log((theta + exp(log(a))*(exp(b*x)-1)/b)/theta)
#
# which we reparametrize as
#
# log h(x) = b0 + b1*x +
#   (logit(p)-1)*log((exp(c) + exp(b0)*(exp(b1*x)-1)/b1)/exp(c))
#
# where
#
# log(a) = b0 <-> a = exp(b0)
# b = b1
# alpha-1 = logit(p)-1 <-> alpha = logit(p)
# theta = exp(c)
#
# To recover the distribution of hazards h(x|z)
# from the distribution of frailties Z we use the fact that
# if Z ~ Hougaard(alpha, delta, theta) and h0(x)>0,
# the distribution of h0(x)Z is P(alpha, h0(x)^alpha*delta, theta/h0(x))
# see Hougaard (1986), Lemma 1e.

# Gompertz baseline hazard
GompertzHzrd <- function (x, a, b) {
  a*exp(b*x)
}
# Gompertz cumulative baseline hazard
GompertzCumHzrd <- function (x, a, b) {
  a*(exp(b*x)-1)/b
}
# Hougaard-Gompertz population hazard
HougaardGompertzModel <- function (x, a, b, alpha, theta) {
  a*exp(b*x)*((theta + a*(exp(b*x)-1)/b)/theta)^(alpha-1)
}
# Hougaard density of frailties
HougaardDensity <- function (x, alpha, delta, theta, max_k = 10) {
  term1 <-
    -exp(-theta*x + delta*theta^alpha/alpha)*(pi*x)^(-1)
  M <- matrix(NA, nrow = length(x), ncol = max_k)
  for (k in 1:max_k) {
    M[,k] <-
      gamma(k*alpha+1)/factorial(k)*
      (-delta*x^(-alpha)/alpha)^k*
      sin(alpha*k*pi)
  }
  term2 <- rowSums(M)
  return(term1*term2)
}

# Fit Hougaard-Gompertz model -------------------------------------

# fit Hougaard-Gompertz model to population
# infant life-table
hougaard$fit <-
  brm(
    bf(
      D_j ~ b0 + b1*x_j +
        (1/(1+exp(-p))-1)*
        log((exp(c) + exp(b0)*(exp(b1*x_j)-1)/b1)/exp(c)) +
        log(O_j),
      b0 + b1 + p + c ~ 1,
      nl = TRUE
    ),
    family =
      # negative binomial had large phi parameter (~200)
      # and in the brms parametrization that means that
      # there's no evidence for overdispersion
      poisson(link = 'log'),
    prior =
      c(
        prior(normal(-8, 1), nlpar = 'b0'),
        prior(normal(-0.03, 1), nlpar = 'b1'),
        prior(normal(0.5, 1), nlpar = 'p'),
        prior(normal(-6, 1), nlpar = 'c')
      ),
    data =
      lifetab$infant_j,
    chains = 8, cores = 8, iter = 1500, warmup = 500,
    seed = 42,
    control = list(max_treedepth = 15, adapt_delta = 0.95)
  )

summary(hougaard$fit)

# recover the parameters of the Hougaard-Gompertz model
hougaard$paras <-
  hougaard$fit %>%
  {c(
    a = exp(fixef(.)[1,'Estimate']),
    b = fixef(.)[2,'Estimate'],
    alpha = plogis(fixef(.)[3,'Estimate']),
    theta = exp(fixef(.)[4,'Estimate'])
  )}

# Explore frailties -----------------------------------------------

# normalized density of frailties at age 0
# for Hougaard-Gompertz fit
HougaardDensity(
  10^seq(-11, 5, 0.01),
  alpha = hougaard$paras[3],
  delta = hougaard$paras[4]^(1-hougaard$paras[3]),
  theta = hougaard$paras[4],
  max_k = 300
) %>% {
  plot(x = 10^seq(-11, 5, 0.01),
       y = ./max(.,na.rm = T), log = 'xy',
       type = 'l')
}

z <- 10^seq(-11, 5, 0.001)
dz <- c(diff(z), 0)
d <- HougaardDensity(
  z,
  alpha = hougaard$paras[3],
  delta = hougaard$paras[4]^(1-hougaard$paras[3]),
  theta = hougaard$paras[4],
  max_k = 300
)
sum(d*dz*z, na.rm = T)

tibble(
  x = 10^seq(-11, 1, 0.001),
  dx = c(diff(x), 0),
  fx = HougaardDensity(
    x,
    alpha = hougaard$paras[3],
    delta = hougaard$paras[4]^(1-hougaard$paras[3]),
    theta = hougaard$paras[4],
    max_k = 300
  ) %>% ifelse(is.nan(.), 0, .),
  Fx =
    cumsum(fx*dx),
  running_expectation =
    cumsum(fx*dx*x)
) %>%
  {plot(
    .$x, .$Fx, log = 'x', type = 'l'
  )}

# Plot Hougaard-Gompertz model ------------------------------------

# Hougaard-Gompertz population hazard
hougaard$pop_hzrd <-
  tibble(
    x = seq(0, 365, 0.1),
    h_x =
      HougaardGompertzModel(
        x,
        a = hougaard$paras[1],
        b = hougaard$paras[2],
        alpha = hougaard$paras[3],
        theta = hougaard$paras[4]
      )
  )
# Hougaard-Gompertz conditional hazards
# for selected frailty levels
hougaard$cond_hzrds <-
  expand_grid(
    x = seq(0, 365, 1),
    z = c(10, 1, 0.1, 0.01, 0.001)
  ) %>%
  mutate(
    h_xz =
      z*GompertzHzrd(
        x,
        a = hougaard$paras[1],
        b = hougaard$paras[2]
      )
  )

# plot Hougaard-Gompertz hazard against life-table death rates
hougaard$plot_hazard <-
  lifetab$infant_j %>%
  ggplot() +
  geom_point(
    aes(x = x_j, y = m_j, size = D_j),
    color = 'grey75', show.legend = FALSE
  ) +
  geom_line(
    aes(x = x, y = h_x),
    size = 1,
    data = hougaard$pop_hzrd
  ) +
  geom_line(
    aes(x = x, y = h_xz, group = z),
    size = 0.2,
    data = hougaard$cond_hzrds
  ) +
  geom_text(
    aes(x = -5, y = h_xz, label = z),
    size = 3, hjust = 1,
    data = hougaard$cond_hzrds %>% filter(x == 0)
  ) +
  annotate(
    'text',
    x = 200, y = 10^(-1),
    label = expression(h[0](x)==ae^{-bx}),
    hjust = 0,
  ) +
  annotate(
    'text',
    x = 200, y = 5*10^(-2),
    label = expression(H[0](x)==a(e^{bx}-1)/b),
    hjust = 0,
  ) +
  annotate(
    'text',
    x = 200, y = 2.3*10^(-2),
    label = expression(h[Z](x)==Zh[0](x)),
    hjust = 0,
  ) +
  annotate(
    'text',
    x = 200, y = 1.1*10^(-2),
    label = expression(Z%~%Hougaard(alpha, theta^{1-alpha}, theta)),
    hjust = 0,
  ) +
  annotate(
    'text',
    x = 200, y = 5*10^(-3),
    label = expression(h(x)==h[0](x)*((theta + H[0](x))/theta)^{alpha-1}),
    hjust = 0,
  ) +
  scale_y_continuous(
    trans = 'log10',
    breaks = figspec$bl$hazard2,
    labels = parse(text = names(figspec$bl$hazard2)),
  ) +
  scale_size_area() +
  coord_cartesian(
    ylim = c(6*10^(-7), 10^(-1))
  ) +
  figspec$MyGGplotTheme(size = 7, ar = 0.7, axis = 'xy', grid = '') +
  labs(
    x = 'Days since birth',
    y = 'Deaths per person-day'
  )

hougaard$plot_hazard

# Export ----------------------------------------------------------

figspec$ExportFigure(
  hougaard$plot_hazard,
  filename = 'hougaard',
  path = path$out, device = 'pdf'
)
