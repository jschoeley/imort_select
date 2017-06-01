###############################
# CALCULATE INFANT LIFETABLES #
###############################

# Based on the individual level data on births and infant deaths we calculate
# infant lifetables for various groups

library(dplyr)
library(readr)
library(ggplot2)

load("./priv/data/02-harmonized/ideath.RData")

# The cohort lifetable function for individual level data -----------------

#' A Lifetable from survival data
LT <- function (time, event) {

  # highest transition age
  max_age = max(time, na.rm = TRUE)
  # integer age vector
  x = 0L:max_age
  # initial cohort size
  N0 = length(time)

  # tabulate deaths and censorings at each point in time
  tab <-
    data.frame(
      x,
      lapply(list(time[event==1], time[event==0]),
             function(x) tabulate(x+1L, nbins = max_age+1L))
    )
  names(tab) <- c("x", "Dx", "Cx")

  ##### ah shit, i need to hack that.
  tab$x_start = tab$x%/%24*24
  tab$x_start[tab$x == 0] = 0
  tab$x_start[tab$x %in% 1:23] = 1
  tab$x_width = 24
  tab$x_width[tab$x_start == 0] = 1
  tab$x_width[tab$x_start %in% 1:23] = 23
  tab = aggregate(tab[-1], list(x = tab$x_start, nx = tab$x_width), sum)[-c(5,6)]
  # number of ages
  J = length(tab$x)
  #####

  # number of survivors at start of age interval x
  # (persons alive at start of interval - persons died in interval -
  #  persons censored in interval)
  Nx = c(N0, N0 - cumsum(tab$Dx + tab$Cx))
  # censoring adjusted number of survivors at start of interval
  # person-time of exposure to risk of event during interval [x,x+nx)
  # assuming censoring or death occoured mid-interval
  Ex = (Nx[-1] + 0.5*tab$Cx + 0.5*tab$Dx)*tab$nx
  # probability of dying in x given survival to x
  qx = tab$Dx/Nx[-(J+1)]
  # adjusted probability of surviving to x+nx given survival to x
  px = 1-qx
  # probability of surviving to x (lifetable survivor function)
  lx = c(1, cumprod(px))
  # density of deaths
  dx = lx[-(J+1)]*qx
  # hazard of death in interval x
  hx = tab$Dx/Ex

  lt = data.frame(x = tab$x, nx = tab$nx,
                  Dx = tab$Dx, Cx = tab$Cx,
                  Nx = Nx[-(J+1)], Ex,
                  qx, px,
                  dx, lx = lx[-(J+1)],
                  hx)

  return(lt)

}

# Calculate infant cohort lifetables --------------------------------------

# whole population
ideath %>%
  filter(date_of_delivery_y %in% 2005:2010) %>%
  do(LT(time = .$age_at_death_or_cens_h, event = .$death)) -> ilt1

# by apgar 5 score
ideath %>%
  filter(date_of_delivery_y %in% 2005:2010) %>%
  group_by(apgar5) %>%
  do(LT(time = .$age_at_death_or_cens_h, event = .$death)) -> ilt2

# Decompose ---------------------------------------------------------------

ilt2 %>%
  group_by(apgar5) %>%
  mutate(dhx = c(diff(hx), NA)) %>%
  group_by(x) %>%
  summarise(
    adaptation = mean(dhx),
    selection  = -var(hx),
    rel        = selection / (selection+adaptation)
  ) -> decomp

# Plot infant lifetables --------------------------------------------------

ilt_x_break = c(0, 1, 24, 24*7, seq(30, 360, 30)*24)
ilt_x_lab   = c("Day\nof\nbirth", "1 Hour", "1 Day", "1 Week",
                "1 Month", rep("", 10),
                "1 Year")
ilt_y_break = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1, 10)

ilt2 %>%
  filter(hx != 0, !is.na(apgar5)) %>%
  ggplot(aes(x+1)) +
  geom_smooth(aes(y = hx, colour = apgar5, group = apgar5),
              size = 0.6, se = FALSE) +
  geom_smooth(aes(y = hx), colour = "red", se = FALSE, data = ilt1, size = 2) +
  scale_y_continuous("Hourly hazard of death", trans = "log10",
                     breaks = ilt_y_break,
                     labels = function (x) x*100000) +
  scale_x_continuous("Age", trans = "log10", breaks = ilt_x_break+1, labels = ilt_x_lab) +
  scale_color_continuous(guide = "none") +
  theme_minimal() +
  theme(aspect.ratio = 0.7,
        panel.grid.minor = element_blank()) -> plot_ilt_dob0510_apgar5

ggsave("./out/fig/plot_ilt_dob0510_apgar5.pdf", plot_ilt_dob0510_apgar5, width = 7, height = 5)
