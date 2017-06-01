##################
# TRANSFORM DATA #
##################

# We merge the data on births, infant- and fetal deaths from various years into
# a single table. This of course requires the variables to be harmonized across
# the years. Non-harmonized variables should be stored under different names,
# e.g. plurality91, plurality94. They can be harmonized after concatenantion.
# The matchin of variables across years is based on the variable names. If a
# variable is available in year y, but not in year z, then it will contain NAs
# for year z.

# Init --------------------------------------------------------------------

# set available memory to 100GB (important only on windows)
memory.limit(size = 100000)

library(dplyr)
library(lubridate)

# Input -------------------------------------------------------------------

# load linked birth-infant death file
load("./priv/data/01-pre_harmonized/us_ideath_1995-2010.RData")

# Merge -------------------------------------------------------------------

# merge into common data frame
ideath <- bind_rows(ideath)

# Transform ---------------------------------------------------------------

# change levels of mothers age variable to purely numerical values
# we plan to left censor at age 14 and right censor at age 49 later on
levels(ideath$age_of_mother_c03)[1]  = 14 # "Under 15" -> 14
levels(ideath$age_of_mother_c04)[1]  = 12 # "10-12" -> 12
levels(ideath$age_of_mother_c04)[39] = 50 # "50-54" -> 50

ideath %>%
  # HARMONIZATION ##############################################################
  mutate(
    # harmonize the age of mother across years
    age_of_mother_y =
      # if age of mother in years is available...
      ifelse(!is.na(age_of_mother_y),
             # use it
             age_of_mother_y,
             # otherwise see if age of mother in the definition of 2003 is available
             ifelse(!is.na(age_of_mother_c03),
                    # if yes then convert it to integer
                    as.integer(as.character(age_of_mother_c03)),
                    # otherwise use the definition from 2004 and convert it to integer
                    as.integer(as.character(age_of_mother_c04)))
      ),
    # left censor at age 14, right censor at age 49
    # (largest range available for all years)
    age_of_mother_y = ifelse(age_of_mother_y <= 14, 14, age_of_mother_y),
    age_of_mother_y = ifelse(age_of_mother_y >= 49, 49, age_of_mother_y)
  ) %>%
  # CONCEPTION COHORT ##########################################################
  mutate(
    # construct conception cohorts
    date_of_delivery_ym =
      ymd(paste(date_of_delivery_y, date_of_delivery_m, "01", sep = "-")),
    date_of_conception_ym =
      round_date(date_of_delivery_ym - weeks(gestation_at_delivery_w), unit = "month"),
    date_of_conception_y =
      year(date_of_conception_ym)
  ) %>%
  # LIFETABLE INDICATORS #######################################################
  mutate(
    # death indicator
    death =
      !is.na(age_at_death_d),
    # age at death in (completed) hours
    # we integrate additional information
    # available for the day of birth
    age_at_death_h =
      ifelse(age_at_death_d > 0,  # if death not at first day
             age_at_death_d*24, # convert age in days to hours
             # otherwise check if death happened in first hour
             # or hour 1-23 and code accordingly
             ifelse(age_at_death_c == "1-23 hours", 1, 0)),
    # age at death in (completed) weeks
    age_at_death_w =
      floor(age_at_death_d/7),
    # by definition an infant death is only registered if it occours
    # within the first 52 weeks of life, therefore we right censor at 365 days
    # age at death or censoring in (completed) days
    age_at_death_or_cens_d =
      ifelse(death, # == TRUE
             age_at_death_d,
             365),
    # age at death or censoring in (completed) hours
    age_at_death_or_cens_h =
      ifelse(death, # == TRUE
             age_at_death_h,
             365*24),
    # categorical age at death or censoring
    age_at_death_or_cens_c =
      cut(age_at_death_or_cens_h,
          # 1h, 1d, 1w, 1m, 2-12m
          breaks = c(0, 1, 24, 24*7, 24*7*4, seq(1344, 8064, 672)),
          right = FALSE, include.lowest = TRUE)
  ) %>%
  # RECODE VARIABLES ###########################################################
  mutate(
    # discrete gestational age at delivery
    gestation_at_delivery_c =
      cut(gestation_at_delivery_w,
          breaks = c(23, 28, 32, 37, 39, 41, 42, 52),
          right = FALSE, include.lowest = TRUE, # actually means to include the highest
          labels = c("extremely preterm [23,28)", "very preterm [28, 32)",
                     "moderate to late preterm [32, 37)", "early term [37, 39)",
                     "full term [39, 41)", "late term [41, 42)",
                     "post term [42, 50)")),
    # discrete age of mother
    age_of_mother_c =
      cut(age_of_mother_y,
          breaks = c(0, 14, 20, 30, 40, Inf),
          labels = c("Child", "Teenager", "[20-30)", "[30-40)", "40+"),
          right = FALSE),
    # recode education of mother
    education_of_mother =
      cut(as.integer(education_of_mother),
          breaks = c(1, 2, 10, 14, Inf),
          labels = c("No education", "Elementary school", "High school", "College"),
          right = FALSE),
    # recode race and hispanic origin of mother
    race_and_hispanic_orig_of_mother =
      cut(as.integer(race_and_hispanic_orig_of_mother),
          breaks = c(1, 6, 7, 8, Inf),
          labels = c("Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Other"),
          right = FALSE),
    # any severe congenital anomalies?
    severe_congenital_anomalies =
      (anencephalus == "Yes" |
         diaphragmatic_hernia == "Yes" |
         other_chromosomal_anomalies == "Yes" |
         renal_agensis == "Yes" |
         heart_malformations == "Yes" |
         omphalocele == "Yes" |
         hydrocephalus == "Yes" |
         microcephalus == "Yes" |
         other_congenital_anomalies == "Yes" |
         other_central_nervous_system_anomalies == "Yes"),
    # any less severe congenital anomalies?
    less_severe_congenital_anomalies =
      (cleft_lip == "Yes" |
         club_foot == "Yes" |
         downs_syndrome == "Yes" |
         other_urogenital_anomalies == "Yes" |
         polydactyly == "Yes" |
         spina_bifida == "Yes" |
         anencephalus == "Yes" |
         tracheo_esophageal_fistula == "Yes" |
         malformed_genitalia == "Yes"),
    # any least severe congenital anomalies?
    least_severe_congenital_anomalies =
      (other_musculoskeletal_anomalies == "Yes" |
         other_circulatory_respiratory_anomalies == "Yes" |
         other_gastrointestinal_anomalies == "Yes" |
         stenosis == "Yes"),
    # existence and severity of congenital anomalies
    congenital_anomalies =
      ifelse(least_severe_congenital_anomalies, 1, 0),
    congenital_anomalies =
      ifelse(less_severe_congenital_anomalies, 2, congenital_anomalies),
    congenital_anomalies =
      ifelse(severe_congenital_anomalies, 3, congenital_anomalies),
    # if no anomalies have been registred we assume none
    congenital_anomalies = ifelse(is.na(congenital_anomalies), "None", congenital_anomalies),
    congenital_anomalies = factor(congenital_anomalies,
                                  levels = 0:3, labels = c("None", "Least severe",
                                                           "Less severe", "Severe"))
  ) %>%
  # PRETTIFY
  select(
    # information on delivery
    date_of_delivery_ym,
    date_of_delivery_y,
    date_of_delivery_m,
    weekday_of_delivery,
    method_of_delivery,
    prolonged_labor,
    # information on child
    sex,
    gestation_at_delivery_w,
    gestation_at_delivery_c,
    date_of_conception_ym,
    date_of_conception_y,
    birthweight_g,
    apgar5,
    plurality,
    birth_injury,
    # congenital anomalies
    congenital_anomalies,
    anencephalus,
    spina_bifida,
    hydrocephalus,
    microcephalus,
    other_central_nervous_system_anomalies,
    heart_malformations,
    other_circulatory_respiratory_anomalies,
    stenosis,
    tracheo_esophageal_fistula,
    omphalocele,
    other_gastrointestinal_anomalies,
    malformed_genitalia,
    renal_agensis,
    other_urogenital_anomalies,
    cleft_lip,
    polydactyly,
    club_foot,
    diaphragmatic_hernia,
    other_musculoskeletal_anomalies,
    downs_syndrome,
    other_chromosomal_anomalies,
    other_congenital_anomalies,
    # information on mother
    age_of_mother_y,
    age_of_mother_c,
    resident_status_of_mother,
    race_and_hispanic_orig_of_mother,
    martial_status_of_mother,
    education_of_mother,
    tobacco_use_during_pregnancy,
    alcohol_use_during_pregnancy,
    # information on medical care
    time_prenatal_care_began,
    # information on survival
    death,
    age_at_death_w,
    age_at_death_d,
    age_at_death_h,
    age_at_death_or_cens_d,
    age_at_death_or_cens_h,
    age_at_death_or_cens_c
  ) %>%
  arrange(date_of_delivery_ym, age_at_death_h) -> ideath

ideath %>% filter(death == TRUE) %>% summary
ideath %>% filter(death == FALSE) %>% summary

save(ideath, file = "./priv/data/02-harmonized/ideath.RData")

haven::write_dta(
  ideath %>% rename(othr_cntrl_nervs_sys_anomalies = other_central_nervous_system_anomalies,
                    othr_circ_resp_anomalies = other_circulatory_respiratory_anomalies),
  path = "./priv/data/02-harmonized/ideath.dta")
