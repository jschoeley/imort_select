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

library(tidyverse)
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

ideath <-
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
    death = !is.na(age_at_death_d),
    # age at death in (completed) hours
    # we integrate additional information
    # available for the day of birth
    age_at_death_h =
      ifelse(age_at_death_d > 0,  # if death not at first day
             age_at_death_d*24, # convert age in days to hours
             # otherwise check if death happened in first hour
             # or hour 1-23 and code accordingly
             ifelse(age_at_death_c == "1-23 hours", 1, 0)),
    # so now we have interval censored data on the age at death in hours,
    # to deal with it we add the width of the age interval of death [x, x+nx)
    age_at_death_h_width = 24,
    age_at_death_h_width = ifelse(age_at_death_h == 0, 1, age_at_death_h_width),
    age_at_death_h_width = ifelse(age_at_death_h == 1, 23, age_at_death_h_width),
    # by definition an infant death is only registered if it occours within
    # the first 52 weeks of life, therefore we right censor at 365 days or
    # 8767 hours (Gregorian solar calendar: 1 Month = 30.44 days = 730.56 hours)
    survtime_d = ifelse(death == TRUE, age_at_death_d, 365),
    # age at death or censoring in (completed) hours
    survtime_h = ifelse(death == TRUE, age_at_death_h, 8767),
    survtime_h_width = 24,
    survtime_h_width = ifelse(survtime_h == 0, 1, survtime_h_width),
    survtime_h_width = ifelse(survtime_h == 1, 23, survtime_h_width),
    # grouped age at death or censoring, categorical
    survtime_c =
      cut(survtime_h,
          # 1h-1d-1w-4w-1m-12m in hours
          breaks = c(0, 1, seq(24, 144, 24), 168, 336, 504, ceiling(30.44*24*(1:12))),
          labels = c("[Birth, 1h)",
                     "[1h, 1d)",
                     "[1d, 2d)", "[2d, 3d)", "[3d, 4d)", "[4d, 5d)", "[6d, 1w)",
                     "[1w, 2w)", "[2w, 3w)", "[3w, 1m)",
                     "[1m, 2m)", "[2m, 3m)", "[3m, 4m)", "[4m, 5m)", "[5m, 6m)",
                     "[6m, 7m)", "[7m, 8m)", "[8m, 9m)", "[9m, 10m)", "[10m, 11m)",
                     "[11m, 1y)", "[1y+"),
          # actually means to include the highest
          right = FALSE, include.lowest = TRUE)
  ) %>%
  # RECODE VARIABLES ###########################################################
  mutate(
    # discrete gestational age at delivery
    gestation_at_delivery_c =
      cut(gestation_at_delivery_w,
          breaks = c(0, 28, 32, 37, 39, 41, 42, 52),
          labels = c("Extremely preterm <28", "Very preterm [28, 32)",
                     "Moderate to late preterm [32, 37)", "Early term [37, 39)",
                     "Full term [39, 41)", "Late term [41, 42)",
                     "Post term [42, 50)"),
          # actually means to include the highest
          right = FALSE, include.lowest = TRUE),
    # discrete birthweight
    birthweight_c =
      cut(birthweight_g,
          breaks = c(0, 1000, 1500, 2500, 4200, Inf),
          label = c("Extremely low", "Very low", "Low", "Regular", "High"),
          right = FALSE),
    # discrete age of mother
    age_of_mother_c =
      cut(age_of_mother_y,
          breaks = c(0, 16, 20, 30, 40, Inf),
          labels = c("<16", "[16, 20)", "[20, 30)", "[30, 40)", "40+"),
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
    congenital_anomalies = ifelse(least_severe_congenital_anomalies, 1, 0),
    congenital_anomalies = ifelse(less_severe_congenital_anomalies, 2, congenital_anomalies),
    congenital_anomalies = ifelse(severe_congenital_anomalies, 3, congenital_anomalies),
    # if no anomalies have been registred (even due to NA) we assume none
    congenital_anomalies = ifelse(is.na(congenital_anomalies), 0, congenital_anomalies),
    congenital_anomalies = factor(congenital_anomalies,
                                  levels = 0:3,
                                  labels = c("None", "Least severe", "Less severe", "Severe"))
  ) %>%
  # SET REFERENCE CATEGORIES FOR FACTOR VARIABLES
  mutate(
    method_of_delivery = relevel(as.factor(method_of_delivery), "Vaginal (no previous C-section)"),
    gestation_at_delivery_c = relevel(gestation_at_delivery_c, "Full term [39, 41)"),
    birthweight_c = relevel(birthweight_c, "Regular"),
    age_of_mother_c = relevel(age_of_mother_c, "[20, 30)"),
    education_of_mother = relevel(education_of_mother, "High school"),
    race_and_hispanic_orig_of_mother = relevel(race_and_hispanic_orig_of_mother, "Non-Hispanic White")
  ) %>%
  mutate_at(
    vars(
      prolonged_labor,
      birth_injury,
      tobacco_use_during_pregnancy,
      alcohol_use_during_pregnancy,
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
      other_congenital_anomalies
    ),
    funs(relevel(., "No"))
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
    birthweight_c,
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
    age_at_death_d,
    age_at_death_h,
    age_at_death_h_width,
    survtime_h,
    survtime_h_width,
    survtime_c
  ) %>%
  arrange(date_of_delivery_ym, age_at_death_h)

#ideath %>% filter(death == TRUE) %>% summary
#ideath %>% filter(death == FALSE) %>% summary

# Save --------------------------------------------------------------------

# save all of it
save(ideath, file = "./priv/data/02-harmonized/ideath.RData")

# save a small subset
ideath2004_small <-
  ideath %>%
  filter(date_of_delivery_y == 2004) %>%
  select(sex, age_of_mother_c, death, survtime_h, survtime_h_width,
         birthweight_c, gestation_at_delivery_c,
         congenital_anomalies, apgar5, education_of_mother,
         race_and_hispanic_orig_of_mother, martial_status_of_mother)
save(ideath2004_small, file = "./priv/data/02-harmonized/ideath2004_small.RData")
haven::write_dta(ideath2004_small, path = "./priv/data/02-harmonized/ideath2004_small.dta")
