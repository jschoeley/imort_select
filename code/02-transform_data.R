##################
# TRANSFORM DATA #
##################

# We harmonize variables across years, and recode some variables. Non-harmonized
# variables should be stored under different names, e.g. plurality91,
# plurality94. They can be harmonized after concatenantion. The matching of
# variables across years is based on the variable names. If a variable is
# available in year y, but not in year z, then it will contain NAs for year z.

# Init --------------------------------------------------------------------

# set available memory to 100GB (important only on windows)
memory.limit(size = 100000)

library(tidyverse)
library(lubridate)

# Input -------------------------------------------------------------------

# load linked birth-infant death file
load("./priv/data/01-pre_harmonized/us_ideath_1995-2010.RData")

# Merge -------------------------------------------------------------------

# merge different years into common data frame
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
  # RECODE VARIABLES ###########################################################
  mutate(
    # date of delivery
    date_of_delivery_ym =
      ymd(paste(date_of_delivery_y, date_of_delivery_m, "01", sep = "-")),
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
    # recode plurality
    plurality =
      recode_factor(plurality,
                    "Single" = "Single", "Twin" = "Twin", "Triplet" = "Triplet",
                    "Quadruplet" = "Quadruplet or higher",
                    "Quintruplet or higher" = "Quadruplet or higher"),
    # merge alcohol and tobacco use during pregnancy
    alcto_use_during_pregnancy = NA,
    alcto_use_during_pregnancy = ifelse(
      tobacco_use_during_pregnancy == 'Yes' |
        alcohol_use_during_pregnancy == 'Yes',
      'Yes', alcto_use_during_pregnancy),
    alcto_use_during_pregnancy = ifelse(
      tobacco_use_during_pregnancy == 'No' &
        alcohol_use_during_pregnancy == 'No',
      'No', alcto_use_during_pregnancy),
    alcto_use_during_pregnancy = as.factor(alcto_use_during_pregnancy),
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
    method_of_delivery =
      relevel(as.factor(method_of_delivery), "Vaginal (no previous C-section)"),
    gestation_at_delivery_c =
      relevel(gestation_at_delivery_c, "Full term [39, 41)"),
    birthweight_c =
      relevel(birthweight_c, "Regular"),
    age_of_mother_c =
      relevel(age_of_mother_c, "[20, 30)"),
    education_of_mother =
      relevel(education_of_mother, "High school"),
    race_and_hispanic_orig_of_mother =
      relevel(race_and_hispanic_orig_of_mother, "Non-Hispanic White")
  ) %>%
  mutate_at(
    vars(
      prolonged_labor,
      birth_injury,
      tobacco_use_during_pregnancy,
      alcohol_use_during_pregnancy,
      alcto_use_during_pregnancy,
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
    birthweight_g,
    birthweight_c,
    apgar5,
    plurality,
    birth_injury,
    congenital_anomalies,
    # information on mother
    age_of_mother_y,
    age_of_mother_c,
    resident_status_of_mother,
    race_and_hispanic_orig_of_mother,
    martial_status_of_mother,
    education_of_mother,
    alcto_use_during_pregnancy,
    # information on medical care
    time_prenatal_care_began,
    # information on survival
    age_at_death_d,
    age_at_death_c
  ) %>%
  arrange(date_of_delivery_ym, age_at_death_d)

#ideath %>% filter(death == TRUE) %>% summary
#ideath %>% filter(death == FALSE) %>% summary

# Save --------------------------------------------------------------------

# save all of it
save(ideath, file = paste0("./priv/data/02-harmonized/", Sys.Date(), "-ideath.RData"))