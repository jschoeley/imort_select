##################
# TRANSFORM DATA #
##################

# We harmonize variables across years, and recode/discretize most variables.
# Non-harmonized variables should be stored under different names, e.g.
# plurality91, plurality94. They can be harmonized after concatenantion. The
# matching of variables across years is based on the variable names. If a
# variable is available in year y, but not in year z, then it will contain NAs
# for year z.

# Init --------------------------------------------------------------------

# set available memory to 100GB (important only on windows)
memory.limit(size = 100000)

library(tidyverse)
library(lubridate)

# Input -------------------------------------------------------------------

# load linked birth-infant death file
load('./priv/data/01-pre_harmonized/us_ideath_1995-2010.RData')

# Merge -------------------------------------------------------------------

# merge different years into common data frame
ideath <- bind_rows(ideath)

# Transform ---------------------------------------------------------------

# change levels of mothers age variable to purely numerical values
# we plan to left censor at age 14 and right censor at age 49 later on
levels(ideath$age_of_mother_c03)[1]  = 14 # 'Under 15' -> 14
levels(ideath$age_of_mother_c04)[1]  = 12 # '10-12' -> 12
levels(ideath$age_of_mother_c04)[39] = 50 # '50-54' -> 50

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
      ymd(paste(date_of_delivery_y, date_of_delivery_m, '01', sep = '-')),
    # discrete gestational age at delivery
    gestation_at_delivery_c7 =
      cut(gestation_at_delivery_w,
          breaks = c(0, 28, 32, 37, 39, 41, 42, Inf),
          labels = c('Extremely preterm <28w', 'Very preterm [28,32)w',
                     'Moderate to late preterm [32,37)w', 'Early term [37,39)w',
                     'Full term [39,41)w', 'Late term [41,42)w',
                     'Post term [42,50)w'),
          right = FALSE),
    gestation_at_delivery_c4 =
      cut(gestation_at_delivery_w,
          breaks = c(0, 28, 32, 37, Inf),
          labels = c('Extremely preterm <28w', 'Very preterm [28,32)w',
                     'Moderate to late preterm [32,37)w', 'Term 37w+'),
          right = FALSE),
    # discrete birthweight
    birthweight_c5 =
      cut(birthweight_g,
          breaks = c(0, 1000, 1500, 2500, 4200, Inf),
          label = c('Extremely low [0,1000)g', 'Very low [1000,1500)g',
                    'Low [1500,2500)g', 'Regular [2500,4200)g', 'High 4200g+'),
          right = FALSE),
    birthweight_c3 =
      cut(birthweight_g,
          breaks = c(0, 1500, 2500, Inf),
          label = c('Very low [0,1500)g', 'Low [1500,2500)g', 'Regular 2500g+'),
          right = FALSE),
    # discrete apgar score
    apgar5_c3 =
      cut(apgar5,
          breaks = c(0, 5, 9, Inf),
          label = c('Very low [0,5)', 'Low [5,9)', 'Regular 9+'),
          right = FALSE),
    # discrete age of mother
    age_of_mother_c5 =
      cut(age_of_mother_y,
          breaks = c(0, 16, 20, 30, 40, Inf),
          labels = c('<16y', '[16, 20)y', '[20, 30)y', '[30, 40)y', '40y+'),
          right = FALSE),
    age_of_mother_c3 =
      cut(age_of_mother_y,
          breaks = c(0, 16, 40, Inf),
          labels = c('<16y', '[16, 40)y', '40y+'),
          right = FALSE),
    # recode education of mother
    education_of_mother_c4 =
      cut(as.integer(education_of_mother),
          breaks = c(1, 2, 10, 14, Inf),
          labels = c('No education', 'Elementary school', 'High school', 'College'),
          right = FALSE),
    education_of_mother_c2 =
      cut(as.integer(education_of_mother),
          breaks = c(1, 10, Inf),
          labels = c('No or elementary education', 'High school or College'),
          right = FALSE),
    # recode race and hispanic origin of mother
    race_and_hispanic_orig_of_mother_c4 =
      cut(as.integer(race_and_hispanic_orig_of_mother),
          breaks = c(1, 6, 7, 8, Inf),
          labels = c('Hispanic', 'Non-Hispanic White', 'Non-Hispanic Black', 'Other'),
          right = FALSE),
    race_and_hispanic_orig_of_mother_c2 =
      recode_factor(race_and_hispanic_orig_of_mother_c4,
                    'Hispanic' = 'Other',
                    'Non-Hispanic White' = 'Other',
                    'Non-Hispanic Black' = 'Non-Hispanic Black',
                    'Other' = 'Other'),
    # recode plurality
    plurality_c4 =
      recode_factor(plurality,
                    'Single' = 'Single', 'Twin' = 'Twin', 'Triplet' = 'Triplet',
                    'Quadruplet' = 'Quadruplet or higher',
                    'Quintruplet or higher' = 'Quadruplet or higher'),
    plurality_c2 =
      recode_factor(plurality_c4,
                    'Single' = 'Single',
                    'Twin' = 'Twin or higher',
                    'Triplet' = 'Twin or higher',
                    'Quadruplet or higher' = 'Twin or higher'),
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
      (anencephalus == 'Yes' |
         diaphragmatic_hernia == 'Yes' |
         other_chromosomal_anomalies == 'Yes' |
         renal_agensis == 'Yes' |
         heart_malformations == 'Yes' |
         omphalocele == 'Yes' |
         hydrocephalus == 'Yes' |
         microcephalus == 'Yes' |
         other_congenital_anomalies == 'Yes' |
         other_central_nervous_system_anomalies == 'Yes'),
    # any less severe congenital anomalies?
    less_severe_congenital_anomalies =
      (cleft_lip == 'Yes' |
         club_foot == 'Yes' |
         downs_syndrome == 'Yes' |
         other_urogenital_anomalies == 'Yes' |
         polydactyly == 'Yes' |
         spina_bifida == 'Yes' |
         anencephalus == 'Yes' |
         tracheo_esophageal_fistula == 'Yes' |
         malformed_genitalia == 'Yes'),
    # any least severe congenital anomalies?
    least_severe_congenital_anomalies =
      (other_musculoskeletal_anomalies == 'Yes' |
         other_circulatory_respiratory_anomalies == 'Yes' |
         other_gastrointestinal_anomalies == 'Yes' |
         stenosis == 'Yes'),
    # existence and severity of congenital anomalies
    congenital_anomalies_c4 = ifelse(least_severe_congenital_anomalies, 1, 0),
    congenital_anomalies_c4 = ifelse(less_severe_congenital_anomalies, 2, congenital_anomalies_c4),
    congenital_anomalies_c4 = ifelse(severe_congenital_anomalies, 3, congenital_anomalies_c4),
    # if no anomalies have been registred (even due to NA) we assume none
    congenital_anomalies_c4 = ifelse(is.na(congenital_anomalies_c4), 0, congenital_anomalies_c4),
    congenital_anomalies_c4 = factor(congenital_anomalies_c4,
                                     levels = 0:3,
                                     labels = c('None', 'Least severe', 'Less severe', 'Severe')),
    congenital_anomalies_c3 = recode_factor(congenital_anomalies_c4,
                                            'None' = 'None',
                                            'Least severe' = 'Less severe',
                                            'Less severe' = 'Less severe',
                                            'Severe' = 'Severe')
  ) %>%
  # SET REFERENCE CATEGORIES FOR FACTOR VARIABLES
  mutate(
    method_of_delivery =
      relevel(as.factor(method_of_delivery), 'Vaginal (no previous C-section)'),
    gestation_at_delivery_c7 =
      relevel(gestation_at_delivery_c7, 'Full term [39,41)w'),
    gestation_at_delivery_c4 =
      relevel(gestation_at_delivery_c4, 'Term 37w+'),
    birthweight_c5 =
      relevel(birthweight_c5, 'Regular [2500,4200)g'),
    birthweight_c3 =
      relevel(birthweight_c3, 'Regular 2500g+'),
    apgar5_c3 =
      relevel(apgar5_c3, 'Regular 9+'),
    age_of_mother_c5 =
      relevel(age_of_mother_c5, '[20, 30)y'),
    age_of_mother_c3 =
      relevel(age_of_mother_c3, '[16, 40)y'),
    education_of_mother_c4 =
      relevel(education_of_mother_c4, 'College'),
    education_of_mother_c2 =
      relevel(education_of_mother_c2, 'High school or College'),
    race_and_hispanic_orig_of_mother_c4 =
      relevel(race_and_hispanic_orig_of_mother_c4, 'Non-Hispanic White'),
    race_and_hispanic_orig_of_mother_c2 =
      relevel(race_and_hispanic_orig_of_mother_c2, 'Other')
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
    funs(relevel(., 'No'))
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
    gestation_at_delivery_c7, gestation_at_delivery_c4,
    birthweight_g, birthweight_c3, birthweight_c5,
    apgar5, apgar5_c3,
    plurality, plurality_c4, plurality_c2,
    birth_injury,
    congenital_anomalies_c4, congenital_anomalies_c3,
    # information on mother
    age_of_mother_y, age_of_mother_c5, age_of_mother_c3,
    resident_status_of_mother,
    race_and_hispanic_orig_of_mother_c4, race_and_hispanic_orig_of_mother_c2,
    martial_status_of_mother,
    education_of_mother_c4, education_of_mother_c2,
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
save(ideath, file = paste0('./priv/data/02-harmonized/', Sys.Date(), '-ideath.RData'))
