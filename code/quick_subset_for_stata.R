library(dplyr)

load("./priv/data/02-harmonized/ideath.RData")

congenital_disorder_varnames <-
  c("anencephalus",
    "spina_bifida",
    "hydrocephalus",
    "microcephalus",
    "other_central_nervous_system_anomalies",
    "heart_malformations",
    "other_circulatory_respiratory_anomalies",
    "stenosis",
    "tracheo_esophageal_fistula",
    "omphalocele",
    "other_gastrointestinal_anomalies",
    "malformed_genitalia",
    "renal_agensis",
    "other_urogenital_anomalies",
    "cleft_lip",
    "polydactyly",
    "club_foot",
    "diaphragmatic_hernia",
    "other_musculoskeletal_anomalies",
    "downs_syndrome",
    "other_chromosomal_anomalies",
    "other_congenital_anomalies")

ideath2004 <-
  ideath %>%
  filter(date_of_delivery_y == 2004) %>%
  select(sex, age_of_mother_y, death, age_at_death_or_cens_h,
         birthweight_g, congenital_anomalies, apgar5, education_of_mother,
         race_and_hispanic_orig_of_mother) %>%
  mutate(
    age_of_mother_c =
      cut(age_of_mother_y,
          breaks = c(0, 14, 20, 30, 40, Inf),
          labels = c("Child", "Teenager", "[20-30)", "[30-40)", "40+"),
          right = FALSE)
  ) %>% select(-age_of_mother_y)

haven::write_dta(ideath2004, path = "./priv/data/02-harmonized/ideath2004.dta")
