#################################################
# Which congenital disorders are the deadliest? #
#################################################

# Init --------------------------------------------------------------------

library(dplyr)
library(effects)

# Input -------------------------------------------------------------------

load("./priv/data/02-harmonized/ideath.RData")

# Transform ---------------------------------------------------------------

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

ideath %>%
  filter(date_of_delivery_y == 2004) %>%
  select(one_of(congenital_disorder_varnames),
         sex, plurality, gestation_at_delivery_c, death) %>%
  mutate_at(congenital_disorder_varnames, relevel, ref = "No") %>%
  mutate(plurality = as.factor(plurality)) -> ideath_sub

summary(ideath_sub)
rm(ideath)

# Logit regression --------------------------------------------------------

fit_ideath <-
  glm(death ~
      anencephalus + spina_bifida + hydrocephalus + microcephalus +
      other_central_nervous_system_anomalies + heart_malformations +
      other_circulatory_respiratory_anomalies + stenosis + tracheo_esophageal_fistula +
      omphalocele + other_gastrointestinal_anomalies + malformed_genitalia + renal_agensis +
      other_urogenital_anomalies + cleft_lip + polydactyly + club_foot + diaphragmatic_hernia +
      other_musculoskeletal_anomalies + downs_syndrome + other_chromosomal_anomalies +
      other_congenital_anomalies +
      sex + plurality + gestation_at_delivery_c,
      family = binomial, data = ideath_sub)

marginal_fx_ideath <- effect(congenital_disorder_varnames, fit_ideath)

names(ideath_sub) <- abbreviate(names(ideath_sub), minlength = 10)

haven::write_dta(ideath_sub, path = "./priv/data/02-harmonized/ideath.dta")
