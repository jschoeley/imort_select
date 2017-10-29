############################################
# WHICH CONGENITAL DISEASES ARE DEADLIEST? #
############################################

# Init --------------------------------------------------------------------

library(tidyverse)

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

ideath_sub <-
  ideath %>%
  # coding of congenital disorders changes after 2004 and I haven't
  # transcribed the new coding yet. the variable with the old coding
  # is still available post 2004 but has lots of missings.
  filter(date_of_delivery_y == 2004) %>%
  select(one_of(congenital_disorder_varnames),
         sex, plurality, birthweight_c, gestation_at_delivery_c,
         martial_status_of_mother, race_and_hispanic_orig_of_mother,
         education_of_mother, death) %>%
  mutate(plurality = factor(plurality, ordered = FALSE))

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
      sex + plurality + gestation_at_delivery_c + birthweight_c +
      martial_status_of_mother + race_and_hispanic_orig_of_mother +
      education_of_mother,
      family = binomial(link=logit), data = ideath_sub)

# ranking in odds ratio of infant death by congenital disorder
odds_ratio_ideath <-
  fit_ideath %>%
  broom::tidy() %>%
  mutate(odds_ratio = exp(estimate)) %>%
  filter(grepl(pattern = paste0(congenital_disorder_varnames, collapse = "|"), term)) %>%
  mutate(condition = gsub("Yes", "", term)) %>% select(-term) %>%
  arrange(desc(odds_ratio))


# IMR calculation ---------------------------------------------------------

# probablity of infant death by type of congenital condition
imr_ideath = data.frame(condition = congenital_disorder_varnames, N = NA, D = NA, imr = NA)
for (i in congenital_disorder_varnames) {
  exposed = ideath_sub[ideath_sub[,i] == "Yes" & !is.na(ideath_sub[,i]),]
  N = nrow(exposed)
  D = nrow(exposed[exposed$death == TRUE,])
  imr_ideath[imr_ideath$condition == i, "N"] = N
  imr_ideath[imr_ideath$condition == i, "D"] = D
  imr_ideath[imr_ideath$condition == i, "imr"] = D/N
  rm(exposed, N, D)
}

# Summary -----------------------------------------------------------------

# a summary table
ideath_by_condition <-
  imr_ideath %>%
  mutate(rel_to_total_imr = imr/(sum(ideath_sub$death)/nrow(ideath_sub))) %>%
  arrange(desc(imr)) %>%
  left_join(odds_ratio_ideath, by = "condition")