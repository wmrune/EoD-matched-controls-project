rm(list = ls())

library(DBI)
library(RClickhouse)
library(askpass)
library(tidyverse)
library(dbplyr)
library(readxl)
library(fixest)

 
# ---- Connection 
con <- dbConnect(clickhouse(),
                 host = "h1cogbase01.nvs.ki.se",
                 port = 9000,
                 user = "william_rune",
                 password = askpass::askpass("Please enter clickhouse password: ")
)
 
# ---- Inputs ----
diag_dementia <- c("F00", "F01", "F02", "F03", "F051", "G30", "G31")
diag_mci <- c("F067")

#Here, I create a table for ICD-codes for the report
library(knitr)
library(kableExtra)

icd_codes <- data.frame(
  ICD10 = c(
    "F00", "F01", "F02", "F03", "F051",
    "G30", "G31", "F067"
  ),
  Diagnos = c(
    "Dementia in Alzheimer's Disease",
    "Vascular Dementia",
    "Dementia in other specified diseases classified elsewhere",
    "Unspecified Dementia",
    "Delirium superimposed on dementia",
    "Alzheimer's Disease",
    "Other specified degenerative diseases of nervous system",
    "Mild cognitive impairment (MCI)"
  )
)

kable(icd_codes, format = "html") %>%
  kable_styling() %>%
  pack_rows("Dementia codes", 1, 7) %>%
  pack_rows("MCI code", 8, 8)



# Filter for diagnoses (dementia and MCI)
rx_vec <- paste0("^", c(diag_dementia, diag_mci))
rx_sql <- paste0("'", paste(rx_vec, collapse = "|"), "'")
diag_pred_sql <- paste0("match(diagnos, ", rx_sql, ")")
 
# MCI filter for diag_group
mci_rx_sql <- paste0("'", paste(paste0("^", diag_mci), collapse = "|"), "'")
mci_pred_sql <- paste0("match(diagnos, ", mci_rx_sql, ")")

# ---- Tables ----
diagnoses <- tbl(con, in_schema("cognet", "diagnoses_long")) # diagnoses table from par and regional pc data
scb_lisa <- tbl(con, in_schema("cognet", "scb_lisa"))
pop_desc_tbl <- tbl(con, in_schema("cognet", "scb_population")) # Population description table age
cod <- tbl(con, in_schema("cognet", "cod")) # Cause of death
case_control <- tbl(con, in_schema("cognet", "case_control"))



# Skapa studiepopulation
# NOTE (24 Aug): fixed a bug where the age<65 filter was applied per diagnosis
# row (before taking the earliest date) instead of to the true first-diagnosis
# age. That could silently compute the index from a *later* dementia code if
# an earlier one fell outside the age window, understating true diagnosis age.
# Now: take each person's first Dementia-coded event (MCI excluded, per
# protocol), THEN filter on age at that true first event.
t_pop <- Sys.time()  # timing: this scans diagnoses_long (~303M rows), likely the slowest single step
pop <- diagnoses %>%
  select(lopnr, indatuma, setting, diagnos) %>%
  left_join(pop_desc_tbl |> select(lopnr, fodelsear), by = "lopnr") |>
  filter(sql(diag_pred_sql)) %>%
  mutate(
    lopnr = sql("CAST(lopnr AS Int64)"),
    diag_group = sql(paste0("if(", mci_pred_sql, ", 'MCI', 'Dementia')")),
    setting = if_else(setting != "pc", "sc/inpat", setting),
    year_diag = year(indatuma)
  ) |>
  collect() |>
  filter(diag_group == "Dementia") |>            # study scope is dementia; MCI excluded (protocol §5.1)
  group_by(lopnr) |>
  summarise(
    index = min(indatuma),
    age_index = year(min(indatuma)) - min(fodelsear)
  ) |>
  ungroup() |>
  filter(
    age_index >= 40, age_index < 65,              # EoD definition: diagnosis age 40-64 (protocol §5.1)
    year(index) >= 2015, year(index) <= 2023       # study ascertainment period
  )
cat("pop query took:", round(difftime(Sys.time(), t_pop, units = "mins"), 2), "min\n")

ggplot(pop |> filter(age_index > 0), aes(age_index)) +
  geom_histogram(bins=30)

#Tar ut kontoller
t_cc <- Sys.time()  # timing: case_control is ~1.69M rows, should be quick
controls <- case_control |>
  collect() |>
  filter(lopnr_fall %in% pop$lopnr) |>
  mutate(cc_indexdatum = as.Date(indexdatum))
cat("case_control query took:", round(difftime(Sys.time(), t_cc, units = "mins"), 2), "min\n")

# Safeguard against immortal-time/survivorship bias (cognet.md §13.2): a
# control was only validated as alive & dementia-free up to case_control's
# OWN (sometimes wrong) index date. Since we use our own diagnoses_long-based
# index instead, drop any pair whose original case_control index postdates it.
controls <- controls |>
  left_join(pop |> select(lopnr_fall = lopnr, study_index = index), by = "lopnr_fall") |>
  filter(cc_indexdatum <= study_index)

# Cap at 2 controls per case (protocol §5.2); lowest lopnr_kontroll kept for
# reproducibility when more than 2 matches are available.
controls <- controls |>
  group_by(lopnr_fall) |>
  slice_min(lopnr_kontroll, n = 2, with_ties = FALSE) |>
  ungroup()

#skapa kontrollpopulation
controls <- controls %>%
  select(lopnr = lopnr_kontroll, lopnr_fall) |>
  left_join(pop_desc_tbl |> select(lopnr, fodelsear) |> collect(), by = "lopnr") |>
  distinct(lopnr, lopnr_fall, fodelsear)


# Totala studiepopulationen:
# ett unikt lopnr per individ
tot_lopnr <- bind_rows(
  pop |> select(lopnr),
  controls |> select(lopnr)
) |>
  distinct(lopnr)

#Beräkna inkomst per år innan och efter index
# FIX (24 Aug): the old version filtered scb_lisa by year only, collect()ed
# the WHOLE table for those years (~42M rows, i.e. every person in Sweden,
# not just this cohort), and only THEN joined down to tot_lopnr. The
# commented-out `filter(lopnr %in% tot_lopnr$lopnr)` below was the original
# attempt to avoid that -- it doesn't push down to SQL through dbplyr against
# ClickHouse, so it silently did nothing useful and was replaced with the
# post-collect join instead of being fixed. This is almost certainly the
# "43 million rows" you ran into before. Fixed the same way the CCI section
# already does it further down: batch lopnr into IN (...) lists and filter
# server-side, so only cohort rows ever cross the network.
t_lisa <- Sys.time()
lopnr_vec_income <- tot_lopnr |> pull(lopnr) |> unique()
lisa_batch_size <- 5000
lopnr_batches_income <- split(lopnr_vec_income, ceiling(seq_along(lopnr_vec_income) / lisa_batch_size))

income_list <- lapply(lopnr_batches_income, function(ids) {
  ids_sql <- paste(ids, collapse = ",")
  scb_lisa |>
    filter(sql(paste0("lopnr IN (", ids_sql, ")"))) |>
    filter(year > 2000, year <= 2024) |>              # widened: LISA runs to 2024, needed for +3y follow-up
    select(lopnr, year, raks_forvink, inkkalla1_forv, sun2020niva, sun2000niva) |>  # both SUN vintages: sun2020niva alone had ~47% NA (only populated from its introduction year, per cognet.md), sun2000niva is the fallback
    collect()
})

income_raw <- bind_rows(income_list) |>
  mutate(
    income = ifelse(year > 2021, inkkalla1_forv, raks_forvink)
  )
cat("scb_lisa query took:", round(difftime(Sys.time(), t_lisa, units = "mins"), 2), "min\n")
cat("income_raw rows pulled:", nrow(income_raw), "\n")

# Outlier handling (protocol §8): winsorise at the 99th percentile of the
# pooled case+control income distribution, rather than an arbitrary round-number
# cutoff. Values above the cap are pulled down to it, not dropped, so nobody
# loses a person-year of follow-up over one extreme value.
income_cap <- quantile(income_raw$income, 0.99, na.rm = TRUE)
income_raw <- income_raw |>
  mutate(income = pmin(income, income_cap))


# Skapa totalpopoulation med lopnr, exposure, index, ålder
pop_tot <- bind_rows(
  pop |> 
    mutate(dementia = 1) |> 
    select(lopnr, dementia, index, age_index),

  controls |> 
    left_join(
      pop |> select(lopnr, index),
      by = c("lopnr_fall" = "lopnr")
    ) |> 
    mutate(
      dementia = 0,
      age_index = year(index) - fodelsear
    ) |> 
    select(lopnr, dementia, index, age_index)
) %>%
  distinct(lopnr, .keep_all = TRUE)

# Joina inkomst med totalpopulation
income_by_year <- income_raw |> left_join(pop_tot, by= "lopnr")

# Skapa år förhållande till index
# Alignment rule (protocol §6): LISA income is annual, diagnosdatum is not, so
# a mid-year cutoff decides which calendar year is the "last full pre-year".
# Diagnosis before July -> that calendar year already counts as post; after
# June -> that year is still counted as pre, post starts the following year.
year_to_index <- income_by_year |>
  mutate(
    diagnosis_month = month(index),
    anchor_year = if_else(diagnosis_month <= 6, year(index), year(index) + 1),
    t = year - anchor_year,      # relative year: t = 0 is the first full post-diagnosis year
    post = if_else(t >= 0, 1, 0)
  ) |>
  filter(t >= -5, t <= 3)         # analysis window per protocol §3

# Första plot (this is the RQ2 descriptive check: does the gap already show
# before t = 0, i.e. before diagnosis?)
p1 <- year_to_index |>
  mutate(group_label = if_else(dementia == 1, "EoD", "Control")) |>
  group_by(group_label, t) |>
  summarise(
    mean_income = mean(income, na.rm = T),
    .groups = "drop"
  ) |>
  ggplot(aes(x = t, y = mean_income, group = group_label, colour = group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  labs(
    x = "Years relative to diagnosis",
    y = "Mean income (SEK, LISA units)",
    colour = NULL,
    title = "Income trajectory, EoD vs. matched controls"
  ) +
  theme_minimal()

p1

#Charlson Comorbidity Index
# FLAGGED (24 Aug): this scans + aggregates ALL of diagnoses_long (~303M rows,
# every person in the database) with no lopnr or date filter, then collects
# the result -- this is very likely the real "43 million rows" incident, or a
# bigger version of it. `diag` is not referenced anywhere later in this
# script; the CCI computation below uses the properly-filtered `patients`
# object instead. Commented out rather than deleted -- uncomment only if you
# remember why this was here and actually need it, and add a lopnr/date
# filter first if you do.
# diag <- diagnoses |>
#   group_by(lopnr, indatuma) |>
#   summarise(
#     icd10 = sql("arrayStringConcat(arraySort(groupUniqArray(diagnos)), ' ')"),
#     .groups = "drop"
#   ) |>
#   collect()

# REMOVED (24 Aug): this used to re-join `index` onto pop_tot a second time,
# by plain lopnr. pop_tot already has a correct `index` column from when it
# was first built above (cases get their own index; controls get their
# case's index via the lopnr_fall join) -- joining again here collided the
# column into index.x/index.y (hence "column index doesn't exist" below), and
# even without the naming clash it would have been wrong: pop's lopnr only
# covers cases, so matching it against pop_tot's lopnr (cases + controls)
# would silently NA out every control's index anyway.

# beräkna CCI på det jag får kvar enligt GitHub-koden, tabell med lopnr, indatuma och spacesep diagnoser

## Help

# Import patient file in long format. The patient file should contain three columns:
# One column with patient/ID numbers (blir lopnr),
# one column with ICD codes (blir diagnos), 
# and one column with dates (indatuma). As of version 4 of the script, it is recommended to use the OUTDATE to look for ICD codes, and the INDATE to define disease onset. 
# The codes should be a string variable with different ICD codes separated by a whitespace. 
# ICD10 codes should not have a dot.
# The date column should be formatted as yyyymmdd.

# Required packages
library(dplyr)

#reformaterar datumraden utan bindestreck
# Import patient file and call it "patients"
# Hämta diagnoser för studiepopulationen
# Använd en temporär tabell i ClickHouse i stället för
# att skapa en mycket lång SQL-fråga med IN (...).
# Hämta diagnoser i mindre batchar.
# Ingen skrivbehörighet till ClickHouse behövs.

lopnr_vec <- tot_lopnr |>
  pull(lopnr) |>
  unique()

batch_size <- 5000

lopnr_batches <- split(
  lopnr_vec,
  ceiling(seq_along(lopnr_vec) / batch_size)
)

t_cci_pull <- Sys.time()  # timing: per-batch pull of each cohort member's diagnosis history

# Date-bounded (24 Aug): index dates only span 2015-2023 and CCI looks back 5
# years from each person's index, so no diagnosis before 2010 or after 2023
# is ever needed. Adding this bound server-side avoids pulling each person's
# full lifetime diagnosis history (back to 1997) just to discard most of it
# in the later `filter(datum_date >= index - years(5), datum_date < index)`
# step -- same class of problem as the scb_lisa fix above, smaller scale here.
patients_list <- lapply(
  lopnr_batches,
  function(ids) {

    ids_sql <- paste(ids, collapse = ",")

    diagnoses |>
      filter(
        sql(
          paste0(
            "lopnr IN (",
            ids_sql,
            ") AND indatuma >= toDate('2010-01-01') AND indatuma < toDate('2024-01-01')"
          )
        )
      ) |>
      mutate(
        indatuma_yyyymmdd = dbplyr::sql(
          "formatDateTime(indatuma, '%Y%m%d')"
        )
      ) |>
      select(
        LopNr = lopnr,
        OUTDATE = indatuma_yyyymmdd,
        diagnos
      ) |>
      collect()
  }
)

patients <- bind_rows(patients_list)
cat("CCI diagnosis pull took:", round(difftime(Sys.time(), t_cci_pull, units = "mins"), 2), "min\n")
cat("patients rows pulled:", nrow(patients), "\n")

# ---- Förbered diagnoser för CCI ----

# Byt namn så att CCI-scriptet känner igen kolumnerna
patients <- patients %>%
  rename(
    group = LopNr,
    datum = OUTDATE
  ) %>%
  select(group, datum, diagnos)


# ---- Lägg till indexdatum och filtrera 5 år före index ----

# ---- Lägg till indexdatum och filtrera 5 år före index ----

# Säkerställ en rad per individ
pop_tot_cci <- pop_tot %>%
  select(lopnr, index) %>%
  distinct(lopnr, .keep_all = TRUE)


patients <- patients %>%
  left_join(
    pop_tot_cci,
    by = c("group" = "lopnr")
  ) %>%
  mutate(
    datum_date = as.Date(datum, format = "%Y%m%d")
  ) %>%
  filter(
    datum_date >= index - years(5),
    datum_date < index
  ) %>%
  select(group, datum, diagnos)

# ---- Starta CCI ----

Matrix <- distinct(patients, group)

# Myocardial_infarction
icd7  <- "\\<420,1"
icd8  <- "\\<410|\\<411|\\<412,01|\\<412,91"
icd9  <- "\\<410|\\<412"
icd10 <- "\\<I21|\\<I22|\\<I252"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9  <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10 <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Myocardial_infarction=datum,diagnos.Myocardial_infarction=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Myocardial_infarction=if_else(!is.na(date.Myocardial_infarction),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Myocardial infarction klar\n")

# Congestive_heart_failure
icd7  <- "\\<422,21|\\<422,22|\\<434,1|\\<434,2"
icd8  <- "\\<425,08|\\<425,09|\\<427,0|\\<427,1|\\<428"
icd9  <- paste(c("\\<402A", "402B", "402X", "404A","404B","404X","425E","425F","425H","425W","425X","428"),collapse="|\\<")
icd10 <- "\\<I110|\\<I130|\\<I132|\\<I255|\\<I420|\\<I426|\\<I427|\\<I428|\\<I429|\\<I43|\\<I50"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10 <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Congestive_heart_failure=datum,diagnos.Congestive_heart_failure=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Congestive_heart_failure=if_else(!is.na(date.Congestive_heart_failure),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Congestive heart failure klar\n")

# Peripheral_vascular_disease
icd7  <- "\\<450,1|\\<451|\\<453"
icd8  <- "\\<440|\\<441|\\<443,1|\\<443,9"
icd9  <- "\\<440|\\<441|\\<443B|\\<443X|\\<447B|\\<557"
icd10 <- "\\<I70|\\<I71|\\<I731|\\<I738|\\<I739|\\<I771|\\<I790|\\<I792|\\<K55"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Peripheral_vascular_disease=datum,diagnos.Peripheral_vascular_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Peripheral_vascular_disease=if_else(!is.na(date.Peripheral_vascular_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Peripheral vascular disease klar\n")

# Cerebrovascular_disease
icd7  <- paste(c("\\<330",331:334),collapse="|\\<")
icd8  <- "\\<430|\\<431|\\<432|\\<433|\\<434|\\<435|\\<436|\\<437|\\<438"
icd9  <- "\\<430|\\<431|\\<432|\\<433|\\<434|\\<435|\\<436|\\<437|\\<438"
icd10 <- "\\<G45|\\<I60|\\<I61|\\<I62|\\<I63|\\<I64|\\<I67|\\<I69"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Cerebrovascular_disease=datum,diagnos.Cerebrovascular_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Cerebrovascular_disease=if_else(!is.na(date.Cerebrovascular_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Cerebrovascular_disease klar\n")

# Chronic_obstructive_pulmonary_disease
icd7  <- "\\<502|\\<527,1"
icd8  <- "\\<491|\\<492"
icd9  <- "\\<491|\\<492|\\<496"
icd10 <- "\\<J43|\\<J44"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Chronic_obstructive_pulmonary_disease=datum,diagnos.Chronic_obstructive_pulmonary_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Chronic_obstructive_pulmonary_disease=if_else(!is.na(date.Chronic_obstructive_pulmonary_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Chronic_obstructive_pulmonary_disease klar\n")

# Chronic_other_pulmonary_disease
icd7  <- paste(c("\\<241",501,523:526),collapse="|\\<")
icd8  <- paste(c("\\<490",493,515:518),collapse="|\\<")
icd9  <- paste(c("\\<490",493:495,500:508,516,517),collapse="|\\<")
icd10 <- paste(c("\\<J45",41,42,46,47,60:70),collapse="|\\<J")

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Chronic_other_pulmonary_disease=datum,diagnos.Chronic_other_pulmonary_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Chronic_other_pulmonary_disease=if_else(!is.na(date.Chronic_other_pulmonary_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Chronic_other_pulmonary_disease klar\n")

# Rheumatic_disease
icd7  <- paste(c("\\<722,00","722,01","722,10","722,20","722,23","456,0","456,1","456,2","456,3"),collapse="|\\<")
icd8  <- paste(c("\\<446",696,"712,0","712,1","712,2","712,3","712,5", 716, "734,0", "734,1", "734,9"),collapse="|\\<")
icd9  <- paste(c("\\<446","696A","710A","710B","710C","710D","710E",714,"719D",720,725),collapse="|\\<")
icd10 <- paste(c("\\<M05","06",123,"070","071","072","073","08",13,30,313:316,32:34,350:351,353,45:46),collapse="|\\<M")

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Rheumatic_disease=datum,diagnos.Rheumatic_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Rheumatic_disease=if_else(!is.na(date.Rheumatic_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Rheumatic_disease klar\n")

# Dementia
icd7  <- "\\<304|\\<305"
icd8  <- "\\<290"
icd9  <- "\\<290|\\<294B|\\<331A|\\<331B|\\<331C|\\<331X"
icd10 <- "\\<F00|\\<F01|\\<F02|\\<F03|\\<F051|\\<G30|\\<G311|\\<G319"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Dementia=datum,diagnos.Dementia=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Dementia=if_else(!is.na(date.Dementia),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Dementia klar\n")

# Hemiplegia
icd7 <- "\\<351|\\<352|\\<357,00"
icd8 <- "\\<343|\\<344"
icd9 <- "\\<342|\\<343|\\<344A|\\<344B|\\<344C|\\<344D|\\<344E|\\<344F"
icd10 <- "\\<G114|\\<G80|\\<G81|\\<G82|\\<G830|\\<G831|\\<G832|\\<G833|\\<G838"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Hemiplegia=datum,diagnos.Hemiplegia=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Hemiplegia=if_else(!is.na(date.Hemiplegia),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Hemiplegia klar\n")

# Diabetes_without_chronic_complication
icd7 <- "\\<260,09"
icd8 <-  "\\<250,00|\\<250,07|\\<250,08"
icd9 <- "\\<250A|\\<250B|\\<250C"
icd10 <- "\\<E100|\\<E101|\\<E106|\\<E108|\\<E109|\\<E110|\\<E111|\\<E118|\\<E119|\\<E120|\\<E121|\\<E129|\\<E130|\\<E131|\\<E139|\\<E140|\\<E141|\\<E149"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Diabetes_without_chronic_complication=datum,diagnos.Diabetes_without_chronic_complication=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Diabetes_without_chronic_complication=if_else(!is.na(date.Diabetes_without_chronic_complication),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Diabetes_without_chronic_complication klar\n")

# Diabetes_with_chronic_complication
icd7 <- "\\<260,2|\\<260,21|\\<260,29|\\<260,3|\\<260,4|\\<260,49|\\<260,99"
icd8 <- "\\<250,01|\\<250,02|\\<250,03|\\<250,04|\\<250,05"
icd9 <- "\\<250D|\\<250E|\\<250F|\\<250G"
icd10 <- "\\<E102|\\<E103|\\<E104|\\<E105|\\<E107|\\<E112|\\<E113|\\<E114|\\<E115|\\<E116|\\<E117|\\<E122|\\<E123|\\<E124|\\<E125|\\<E126|\\<E127|\\<E132|\\<E133|\\<E134|\\<E135|\\<E136|\\<E137|\\<E142|\\<E143|\\<E144|\\<E145|\\<E146|\\<E147"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Diabetes_with_chronic_complication=datum,diagnos.Diabetes_with_chronic_complication=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Diabetes_with_chronic_complication=if_else(!is.na(date.Diabetes_with_chronic_complication),1,0,missing=0))
Matrix <- Matrix %>% mutate(Diabetes_without_chronic_complication=if_else(Diabetes_with_chronic_complication==1,0,Diabetes_without_chronic_complication))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Diabetes_with_chronic_complication klar\n")

# Renal_disease
icd7 <- "\\<592|\\<593|\\<792"
icd8 <- "\\<582|\\<583|\\<584|\\<792|\\<593|\\<403,99|\\<404,99|\\<792,99|\\<Y29,01"
icd9 <- "\\<403A|\\<403B|\\<403X|\\<582|\\<583|\\<585|\\<586|\\<588A|\\<V42A|\\<V45B|\\<V56"
icd10 <- "\\<I120|\\<I131|\\<N032|\\<N033|\\<N034|\\<N035|\\<N036|\\<N037|\\<N052|\\<N053|\\<N054|\\<N055|\\<N056|\\<N057|\\<N11|\\<N18|\\<N19|\\<N250|\\<Q611|\\<Q612|\\<Q613|\\<Q614|\\<Z49|\\<Z940|\\<Z992"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Renal_disease=datum,diagnos.Renal_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Renal_disease=if_else(!is.na(date.Renal_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Renal_disease klar\n")

# Mild_liver_disease
icd7  <- "\\<581"
icd8  <- "\\<070|\\<571|\\<573"
icd9 <-  "\\<070|\\<571C|\\<571E|\\<571F|\\<573"
icd10 <- "\\<B15|\\<B16|\\<B17|\\<B18|\\<B19|\\<K703|\\<K709|\\<K73|\\<K746|\\<K754"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Mild_liver_disease=datum,diagnos.Mild_liver_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Mild_liver_disease=if_else(!is.na(date.Mild_liver_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Mild_liver_disease klar\n")

# liver special
icd8  <- "\\<785,3"
icd9 <- "\\<789F"
icd10 <- "\\<R18"

ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.liver_special=datum,diagnos.liver_special=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Liver_special=if_else(!is.na(date.liver_special),1,0,missing=0))
rm(icd8, icd9, icd10, ICD8, ICD9, ICD10, ptnts)
cat("liver special klar\n")

# moderate severe liver disease
icd7 <- "\\<462,1"
icd8 <- "\\<456,0|\\<571,9|\\<573,02"
icd9 <- "\\<456A|\\<456B|\\<456C|\\<572C|\\<572D|\\<572E"
icd10 <-  "\\<I850|\\<I859|\\<I982|\\<I983"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Severe_liver_disease=datum,diagnos.Severe_liver_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Severe_liver_disease=if_else(!is.na(date.Severe_liver_disease),1,0,missing=0))
Matrix <- Matrix %>% mutate(Severe_liver_disease=if_else(Mild_liver_disease==1 & Liver_special==1,1,Severe_liver_disease))
Matrix <- Matrix %>% mutate(Mild_liver_disease=if_else(Severe_liver_disease==1,0,Mild_liver_disease))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("moderate severe liver disease klar\n")

# Peptic_ulcer_disease
icd7  <- "\\<540|\\<541|\\<542"
icd8  <- "\\<531|\\<532|\\<533|\\<534"
icd9 <- "\\<531|\\<532|\\<533|\\<534"
icd10 <-"\\<K25|\\<K26|\\<K27|\\<K28"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Peptic_ulcer_disease=datum,diagnos.Peptic_ulcer_disease=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Peptic_ulcer_disease=if_else(!is.na(date.Peptic_ulcer_disease),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("peptic ulcer disease klar\n")


# Malignancy
icd7   <- paste(paste("\\<",paste(140:190,collapse = "|\\<"),sep=""), paste("|\\<",paste(192:197,collapse = "|\\<"),sep=""), paste("|\\<",paste(200:204,collapse = "|\\<"),sep=""),sep="")
icd8   <- paste(paste("\\<",paste(c(140:172,174),collapse = "|\\<"),sep=""), paste("|\\<",paste(c(180:207,209),collapse = "|\\<"),sep=""),sep="")
icd9   <- paste(paste("\\<",paste(140:172,collapse = "|\\<"),sep=""), paste("|\\<",paste(174:208,collapse = "|\\<"),sep=""),sep="")
icd10  <- paste("\\<C00|\\<C0",paste(1:9,collapse = "|\\<C0",sep=""),paste("|\\<C",paste(c(10:41,43,45:58,60:76,81:86,88:97),collapse = "|\\<C"),sep=""),sep="")

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.malignancy=datum,diagnos.malignancy=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Malignancy=if_else(!is.na(date.malignancy),1,0,missing=0))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("malignancy klar\n")

# Metastatic_cancer
icd7 <- "\\<156,91|\\<198|\\<199"
icd8 <- "\\<196|\\<197|\\<198|\\<199"
icd9 <- "\\<196|\\<197|\\<198|\\<199A|\\<199B"
icd10 <- "\\<C77|\\<C78|\\<C79|\\<C80"

ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Metastatic_solid_tumor=datum,diagnos.Metastatic_solid_tumor=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Metastatic_solid_tumor=if_else(!is.na(date.Metastatic_solid_tumor),1,0,missing=0))
Matrix <- Matrix %>% mutate(Malignancy=if_else(Metastatic_solid_tumor==1,0,Malignancy))
rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)
cat("Metastatic cancer klar\n")

# Aids
icd9  <- "\\<079J|\\<279K"
icd10 <- "\\<B20|\\<B21|\\<B22|\\<B23|\\<B24|\\<F024|\\<O987|\\<R75|\\<Z219|\\<Z717"

ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Aids=datum,diagnos.Aids=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Aids=if_else(!is.na(date.Aids),1,0,missing=0))
rm(icd9, icd10, ICD9, ICD10, ptnts)
cat("Aids klar\n")

# Calculate the unweighted comorbidity index
Matrix$CCIunw <- Matrix$Myocardial_infarction + Matrix$Congestive_heart_failure + Matrix$Peripheral_vascular_disease + 
              Matrix$Cerebrovascular_disease + Matrix$Chronic_obstructive_pulmonary_disease + Matrix$Chronic_other_pulmonary_disease + 
              Matrix$Rheumatic_disease + Matrix$Dementia + Matrix$Hemiplegia + Matrix$Diabetes_without_chronic_complication + 
              Matrix$Diabetes_with_chronic_complication + Matrix$Renal_disease + Matrix$Mild_liver_disease + Matrix$Severe_liver_disease + 
              Matrix$Peptic_ulcer_disease + Matrix$Malignancy + Matrix$Metastatic_solid_tumor + Matrix$Aids 

# Calculate the weighted comorbidity index
Matrix$CCIw <- Matrix$Myocardial_infarction + Matrix$Congestive_heart_failure + Matrix$Peripheral_vascular_disease + 
              Matrix$Cerebrovascular_disease + Matrix$Chronic_obstructive_pulmonary_disease + Matrix$Chronic_other_pulmonary_disease + 
              Matrix$Rheumatic_disease + Matrix$Dementia + 2*Matrix$Hemiplegia + Matrix$Diabetes_without_chronic_complication + 
              2*Matrix$Diabetes_with_chronic_complication + 2*Matrix$Renal_disease + Matrix$Mild_liver_disease + 3*Matrix$Severe_liver_disease + 
              Matrix$Peptic_ulcer_disease + 2*Matrix$Malignancy + 6*Matrix$Metastatic_solid_tumor + 6*Matrix$Aids 

# Delete date and diagnos information in case not needed 
Matrix <- select(Matrix, -contains("."))

# ---- Lägg tillbaka CCI på hela studiepopulationen ----

Matrix <- Matrix %>%
  rename(lopnr = group)

pop_tot_cci <- pop_tot %>%
  left_join(
    Matrix %>% select(lopnr, CCIunw, CCIw),
    by = "lopnr"
  ) %>%
  mutate(
    CCIunw = replace_na(CCIunw, 0),
    CCIw = replace_na(CCIw, 0)
  )

# ---- Sex (protocol §7) ----
pop_tot_cci <- pop_tot_cci |>
  left_join(pop_desc_tbl |> select(lopnr, kon) |> collect(), by = "lopnr")

# ---- Education (protocol §7): SUN, collapsed to 3 levels ----
# FIX (24 Aug): sun2020niva alone gave ~47% missing education, in both groups
# almost identically -- not random missingness, but a coverage gap: per
# cognet.md, sun2020niva is only populated from its introduction year (2019,
# confirmed via income_raw |> group_by(year) |> summarise(pct_na =
# mean(is.na(sun2020niva))*100) -- 100% NA through 2018, 0% from 2019 on), so
# anyone diagnosed before that had zero eligible rows under the "value at or
# before index" rule below, producing a hard NA rather than a real gap.
# sun2000niva is the longer-running classification and covers the earlier
# years sun2020niva can't -- coalesce prefers the newer vintage where
# available, falls back to the older one otherwise.
#
# Snapshot used: each person's most recent SUN level (2020, falling back to
# 2000) at or before their index year (education treated as not decreasing).
education_lookup <- income_raw |>
  mutate(sun_combined = coalesce(sun2020niva, sun2000niva)) |>
  select(lopnr, year, sun_combined) |>
  filter(!is.na(sun_combined)) |>
  left_join(pop_tot_cci |> select(lopnr, index), by = "lopnr") |>
  filter(year <= year(index)) |>
  group_by(lopnr) |>
  slice_max(year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    sun_level_digit = as.integer(substr(as.character(sun_combined), 1, 1)),
    education_3lvl = case_when(
      sun_level_digit %in% 1:2 ~ "low",
      sun_level_digit %in% 3:4 ~ "intermediate",
      sun_level_digit %in% 5:9 ~ "high",
      TRUE ~ NA_character_
    )
  ) |>
  select(lopnr, sun_combined, sun_level_digit, education_3lvl)

pop_tot_cci <- pop_tot_cci |>
  left_join(education_lookup, by = "lopnr")

# ---- Table 1 (protocol §10): unadjusted, case vs. control ----
table1 <- pop_tot_cci |>
  mutate(group_label = if_else(dementia == 1, "EoD", "Control")) |>
  group_by(group_label) |>
  summarise(
    n = n(),
    age_index_mean = mean(age_index, na.rm = TRUE),
    age_index_sd = sd(age_index, na.rm = TRUE),
    pct_female = mean(kon == 2, na.rm = TRUE) * 100,
    cci_unw_mean = mean(CCIunw, na.rm = TRUE),
    cci_w_mean = mean(CCIw, na.rm = TRUE),
    pct_edu_low = mean(education_3lvl == "low", na.rm = TRUE) * 100,
    pct_edu_intermediate = mean(education_3lvl == "intermediate", na.rm = TRUE) * 100,
    pct_edu_high = mean(education_3lvl == "high", na.rm = TRUE) * 100,
    pct_edu_missing = mean(is.na(education_3lvl)) * 100,
    .groups = "drop"
  )

table1

# ---- Region (added 24 Aug): closing a gap left open in the protocol ----
# case_control was originally built to match 2:1 on birth year, age AND
# county (cognet.md) -- age and sex got verified above (both well balanced),
# but county/region was never pulled or checked. scb_rtb is the source
# (annual snapshots, lopnr+year+lan 1968-2024); per cognet.md, county codes
# can change over time, so this uses the same "closest year at or before
# index" snapshot logic already used for education, not "current" county.
scb_rtb <- tbl(con, in_schema("cognet", "scb_rtb"))

t_rtb <- Sys.time()
lan_list <- lapply(lopnr_batches_income, function(ids) {   # reuses the batch list already built for the LISA pull
  ids_sql <- paste(ids, collapse = ",")
  scb_rtb |>
    filter(sql(paste0("lopnr IN (", ids_sql, ")"))) |>
    select(lopnr, year, lan) |>
    collect()
})
lan_raw <- bind_rows(lan_list)
cat("scb_rtb query took:", round(difftime(Sys.time(), t_rtb, units = "mins"), 2), "min\n")

lan_lookup <- lan_raw |>
  filter(!is.na(lan)) |>
  left_join(pop_tot_cci |> select(lopnr, index), by = "lopnr") |>
  filter(year <= year(index)) |>
  group_by(lopnr) |>
  slice_max(year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(lopnr, lan)

pop_tot_cci <- pop_tot_cci |>
  left_join(lan_lookup, by = "lopnr")

cat("Missing county (lan) after lookup:", sum(is.na(pop_tot_cci$lan)), "of", nrow(pop_tot_cci), "\n")

# Pair-level concordance: does each matched control actually share the case's
# county? This directly tests whether the matching worked, rather than just
# comparing marginal distributions between groups (which can look "balanced"
# in aggregate even with badly-matched individual pairs). Uses `controls`,
# which still has the lopnr/lopnr_fall pairing that pop_tot_cci dropped.
pair_lan <- controls |>
  select(lopnr, lopnr_fall) |>
  inner_join(pop_tot_cci |> filter(dementia == 0) |> select(lopnr, lan_control = lan), by = "lopnr") |>
  inner_join(pop_tot_cci |> filter(dementia == 1) |> select(lopnr_fall = lopnr, lan_case = lan), by = "lopnr_fall") |>
  filter(!is.na(lan_case), !is.na(lan_control))

pct_concordant <- mean(pair_lan$lan_case == pair_lan$lan_control) * 100
cat("Pairs with usable county on both sides:", nrow(pair_lan), "\n")
cat("Pct of matched pairs sharing the same county:", round(pct_concordant, 1), "%\n")

# Add region missingness to Table 1, mirroring pct_edu_missing above
table1 <- table1 |>
  left_join(
    pop_tot_cci |>
      mutate(group_label = if_else(dementia == 1, "EoD", "Control")) |>
      group_by(group_label) |>
      summarise(pct_lan_missing = mean(is.na(lan)) * 100, .groups = "drop"),
    by = "group_label"
  )

table1

# ---- RQ1: primary DiD regression (protocol §10) ----
# NOTE: year_to_index/p1 (the RQ2 plot) were built from pop_tot, which never
# got CCI/sex/education attached -- only pop_tot_cci has those. So this
# rebuilds the person-year panel from income_raw + pop_tot_cci, rather than
# reusing year_to_index, to make sure the regression actually has its
# covariates.
model_data <- income_raw |>
  left_join(pop_tot_cci, by = "lopnr") |>
  mutate(
    diagnosis_month = month(index),
    anchor_year = if_else(diagnosis_month <= 6, year(index), year(index) + 1),
    t = year - anchor_year,
    post = if_else(t >= 0, 1, 0),
    # kept as an explicit category rather than dropped, since education
    # missingness after the sun2000niva fallback should be checked in table1
    # first -- if it's small, dropping instead of pooling into "unknown" is
    # also reasonable
    education_3lvl = if_else(is.na(education_3lvl), "unknown", education_3lvl)
  ) |>
  filter(t >= -5, t <= 3)

cat("model_data rows:", nrow(model_data), " | distinct persons:", n_distinct(model_data$lopnr), "\n")

# install.packages("fixest") if not already installed
library(fixest)

model_rq1 <- feols(
  income ~ dementia * post + age_index + factor(kon) + CCIunw + factor(education_3lvl),
  data = model_data,
  cluster = ~lopnr
)

summary(model_rq1)
# The dementia:post coefficient is the DiD estimate: the average within-window
# income loss attributable to EoD status, over and above whatever controls
# also experienced.

# ---- Robustness check: weighted CCI instead of unweighted ----
# Same model, same data, only CCIunw -> CCIw. If dementia:post stays close to
# non-significant with a similar magnitude/sign, the RQ1 result doesn't
# hinge on which comorbidity weighting was used.
model_rq1_cciw <- feols(
  income ~ dementia * post + age_index + factor(kon) + CCIw + factor(education_3lvl),
  data = model_data,
  cluster = ~lopnr
)

summary(model_rq1_cciw)
etable(model_rq1, model_rq1_cciw, headers = c("CCI unweighted", "CCI weighted"))

# ---- Cohort attrition counts, for a flowchart (report and/or presentation) ----
# Reproduces the pop/controls filtering logic step by step, with a distinct-
# person count captured at each stage. All server-side aggregation where
# possible -- no row-level data collected, just counts. Uses new variable
# names throughout, so this doesn't touch pop/controls/model_data etc.

# Step 1: distinct persons with >=1 qualifying dementia diagnosis, any age, any date
step1_dementia_ever <- diagnoses |>
  filter(sql(diag_pred_sql)) |>
  mutate(diag_group = sql(paste0("if(", mci_pred_sql, ", 'MCI', 'Dementia')"))) |>
  filter(diag_group == "Dementia") |>
  summarise(n = sql("uniqExact(lopnr)")) |>
  collect()
step1_dementia_ever

# Step 2: same persons, first-diagnosis index date + age computed
# (aggregation pushed to SQL; only the small one-row-per-person result is collected)
person_index_all <- diagnoses |>
  select(lopnr, indatuma, diagnos) |>
  left_join(pop_desc_tbl |> select(lopnr, fodelsear), by = "lopnr") |>
  filter(sql(diag_pred_sql)) |>
  mutate(
    lopnr = sql("CAST(lopnr AS Int64)"),
    diag_group = sql(paste0("if(", mci_pred_sql, ", 'MCI', 'Dementia')"))
  ) |>
  filter(diag_group == "Dementia") |>
  group_by(lopnr) |>
  summarise(
    index = min(indatuma),
    fodelsear = min(fodelsear)
  ) |>
  ungroup() |>
  collect() |>
  mutate(age_index = year(index) - fodelsear)
cat("Step 2 - same persons, first-diagnosis date/age computed:", n_distinct(person_index_all$lopnr), "\n")

# Step 3: EoD age criterion alone (40-64 at first diagnosis)
step3_age_ok <- person_index_all |> filter(age_index >= 40, age_index < 65)
cat("Step 3 - after EoD age filter (40-64):", nrow(step3_age_ok), "\n")

# Step 4: + diagnosed within study ascertainment period (2015-2023) -- final EoD cohort
step4_final <- step3_age_ok |> filter(year(index) >= 2015, year(index) <= 2023)
cat("Step 4 - after study period filter (2015-2023), FINAL EoD cohort:", nrow(step4_final), "\n")
cat("  (sanity check: should equal nrow(pop) =", nrow(pop), ")\n")

# ---- Controls funnel ----
raw_controls <- case_control |>
  collect() |>
  filter(lopnr_fall %in% step4_final$lopnr) |>
  mutate(cc_indexdatum = as.Date(indexdatum))
cat("Control step A - raw case_control matches, pairs:", nrow(raw_controls),
    "| distinct control persons:", n_distinct(raw_controls$lopnr_kontroll), "\n")

safeguarded_controls <- raw_controls |>
  left_join(step4_final |> select(lopnr_fall = lopnr, study_index = index), by = "lopnr_fall") |>
  filter(cc_indexdatum <= study_index)
cat("Control step B - after survivorship-bias safeguard, pairs:", nrow(safeguarded_controls),
    "| distinct control persons:", n_distinct(safeguarded_controls$lopnr_kontroll), "\n")

capped_controls <- safeguarded_controls |>
  group_by(lopnr_fall) |>
  slice_min(lopnr_kontroll, n = 2, with_ties = FALSE) |>
  ungroup()
cat("Control step C - after capping at 2 per case, FINAL, pairs:", nrow(capped_controls),
    "| distinct control persons:", n_distinct(capped_controls$lopnr_kontroll), "\n")
cat("  (sanity check: should equal nrow(controls) =", nrow(controls), ")\n")

# ---- Cohort selection flowchart (figure, for report/presentation) ----
# Built with ggplot2 (already loaded via tidyverse) -- no new dependencies.
# All counts pulled from the objects computed above, not typed in, so this
# stays correct if any upstream filter changes.

n1 <- as.numeric(step1_dementia_ever$n)
n2 <- nrow(step3_age_ok)
n3 <- nrow(step4_final)                          # == nrow(pop) == 7,039
n4 <- nrow(raw_controls)
n5 <- nrow(safeguarded_controls)
n6 <- nrow(capped_controls)
n7 <- sum(pop_tot_cci$dementia == 0)             # final analytic controls, after case/control overlap resolution == 11,941
n_overlap <- n6 - n7                             # people reclassified as cases

fc <- function(n) format(n, big.mark = ",")

boxes <- tibble::tribble(
  ~y,   ~title,                            ~subtitle,           ~group,
  8,    "Dementia diagnosis, ever",        paste("n =", fc(n1)), "cases",
  6.5,  "Diagnosed at age 40-64",          paste("n =", fc(n2)), "cases",
  5,    "EoD cohort (cases)",              paste("n =", fc(n3)), "cases",
  3.5,  "Raw matched controls",            paste("n =", fc(n4)), "controls",
  2,    "Survivorship safeguard applied",  paste("n =", fc(n5)), "controls",
  0.5,  "Capped at 2 per case",            paste("n =", fc(n6)), "controls",
  -1,   "Matched controls (final)",        paste("n =", fc(n7)), "controls"
) |>
  mutate(x = 0, w = 3.4, h = 0.8)

arrows <- boxes |>
  arrange(desc(y)) |>
  mutate(y_from = y - h / 2, y_to = lead(y) + lead(h) / 2) |>
  filter(!is.na(y_to))

annotations <- tibble::tribble(
  ~y,     ~label,
  7.25,   paste0("Excl: age <40 or ≥65\nn = ", fc(n1 - n2)),
  5.75,   paste0("Excl: outside 2015-2023\nn = ", fc(n2 - n3)),
  4.25,   "Case_control match, ~2:1",
  2.75,   paste0("Excl: unvalidated window\nn = ", fc(n4 - n5)),
  1.25,   paste0("Excl: capped at 2/case\nn = ", fc(n5 - n6)),
  -0.25,  paste0("Reclassified as cases\nn = ", fc(n_overlap))
)

p_flowchart <- ggplot() +
  geom_rect(
    data = boxes,
    aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = group),
    color = "grey30"
  ) +
  geom_text(data = boxes, aes(x = x, y = y + 0.12, label = title), fontface = "bold", size = 3.3) +
  geom_text(data = boxes, aes(x = x, y = y - 0.18, label = subtitle), size = 3) +
  geom_segment(
    data = arrows,
    aes(x = x, xend = x, y = y_from, yend = y_to),
    arrow = arrow(length = unit(0.15, "cm")), color = "grey30"
  ) +
  geom_text(
    data = annotations,
    aes(x = 2.1, y = y, label = label),
    hjust = 0, size = 2.8, lineheight = 0.9
  ) +
  scale_fill_manual(values = c(cases = "#B5D4F4", controls = "#9FE1CB"), guide = "none") +
  coord_cartesian(xlim = c(-2, 4.8), ylim = c(-1.7, 8.7), clip = "off") +
  theme_void()

p_flowchart

# Export for the presentation/report -- adjust path as needed.
# IMPORTANT: judge the figure from this saved file, not from the small Plots
# pane preview -- the pane can crop both edges of a wide plot like this one
# even when the underlying plot is complete. Open the PNG directly afterward.
ggsave("cohort_flowchart.png", p_flowchart, width = 7, height = 10, dpi = 300)

# ---- Pre-format table(s) for the report ----
# kable() normally auto-detects LaTeX/HTML/etc. from the live Quarto render
# context -- that detection isn't available when this runs as a plain script,
# so format is set explicitly to "latex" since the report's YAML is fixed to
# format: pdf. If the report ever switches output format, this needs to
# change (or move the kable() call back into the .qmd).
table1_kable <- kable(table1, digits = 2, format = "latex", booktabs = TRUE)

# ---- Save results for the report (Project Report Foläk 2 WIP_WR.qmd) ----
# The .qmd renders in its own fresh R session -- it can't see table1,
# model_rq1, etc. from this interactive console. Save everything the report
# references here; the .qmd's setup chunk loads this file instead of needing
# a live DB connection (and your password) just to render. Re-run this line
# whenever an upstream number changes, then re-render the report to match.
saveRDS(
  list(
    table1 = table1,
    table1_kable = table1_kable,
    model_rq1 = model_rq1,
    model_rq1_cciw = model_rq1_cciw,
    p1 = p1,
    p_flowchart = p_flowchart,
    n1 = n1, n2 = n2, n3 = n3, n4 = n4, n5 = n5, n6 = n6, n7 = n7, n_overlap = n_overlap
  ),
  "report_results.rds"
)