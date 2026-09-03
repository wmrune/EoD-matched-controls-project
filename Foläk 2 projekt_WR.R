rm(list = ls())

library(tidyverse)
library(dbplyr)
library(readxl)
library(fixest)
library(arrow)
library(knitr)   

# ---- macOS Unicode-normalization workaround----
feather_path <- function(filename) {
  wd_nfc <- iconv(getwd(), from = "UTF-8-MAC", to = "UTF-8")
  file.path(wd_nfc, filename)
}

# ---- DB toggle----
RUN_FROM_DATABASE <- FALSE

if (RUN_FROM_DATABASE) {

library(DBI)
library(RClickhouse)
library(askpass)

# ---- Connection ----
# Insert your own connection chunk here
 
# ---- Inputs ----
diag_dementia <- c("F00", "F01", "F02", "F03", "F051", "G30", "G31")
diag_mci <- c("F067")

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



# Study population
t_pop <- Sys.time()  
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
  filter(diag_group == "Dementia") |>            
  group_by(lopnr) |>
  summarise(
    index = min(indatuma),
    age_index = year(min(indatuma)) - min(fodelsear)
  ) |>
  ungroup() |>
  filter(
    age_index >= 40, age_index < 65,              
    year(index) >= 2015, year(index) <= 2023       
  )
cat("pop query took:", round(difftime(Sys.time(), t_pop, units = "mins"), 2), "min\n")

ggplot(pop |> filter(age_index > 0), aes(age_index)) +
  geom_histogram(bins=30)

#Controls
t_cc <- Sys.time()  
controls <- case_control |>
  collect() |>
  filter(lopnr_fall %in% pop$lopnr) |>
  mutate(cc_indexdatum = as.Date(indexdatum))
cat("case_control query took:", round(difftime(Sys.time(), t_cc, units = "mins"), 2), "min\n")

# Safeguard against immortal-time/survivorship bias
controls <- controls |>
  left_join(pop |> select(lopnr_fall = lopnr, study_index = index), by = "lopnr_fall") |>
  filter(cc_indexdatum <= study_index)

# Cap at 2 controls per case 
controls <- controls |>
  group_by(lopnr_fall) |>
  slice_min(lopnr_kontroll, n = 2, with_ties = FALSE) |>
  ungroup()

#Study population
controls <- controls %>%
  select(lopnr = lopnr_kontroll, lopnr_fall) |>
  left_join(pop_desc_tbl |> select(lopnr, fodelsear) |> collect(), by = "lopnr") |>
  distinct(lopnr, lopnr_fall, fodelsear)
  
#One lopnr per individual
tot_lopnr <- bind_rows(
  pop |> select(lopnr),
  controls |> select(lopnr)
) |>
  distinct(lopnr)

#Calculated annual income before and after index
t_lisa <- Sys.time()
lopnr_vec_income <- tot_lopnr |> pull(lopnr) |> unique()
lisa_batch_size <- 5000
lopnr_batches_income <- split(lopnr_vec_income, ceiling(seq_along(lopnr_vec_income) / lisa_batch_size))

income_list <- lapply(lopnr_batches_income, function(ids) {
  ids_sql <- paste(ids, collapse = ",")
  scb_lisa |>
    filter(sql(paste0("lopnr IN (", ids_sql, ")"))) |>
    filter(year > 2000, year <= 2024) |>              
    select(lopnr, year, raks_forvink, inkkalla1_forv, sun2020niva, sun2000niva) |>  
    collect()
})

income_raw <- bind_rows(income_list) |>
  mutate(
    income = ifelse(year > 2021, inkkalla1_forv, raks_forvink)
  )
cat("scb_lisa query took:", round(difftime(Sys.time(), t_lisa, units = "mins"), 2), "min\n")
cat("income_raw rows pulled:", nrow(income_raw), "\n")

# Outlier handling (Winsorization)
income_cap <- quantile(income_raw$income, 0.99, na.rm = TRUE)
income_raw <- income_raw |>
  mutate(income = pmin(income, income_cap))


# Creating study population
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

#Charlson Comorbidity Index
# Required packages
library(dplyr)

lopnr_vec <- tot_lopnr |>
  pull(lopnr) |>
  unique()

batch_size <- 5000

lopnr_batches <- split(
  lopnr_vec,
  ceiling(seq_along(lopnr_vec) / batch_size)
)

t_cci_pull <- Sys.time()  # timing: per-batch pull of each cohort member's diagnosis history

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

# ---- Preparing diagnoses for CCI ----

patients <- patients %>%
  rename(
    group = LopNr,
    datum = OUTDATE
  ) %>%
  select(group, datum, diagnos)

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

# ---- Start CCI ----

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
cat("Myocardial infarction done\n")

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
cat("Congestive heart failure done\n")

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
cat("Peripheral vascular disease done\n")

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
cat("Cerebrovascular_disease done\n")

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
cat("Chronic_obstructive_pulmonary_disease done\n")

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
cat("Chronic_other_pulmonary_disease done\n")

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
cat("Rheumatic_disease done\n")

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
cat("Dementia done\n")

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
cat("Hemiplegia done\n")

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
cat("Diabetes_without_chronic_complication done\n")

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
cat("Diabetes_with_chronic_complication done\n")

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
cat("Renal_disease done\n")

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
cat("Mild_liver_disease done\n")

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
cat("liver special done\n")

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
cat("moderate severe liver disease done\n")

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
cat("peptic ulcer disease done\n")

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
cat("malignancy done\n")

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
cat("Metastatic cancer done\n")

# Aids
icd9  <- "\\<079J|\\<279K"
icd10 <- "\\<B20|\\<B21|\\<B22|\\<B23|\\<B24|\\<F024|\\<O987|\\<R75|\\<Z219|\\<Z717"

ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
ptnts <- bind_rows(ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Aids=datum,diagnos.Aids=diagnos) 
Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
Matrix <- Matrix %>% mutate(Aids=if_else(!is.na(date.Aids),1,0,missing=0))
rm(icd9, icd10, ICD9, ICD10, ptnts)
cat("AIDS done\n")

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

Matrix <- select(Matrix, -contains("."))

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

library(arrow)
write_feather(income_raw, feather_path("income_raw.feather"))
write_feather(pop_tot_cci, feather_path("pop_tot_cci.feather"))

pop_tot_cci <- pop_tot_cci |>
  left_join(pop_desc_tbl |> select(lopnr, kon) |> collect(), by = "lopnr")

# ---- Education SUN, collapsed to 3 levels ----
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
  
write_feather(income_raw, feather_path("income_raw.feather"))
write_feather(pop_tot, feather_path("pop_tot.feather"))
write_feather(pop_tot_cci, feather_path("pop_tot_cci_full.feather"))
write_feather(controls, feather_path("controls.feather"))

} 

income_raw   <- read_feather(feather_path("income_raw.feather"))
pop_tot      <- read_feather(feather_path("pop_tot.feather"))
pop_tot_cci  <- read_feather(feather_path("pop_tot_cci_full.feather"))
controls     <- read_feather(feather_path("controls.feather"))

# ---- RQ2 plot  ----

income_by_year <- income_raw |> left_join(pop_tot, by = "lopnr")

year_to_index <- income_by_year |>
  mutate(
    diagnosis_month = month(index),
    anchor_year = if_else(diagnosis_month <= 6, year(index), year(index) + 1),
    t = year - anchor_year,      # relative year: t = 0 is the first full post-diagnosis year
    post = if_else(t >= 0, 1, 0)
  ) |>
  filter(t >= -5, t <= 3)         # analysis window per protocol §3

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

# Pair-level concordance
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
model_data <- income_raw |>
  left_join(pop_tot_cci, by = "lopnr") |>
  mutate(
    diagnosis_month = month(index),
    anchor_year = if_else(diagnosis_month <= 6, year(index), year(index) + 1),
    t = year - anchor_year,
    post = if_else(t >= 0, 1, 0),
    
   
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

# ---- Robustness check: weighted CCI instead of unweighted ----
model_rq1_cciw <- feols(
  income ~ dementia * post + age_index + factor(kon) + CCIw + factor(education_3lvl),
  data = model_data,
  cluster = ~lopnr
)

summary(model_rq1_cciw)
etable(model_rq1, model_rq1_cciw, headers = c("CCI unweighted", "CCI weighted"))

# ---- Event-study (pre-trends test) ----
model_data <- model_data |>
  mutate(income_sek = income * 100)

model_event <- feols(
  income_sek ~ i(t, dementia, ref = -1) + factor(t) + dementia + age_index + factor(kon) + CCIw + factor(education_3lvl),
  data = model_data,
  cluster = ~lopnr
)

summary(model_event)

# Event-study coefficient plot
tmp_png <- tempfile(fileext = ".png")
png(tmp_png, width = 7, height = 5, units = "in", res = 300)
iplot(model_event,
      xlab = "Years relative to diagnosis",
      ylab = "Extra EoD-control gap vs. t = -1 (SEK/year)",
      main = "Event-study: EoD vs. control income gap, relative to t = -1")
dev.off()
file.copy(tmp_png, "event_study_plot.png", overwrite = TRUE)

# ---- Cohort attrition counts, for a flowchart (report and/or presentation) ----
if (RUN_FROM_DATABASE) {

# Step 1: distinct persons with >=1 qualifying dementia diagnosis, any age, any date
step1_dementia_ever <- diagnoses |>
  filter(sql(diag_pred_sql)) |>
  mutate(diag_group = sql(paste0("if(", mci_pred_sql, ", 'MCI', 'Dementia')"))) |>
  filter(diag_group == "Dementia") |>
  summarise(n = sql("uniqExact(lopnr)")) |>
  collect()
step1_dementia_ever

# Step 2: same persons, first-diagnosis index date + age computed
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

n1 <- as.numeric(step1_dementia_ever$n)
n2 <- nrow(step3_age_ok)
n3 <- nrow(step4_final)                         
n4 <- nrow(raw_controls)
n5 <- nrow(safeguarded_controls)
n6 <- nrow(capped_controls)
n7 <- sum(pop_tot_cci$dementia == 0)             
n_overlap <- n6 - n7                             

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
ggsave("cohort_flowchart.png", p_flowchart, width = 7, height = 10, dpi = 300)

} 

if (!RUN_FROM_DATABASE) {
  
  prev_results <- readRDS("report_results.rds")
  p_flowchart <- prev_results$p_flowchart
  n1 <- prev_results$n1; n2 <- prev_results$n2; n3 <- prev_results$n3
  n4 <- prev_results$n4; n5 <- prev_results$n5; n6 <- prev_results$n6
  n7 <- prev_results$n7; n_overlap <- prev_results$n_overlap
  rm(prev_results)
}

# ---- Pre-format table(s) for the report ----
table1_kable <- kable(table1, digits = 2, format = "latex", booktabs = TRUE)

# ---- Save results for the report----
saveRDS(
  list(
    table1 = table1,
    table1_kable = table1_kable,
    model_rq1 = model_rq1,
    model_rq1_cciw = model_rq1_cciw,
    model_event = model_event,
    p1 = p1,
    p_flowchart = p_flowchart,
    n1 = n1, n2 = n2, n3 = n3, n4 = n4, n5 = n5, n6 = n6, n7 = n7, n_overlap = n_overlap
  ),
  "report_results.rds"
)