# 0. Identification ---------------------------------------------------

# Title: Validation of data from the Junio Clima Social Survey
# Responsible: Andreas Laffert

# Executive Summary: This script contains the code to data preparation and validation
# Date: June 9, 2026

# 1. Packages  -----------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse,
               here,
               car,
               sjmisc,
               sjlabelled,
               rlang,
               haven,
               codebook)

options(scipen=999)
rm(list = ls())

# 2. Data -----------------------------------------------------------------

db <- haven::read_sav(file = here("input/data/UDP Junio 2026 - Datos.sav"))

glimpse(db)

# 2.1 Codebook ----

codebook_db <- codebook::codebook_table(db) |>
  select(name, label, value_labels)

View(codebook_db)
library(writexl)

#write_xlsx(codebook_db, here("output/codebook_junio.xlsx"))
# 3. Processing -----------------------------------------------------------

names(db)

# 3.1 Quotas ----

# sexo
frq(db$resp_gender)

db |> 
  group_by(resp_gender) |> 
  summarise(n = sum(PONDERADOR)) |> 
  ungroup() |> 
  mutate(p = (n/sum(n))*100)
  
# edad
frq(db$QUOTAGERANGE)

frq(db$QUOTAGERANGE, weights = db$PONDERADOR)

# macrozona
frq(db$CL02STATEclone1)

frq(db$CL02STATEclone1, weights = db$PONDERADOR) 

# nse
db <- db |>
  mutate(gse = case_when(
    CL04NSE <= 3 ~ "ABC1",
    CL04NSE == 4 ~ "C2",
    CL04NSE == 5 ~ "C3",
    CL04NSE == 6 ~ "D",
    CL04NSE == 7 ~ "E"
  ))

frq(db$gse)

frq(db$gse, weights = db$PONDERADOR)

# educacion

frq(db$IDP60clone1) # segun estos datos tecnia superior es el 29% de la muestra. De hecho, segun el ppt, las tres categ suman 96% no 100%

frq(db$IDP60clone1, weights = db$PONDERADOR)

# quintiles

frq(db$IDP62)
frq(db$IDP62clone1)
frq(db$IDP62clone1, weights = db$PONDERADOR)

# ponderador
sjmisc::descr(db$PONDERADOR)

cv <- sd(db$PONDERADOR, na.rm = TRUE) / mean(db$PONDERADOR, na.rm = TRUE)
cv

w <- db$PONDERADOR[!is.na(db$PONDERADOR)]
n <- length(w)
deff <- (n * sum(w^2)) / (sum(w))^2
deff

n_eff <- n / deff
n_eff

# 3.1 Coyuntura ----

# p1

frq(db$IDP22, weights = db$PONDERADOR)

total <- db |>
  filter(!is.na(IDP22)) |>
  group_by(IDP22) |>
  summarise(n = sum(PONDERADOR, na.rm = TRUE)) |>
  ungroup() |>
  mutate(TOTAL = round((n / sum(n)) * 100)) |>
  select(IDP22, TOTAL)

por_genero <- db |>
  filter(!is.na(IDP22), !is.na(resp_gender)) |>
  group_by(IDP22, resp_gender) |>
  summarise(n = sum(PONDERADOR, na.rm = TRUE), .groups = "drop") |>
  group_by(resp_gender) |>
  mutate(p = round((n / sum(n)) * 100)) |>
  ungroup() |>
  select(IDP22, resp_gender, p) |>
  pivot_wider(names_from = resp_gender, values_from = p, names_prefix = "genero_")

por_edad <- db |>
  filter(!is.na(IDP22), !is.na(QUOTAGERANGE)) |>
  group_by(IDP22, QUOTAGERANGE) |>
  summarise(n = sum(PONDERADOR, na.rm = TRUE), .groups = "drop") |>
  group_by(QUOTAGERANGE) |>
  mutate(p = round((n / sum(n)) * 100)) |>
  ungroup() |>
  select(IDP22, QUOTAGERANGE, p) |>
  pivot_wider(names_from = QUOTAGERANGE, values_from = p, names_prefix = "edad_")

por_educ <- db |>
  filter(!is.na(IDP22), !is.na(IDP60clone1)) |>
  group_by(IDP22, IDP60clone1) |>
  summarise(n = sum(PONDERADOR, na.rm = TRUE), .groups = "drop") |>
  group_by(IDP60clone1) |>
  mutate(p = round((n / sum(n)) * 100)) |>
  ungroup() |>
  select(IDP22, IDP60clone1, p) |>
  pivot_wider(names_from = IDP60clone1, values_from = p, names_prefix = "educ_")

db |>
  filter(!is.na(IDP22), !is.na(IDP62clone1)) |>
  group_by(IDP22, IDP62clone1) |>
  summarise(n = sum(PONDERADOR, na.rm = TRUE), .groups = "drop") |>
  group_by(IDP62clone1) |>
  mutate(p = round((n / sum(n)) * 100)) |>
  ungroup() |>
  select(IDP22, IDP62clone1, p) |>
  pivot_wider(names_from = IDP62clone1, values_from = p, names_prefix = "q_")

# Test para todas las emociones con corrección de Bonferroni
emociones <- unique(db$IDP22[!is.na(db$IDP22)])

resultados_test <- map_dfr(emociones, function(val) {
  conteos <- db |>
    filter(!is.na(IDP22), !is.na(resp_gender)) |>
    group_by(resp_gender) |>
    summarise(
      exitos = sum(PONDERADOR[IDP22 == val], na.rm = TRUE),
      total  = sum(PONDERADOR, na.rm = TRUE),
      .groups = "drop"
    )

  test <- prop.test(round(conteos$exitos), round(conteos$total))

  tibble(
    IDP22    = val,
    p_hombre = test$estimate[1],
    p_mujer  = test$estimate[2],
    p_valor  = test$p.value
  )
}) |>
  mutate(p_valor_adj = p.adjust(p_valor, method = "bonferroni")) |>
  arrange(p_valor)

resultados_test


# Test para todas la emergencia con corrección de Bonferroni
emergencia <- unique(db$IDP26[!is.na(db$IDP26)])

resultados_test <- map_dfr(emergencia, function(val) {
  conteos <- db |>
    filter(!is.na(IDP26), !is.na(QUOTAGERANGE)) |>
    group_by(QUOTAGERANGE) |>
    summarise(
      exitos = sum(PONDERADOR[IDP26 == val], na.rm = TRUE),
      total  = sum(PONDERADOR, na.rm = TRUE),
      .groups = "drop"
    )

  test <- prop.test(round(conteos$exitos), round(conteos$total))

  tibble(
    IDP26    = val,
    p_1 = test$estimate[1],
    p_2  = test$estimate[2],
      p_3  = test$estimate[3],
    p_valor  = test$p.value
  )
}) |>
  mutate(p_valor_adj = p.adjust(p_valor, method = "bonferroni")) |>
  arrange(p_valor)

resultados_test

# p kast
frq(db$IDP35, weights = db$PONDERADOR)
frq(db$IDP36, weights = db$PONDERADOR)


db |> 
  filter(IDP35 %in% c(1, 2)) |>
  group_by(IDP35, IDP36) |> 
  summarise(n = sum(PONDERADOR)) |> 
  ungroup() |> 
  group_by(IDP35) |> 
  mutate(p = (n/sum(n))*100)


db |> 
  filter(IDP35 %in% c(1, 2)) |>
  group_by(IDP35, IDP37) |> 
  summarise(n = sum(PONDERADOR)) |> 
  ungroup() |> 
  group_by(IDP35) |> 
  mutate(p = (n/sum(n))*100)

# Experimento ----

frq(db$SCREENER4)

db <- db |> 
  rowwise() |> 
  mutate(redis = mean(c(IDP58L431, IDP58L432), na.rm = TRUE)) |>
  ungroup() 

# Filtrar solo grupos relevantes para cada comparación
db_exp <- db |> filter(SCREENER4 %in% c(1, 2, 3))
db_exp <- db_exp |> 
  mutate(treatment = case_when(
    SCREENER4 == 1 ~ "Treatment 1",
    SCREENER4 == 2 ~ "Treatment 2",
    SCREENER4 == 3 ~ "Control"
  ),
      treatment = factor(treatment, levels = c("Control", "Treatment 1", "Treatment 2")))
  
frq(db_exp$treatment)

# --- 1. Comparación global: ANOVA (los tres grupos) ---
lm(redis ~ treatment, data = db_exp) |> 
  car::Anova(type = "III") 

# --- 2. Comparaciones pareadas (t-test) ---

# Tratamiento 1 vs Control (3)
t.test(redis ~ treatment, data = db_exp |> filter(treatment %in% c("Treatment 1", "Control")))

# Tratamiento 2 vs Control (3)
t.test(redis ~ treatment, data = db_exp |> filter(treatment %in% c("Treatment 2", "Control")))

# --- 3. Medias descriptivas por grupo ---
db_exp |>
  group_by(treatment) |>
  summarise(
    n    = n(),
    mean = mean(redis, na.rm = TRUE),
    sd   = sd(redis,   na.rm = TRUE)
  )


lm(redis ~ treatment + QUOTAGERANGE + resp_gender + IDP60clone1 + IDP62clone1, data = db_exp) |>
  car::Anova(type = "III")


db_exp$IDP58L431_dic <- if_else(db_exp$IDP58L431 >= 5, 1, 0)

lm(redis ~ treatment, data = db_exp) |> 
  summary()


# Cruce Maca ------


db |> 
  group_by(IDP37, IDP36) |>
  summarise(n = sum(PONDERADOR)) |> 
  ungroup() |> 
  group_by(IDP37) |> 
  mutate(p = (n/sum(n))*100)

sjPlot::sjt.xtab(db$IDP37, db$IDP36, 
                 weight.by = db$PONDERADOR,
                 show.col.prc = TRUE)


sig_stars <- function(p) case_when(
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  TRUE      ~ "ns"
)

p37 <- unique(db$IDP37[!is.na(db$IDP37)])

resultados_test <- map_dfr(p37, function(val) {
  conteos <- db |>
    filter(!is.na(IDP37), !is.na(IDP36)) |>
    group_by(IDP36) |>
    summarise(
      exitos = sum(PONDERADOR[IDP37 == val], na.rm = TRUE),
      total  = sum(PONDERADOR, na.rm = TRUE),
      .groups = "drop"
    )

  test <- prop.test(round(conteos$exitos), round(conteos$total))
  pw   <- pairwise.prop.test(round(conteos$exitos), round(conteos$total),
                              p.adjust.method = "bonferroni")

  tibble(
    IDP37      = val,
    p_1        = test$estimate[1],
    p_2        = test$estimate[2],
    p_3        = test$estimate[3],
    p_omnibus  = test$p.value,
    p_1vs2     = pw$p.value[1, 1],
    p_1vs3     = pw$p.value[2, 1],
    p_2vs3     = pw$p.value[2, 2]
  )
}) |>
  mutate(
    sig_1vs2 = sig_stars(p_1vs2),
    sig_1vs3 = sig_stars(p_1vs3),
    sig_2vs3 = sig_stars(p_2vs3)
  ) |>
  arrange(p_omnibus)

resultados_test

