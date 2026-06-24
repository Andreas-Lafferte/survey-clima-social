#******************************************************************************************************************************************************
# 0. Identification -------------------------------------------------------
# Title: Code for generate the parameters for the quota design
# Responsible: Andreas Laffert
# Executive Summary: This script contains the code to perform CENSO and CASEN processing and to generate the parameters for the quota design.
# Date: May 14, 2026
#******************************************************************************************************************************************************

options(scipen=999)
rm(list = ls())

 # 1. Load packages -------------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse,
               here,
               magrittr,
               rio,
               sjmisc,
               sjlabelled,
               rstatix,
               psych,
               srvyr,
               Hmisc,
               writexl)

options(scipen=999)

# 2. Load data -------------------------------------------------------

censo <- rio::import(here("input/data/personas_censo2024.csv")) |> as_tibble()
casen <- rio::import(here("input/data/casen_2024.RData")) |> 
  sjlabelled::remove_all_labels() |> 
  as_tibble()

# 3. CENSO processing -------------------------------------------------------

names(censo)

# select

censo_proc <- censo |> 
  dplyr::select(1:3, region, sexo, edad, edad_quinquenal, cine11) |> 

# filter
  
  filter(edad >= 18) |> 

# recode and transform
  
  mutate(edad_tramos = case_when(edad >= 18 & edad <= 29 ~ "18-29",
                                edad >= 30 & edad <= 44 ~ "30-44",
                                edad >= 45 & edad <= 59 ~ "45-59",
                                edad >= 60 ~ "60+",
                                TRUE ~ NA_character_),
         edad_tramos = factor(edad_tramos,
                              levels = c("18-29", "30-44", "45-59", "60+")),
         sexo_f = if_else(sexo == 1, "Hombre", "Mujer"),
         sexo_f = factor(sexo_f, levels = c("Hombre", "Mujer")),
         macrozona = case_when(
                            region %in% c(15,1,2,3,4) ~ "Norte",
                            region %in% c(5,6,7) ~ "Centro",
                            region == 13 ~ "Metropolitana",
                            region %in% c(16,8,9,14,10,11,12) ~ "Sur y Austral",
                            TRUE ~ NA_character_),
         macrozona = factor(macrozona, levels = c("Norte", 
                                                  "Centro", 
                                                  "Metropolitana", 
                                                  "Sur y Austral")),
         nivel_educ = case_when(
                            cine11 %in% c(1,2,3,4,5) ~ "Básica o menos",      # 01, 02, 03, 10, 14, 98
                            cine11 %in% c(6,7)               ~ "Media",        # 24, 25
                            cine11 == 8                        ~ "Técnica superior", # 35
                            cine11 %in% c(9,10,11)          ~ "Universitaria o más", # 46, 56, 64
                            TRUE                             ~ NA_character_)
    )                                       

censo_proc <- censo_proc |> na.omit()

# Quotas
frq(censo_proc$macrozona)
frq(censo_proc$sexo_f)
frq(censo_proc$edad_tramos)
frq(censo_proc$nivel_educ)

# 4. Export frequencies -------------------------------------------------------

freq_table <- function(df, var) {
  total_n <- nrow(df)
  valid_n <- sum(!is.na(df[[deparse(substitute(var))]]))

  tbl <- df |>
    count({{ var }}, .drop = FALSE, name = "n") |>
    rename(categoria = 1) |>
    mutate(
      categoria      = as.character(categoria),
      porc           = round(n / total_n * 100, 2),
      porc_valido    = round(n / valid_n  * 100, 2),
      porc_acumulado = round(cumsum(porc_valido), 2)
    )

  bind_rows(
    tbl,
    tibble(categoria = "Total", n = sum(tbl$n),
           porc = sum(tbl$porc), porc_valido = sum(tbl$porc_valido),
           porc_acumulado = NA_real_)
  )
}

freq_list <- list(
  macrozona   = freq_table(censo_proc, macrozona),
  sexo        = freq_table(censo_proc, sexo_f),
  edad_tramos = freq_table(censo_proc, edad_tramos),
  nivel_educ  = freq_table(censo_proc, nivel_educ)
)

# 4. CASEN processing -------------------------------------------------------

names(casen)

db_proc <- casen %>%

  # select ---

  select(1:3, hogar, nucleo, varunit, varstrat, expr, nse,
         pco1, o15, ytrabajocorh, ypchautcor, ypchtrabcor, ypchtotcor, ypchmonecor, ymonecorh,
         yaut, yauth, yautcor, yautcorh, ytot, ytoth, ytotcor, ytotcorh,
         dau, qaut) %>%

  # filter ---

  filter(nucleo != 0) # exclude SDPA

# Puntos de corte ponderados para ymonecorh (jefes de hogar) ---

hh_idx <- db_proc$pco1 == 1 & !is.na(db_proc$ymonecorh)

cortes_decil <- Hmisc::wtd.quantile(db_proc$ymonecorh[hh_idx], weights = db_proc$expr[hh_idx],
                                     probs = seq(0.1, 0.9, by = 0.1), normwt = TRUE)

cortes_quintil <- Hmisc::wtd.quantile(db_proc$ymonecorh[hh_idx], weights = db_proc$expr[hh_idx],
                                       probs = seq(0.2, 0.8, by = 0.2), normwt = TRUE)

rm(hh_idx)

db_proc <- db_proc %>%
  mutate(
    decil_ymon   = as.integer(cut(ymonecorh, breaks = c(-Inf, cortes_decil, Inf),
                                  labels = 1:10, include.lowest = TRUE)),
    quintil_ymon = as.integer(cut(ymonecorh, breaks = c(-Inf, cortes_quintil, Inf),
                                  labels = 1:5,  include.lowest = TRUE))
  )

db_pond <- db_proc %>%
  as_survey_design(ids = varunit, strata = varstrat, weights = expr)

# 5. Calculations -------------------------------------------------------

db_pond %>%
  filter(pco1 == 1, !is.na(dau)) %>%
  group_by(dau) %>%
  summarise(ing_mone = survey_mean(ymonecorh, na.rm = TRUE, vartype = "ci"))

db_pond %>%
  filter(pco1 == 1, !is.na(qaut)) %>%
  group_by(qaut) %>%
  summarise(ing_mone = survey_mean(ymonecorh, na.rm = TRUE, vartype = "ci"))

db_pond %>%
  filter(pco1 == 1, !is.na(qaut)) %>%
  group_by(qaut) %>%
  summarise(ing_tot = survey_mean(ytotcorh, na.rm = TRUE, vartype = "ci"))

# 6. Export frequencies -------------------------------------------------------

decil_freq <- db_proc |>
  filter(pco1 == 1, !is.na(dau)) |>
  group_by(dau) |>
  summarise(n = round(sum(expr, na.rm = TRUE))) |>
  mutate(
    categoria      = as.character(dau),
    porc           = round(n / sum(n) * 100, 2),
    porc_valido    = round(n / sum(n) * 100, 2),
    porc_acumulado = round(cumsum(porc_valido), 2)
  ) |>
  select(categoria, n, porc, porc_valido, porc_acumulado)

decil_freq <- bind_rows(
  decil_freq,
  tibble(categoria = "Total", n = sum(decil_freq$n),
         porc = sum(decil_freq$porc), porc_valido = sum(decil_freq$porc_valido),
         porc_acumulado = NA_real_)
)

quintil_freq <- db_proc |>
  filter(pco1 == 1, !is.na(qaut)) |>
  group_by(qaut) |>
  summarise(n = round(sum(expr, na.rm = TRUE))) |>
  mutate(
    categoria      = as.character(qaut),
    porc           = round(n / sum(n) * 100, 2),
    porc_valido    = round(n / sum(n) * 100, 2),
    porc_acumulado = round(cumsum(porc_valido), 2)
  ) |>
  select(categoria, n, porc, porc_valido, porc_acumulado)

quintil_freq <- bind_rows(
  quintil_freq,
  tibble(categoria = "Total", n = sum(quintil_freq$n),
         porc = sum(quintil_freq$porc), porc_valido = sum(quintil_freq$porc_valido),
         porc_acumulado = NA_real_)
)

decil_ymon_freq <- db_proc |>
  filter(pco1 == 1, !is.na(decil_ymon)) |>
  group_by(decil_ymon) |>
  summarise(
    n           = round(sum(expr, na.rm = TRUE)),
    min_ingreso = min(ymonecorh, na.rm = TRUE),
    max_ingreso = max(ymonecorh, na.rm = TRUE)
  ) |>
  mutate(
    categoria      = as.character(decil_ymon),
    porc           = round(n / sum(n) * 100, 2),
    porc_valido    = round(n / sum(n) * 100, 2),
    porc_acumulado = round(cumsum(porc_valido), 2)
  ) |>
  select(categoria, n, porc, porc_valido, porc_acumulado, min_ingreso, max_ingreso)

decil_ymon_freq <- bind_rows(
  decil_ymon_freq,
  tibble(categoria = "Total", n = sum(decil_ymon_freq$n),
         porc = sum(decil_ymon_freq$porc), porc_valido = sum(decil_ymon_freq$porc_valido),
         porc_acumulado = NA_real_, min_ingreso = NA_real_, max_ingreso = NA_real_)
)

quintil_ymon_freq <- db_proc |>
  filter(pco1 == 1, !is.na(quintil_ymon)) |>
  group_by(quintil_ymon) |>
  summarise(
    n           = round(sum(expr, na.rm = TRUE)),
    min_ingreso = min(ymonecorh, na.rm = TRUE),
    max_ingreso = max(ymonecorh, na.rm = TRUE)
  ) |>
  mutate(
    categoria      = as.character(quintil_ymon),
    porc           = round(n / sum(n) * 100, 2),
    porc_valido    = round(n / sum(n) * 100, 2),
    porc_acumulado = round(cumsum(porc_valido), 2)
  ) |>
  select(categoria, n, porc, porc_valido, porc_acumulado, min_ingreso, max_ingreso)

quintil_ymon_freq <- bind_rows(
  quintil_ymon_freq,
  tibble(categoria = "Total", n = sum(quintil_ymon_freq$n),
         porc = sum(quintil_ymon_freq$porc), porc_valido = sum(quintil_ymon_freq$porc_valido),
         porc_acumulado = NA_real_, min_ingreso = NA_real_, max_ingreso = NA_real_)
)

freq_list[["decil_autonomo"]]    <- decil_freq
freq_list[["quintil_autonomo"]]  <- quintil_freq
freq_list[["decil_ymonecorh"]]   <- decil_ymon_freq
freq_list[["quintil_ymonecorh"]] <- quintil_ymon_freq

writexl::write_xlsx(freq_list, here("output/frecuencias_cuotas.xlsx"))



# 7. Brecha salarial etnica -------------------------------------------------------

db_etnia <-  casen %>%

  # select ---

  select(1:3, hogar, nucleo, varunit, varstrat, expr, r3,
         pco1, o15, ytrabajocorh, ypchautcor, ypchtrabcor, ypchtotcor, ypchmonecor, ymonecorh,
         yaut, yauth, yautcor, yautcorh, ytot, ytoth, ytotcor, ytotcorh,
         dau, qaut) %>%

  # filter ---

  filter(nucleo != 0) # exclude SDPA


frq(db_etnia$r3)

db_etnia <- db_etnia |> 
  mutate(etnia = if_else(r3 == 12, "No indigena", "Indigena"))

frq(db_etnia$etnia)

db_etnia |> 
  as_survey_design(ids = varunit, strata = varstrat, weights = expr) |> 
  group_by(etnia) |> 
  summarise(salario = survey_mean(yaut, na.rm = TRUE, vartype = "ci"))
