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
               psych,
               poLCA,
               lavaan,
               corrplot)

options(scipen=999)
rm(list = ls())

library(conflicted)

conflict_prefer_all("dplyr")

# 2. Data -----------------------------------------------------------------

db <- haven::read_sav(file = here("input/data/UDP Junio 2026 - Datos.sav"))

glimpse(db)

# 3. Processing -----------------------------------------------------------------

db <- janitor::clean_names(db)

vars_lca <- c("idp38", 
              "idp40l366", "idp40l367", "idp40l368", "idp40l369",
              "idp42l375", "idp42l376", "idp42l377", "idp42l378", "idp42l379", "idp42l380")

db <- db |>
  mutate(gse = case_when(
    cl04nse <= 3 ~ "ABC1",
    cl04nse == 4 ~ "C2",
    cl04nse == 5 ~ "C3",
    cl04nse == 6 ~ "D",
    cl04nse == 7 ~ "E"
  ))

# 4. Analysis -----------------------------------------------------------------

db_proc <- db |> 
  dplyr::select(vars_lca, 
                age_group, 
                resp_gender,
                cl02state,
                idp60,
                idp62,
                idp62clone1,
                idp67,
                idp67clone1,
                idp63,
                idp64,
                idp65,
                idp66,
                idp23,
                idp24,
                idp35,
                idp36,
                idp48l401,
                idp48l402,
                idp50l408,
                idp50l409,
                idp50l410,
                idp50l411,
                idp56l421,
                idp56l422,
                idp56l423,
                idp58l431,
                idp58l432,
                ponderador
              )

# 4.1 Distributions

db_proc |> 
  frq(vars_lca)

db_proc <- db_proc |>
  mutate(
    idp38_r = case_when(
      idp38 %in% c(1, 2) ~ 1,  # Nada o poco
      idp38 == 3 ~ 2,          # Algo
      idp38 == 4 ~ 3,          # Bastante
      idp38 == 5 ~ 4           # Mucho
    ),
    across(
      c(idp40l366, idp40l367, idp40l368, idp40l369),
      ~ case_when(
          .x %in% c(4, 5) ~ 4,  # De acuerdo o muy de acuerdo
          TRUE ~ .x
        ),
      .names = "{.col}_r"
    )
  ) |> 
  mutate(across(starts_with("idp42"), 
                ~ case_when(.x == 1 ~ 1,   # Estado
                            .x == 3 ~ 2,   # Mixto
                            .x == 2 ~ 3),  # Privado
                .names = "{.col}_ord"))

# 4.2 Correlaciones

M <- db_proc |> 
  select(idp38_r, idp40l366_r, idp40l367_r, idp40l368_r, idp40l369_r,
         idp42l375_ord, idp42l376_ord, idp42l377_ord, idp42l378_ord, idp42l379_ord, idp42l380_ord) |> 
  psych::polychoric()

M$rho

diag(M$rho) <- NA

rownames(M$rho) <- c("A. Intervención Estado",
                     "B. Eº grande vs servicios sociales",
                     "C. Crecimiento vs desigualdad",
                     "D. Protección ambiental vs empleo",
                     "E. Dº laborales vs contratación",
                     "F. Provisión: Pensiones",
                     "G. Provisión: Salud",
                     "H. Provisión: Educación escolar",
                     "I. Provisión: Educación universitaria",
                     "J. Provisión: Cuidado niños",
                     "K. Provisión: Cuidado mayores")

#set Column names of the matrix
colnames(M$rho) <-c("(A)", "(B)","(C)","(D)","(E)","(F)","(G)",
                    "(H)", "(I)", "(J)", "(K)")


corrplot::corrplot(M$rho,
                   method = "color",
                   addCoef.col = "black",
                   type = "upper",
                   tl.col = "black",
                   col = colorRampPalette(c("#E16462", "white", "#0D0887"))(12),
                   bg = "white",
                   na.label = "-") 

# 4.3 EFA

db_efa <- db_proc |> 
  select(idp38_r, idp40l366_r, idp40l367_r, idp40l368_r, idp40l369_r,
         idp42l375_ord, idp42l376_ord, idp42l377_ord, idp42l378_ord, idp42l379_ord, idp42l380_ord)

KMO(db_efa)

cortest.bartlett(db_efa)

set.seed(231018) # Resultado reproducible

fa.parallel(db_efa, fa = "fa", 
            fm = "ml", n.iter = 100)

plot(scree(db_efa))

efa_1f <- fa(db_efa, # datos
   nfactors = 1, # n factores
   fm = "ml", # método estimación
   rotate = "oblimin")

efa_2f <- fa(db_efa, # datos
   nfactors = 2, # n factores
   fm = "ml", # método estimación
   rotate = "oblimin")

efa_3f <- fa(db_efa, # datos
   nfactors = 3, # n factores
   fm = "ml", # método estimación
   rotate = "oblimin")

efa_4f <- fa(db_efa, # datos
   nfactors = 4, # n factores
   fm = "ml", # método estimación
   rotate = "oblimin")

fa(db_efa, # datos
   nfactors = 4, # n factores
   fm = "ml", # método estimación
   rotate = "varimax")$Vaccounted

print(efa_1f$loadings, cutoff = 0.30)
print(efa_2f$loadings, cutoff = 0.30)
print(efa_3f$loadings, cutoff = 0.30)
print(efa_4f$loadings, cutoff = 0.30)

# mejor modelo es el de 2 factores

vars_lca_final <- db_proc |> 
  select(idp40l367_r, idp40l368_r, idp40l369_r,
         idp42l375_ord, idp42l376_ord, idp42l377_ord, 
         idp42l378_ord, idp42l379_ord, idp42l380_ord)

KMO(vars_lca_final)

efa_final <- fa(vars_lca_final, nfactors = 2, 
                 rotate = "oblimin", fm = "ml", cor = "poly")

print(efa_final$loadings, cutoff = 0.30)

efa_final$communality

print(efa_final)

# ============================================================
# 4.4 LCA
# ============================================================

# --- 0. Verificación previa de codificación ---------------------------------
# poLCA exige enteros sucesivos partiendo en 1 (no factor, no 0-indexado)

sapply(vars_lca_final, class)
sapply(vars_lca_final, function(x) sort(unique(x)))
sapply(vars_lca_final, function(x) sum(is.na(x)))

vars_lca_final <- vars_lca_final |> sjlabelled::remove_all_labels()

# --- 1. Preparar matriz de indicadores (con id explícito) -------------------

db_proc <- db_proc %>%
  mutate(id_row = row_number())  # id estable, creado una sola vez al inicio

a1 <- vars_lca_final$idp40l367_r
a2 <- vars_lca_final$idp40l368_r
a3 <- vars_lca_final$idp40l369_r
a4 <- vars_lca_final$idp42l375_ord
a5 <- vars_lca_final$idp42l376_ord
a6 <- vars_lca_final$idp42l377_ord
a7 <- vars_lca_final$idp42l378_ord
a8 <- vars_lca_final$idp42l379_ord
a9 <- vars_lca_final$idp42l380_ord

db_lca <- data.frame(
  id_row = db_proc$id_row,
  a1, a2, a3, a4, a5, a6, a7, a8, a9
) %>%
  as_tibble()

# Filtrado explícito de NAs preservando el id (deja el pipeline a prueba
# de cambios futuros, aunque la sección 0 ya confirmó 0 NAs)
n_antes <- nrow(db_lca)
db_lca <- db_lca %>% filter(complete.cases(select(., a1:a9)))
n_despues <- nrow(db_lca)

if (n_antes != n_despues) {
  message(sprintf("Se eliminaron %d filas con NA (%d -> %d)",
                   n_antes - n_despues, n_antes, n_despues))
}

# La fórmula solo referencia a1...a9 — id_row queda como columna extra
# en db_lca sin afectar la estimación, pero disponible para el join final
f <- cbind(a1, a2, a3, a4, a5, a6, a7, a8, a9) ~ 1

# --- 2. Helpers ---------------------------------------------------------------

get_K <- function(fit) {
  if (!is.null(fit$posterior)) return(as.integer(ncol(fit$posterior)))
  if (!is.null(fit$probs) && length(fit$probs) > 0) return(as.integer(nrow(fit$probs[[1]])))
  stop("Cannot determine K from poLCA object (no posterior/probs).")
}

get_npar <- function(fit) {
  if (!is.null(fit$npar)) return(as.numeric(fit$npar))
  K <- get_K(fit)
  Rj <- vapply(fit$probs, function(m) if (is.matrix(m)) ncol(m) else NA_integer_, integer(1))
  if (anyNA(Rj)) stop("Cannot compute npar: missing category counts in fit$probs.")
  as.numeric(sum(K * (Rj - 1)) + (K - 1))
}

get_n_patterns_obs <- function(fit) {
  if (!is.null(fit$predcell)) return(as.integer(length(fit$predcell)))
  if (!is.null(fit$observed) && is.matrix(fit$observed)) return(as.integer(nrow(fit$observed)))
  NA_integer_
}

n_patterns_from_data <- function(data_indicators) {
  as.integer(nrow(dplyr::distinct(data_indicators)))
}

entropy_stats <- function(fit) {
  K <- get_K(fit)
  if (K == 1) {
    return(tibble(entropy = NA_real_, mean_max_p = NA_real_,
                  p10_max_p = NA_real_, avg_diag_min = NA_real_,
                  avg_diag_mean = NA_real_))
  }
  P <- pmax(fit$posterior, 1e-12)
  N <- nrow(P)

  rel_entropy <- 1 - (-sum(P * log(P)) / (N * log(K)))
  class_hat <- max.col(P)
  max_p <- apply(P, 1, max)

  avg_diag <- vapply(1:K, function(k) {
    idx <- which(class_hat == k)
    if (length(idx) == 0) NA_real_ else mean(P[idx, k])
  }, numeric(1))

  tibble(
    entropy = rel_entropy,
    mean_max_p = mean(max_p),
    p10_max_p = unname(quantile(max_p, 0.10)),
    avg_diag_min = suppressWarnings(min(avg_diag, na.rm = TRUE)),
    avg_diag_mean = mean(avg_diag, na.rm = TRUE)
  )
}

min_class_share <- function(fit) {
  K <- get_K(fit)
  if (K == 1) return(tibble(min_class_prior = 1, min_class_map = 1))
  Ppost <- fit$posterior
  N <- nrow(Ppost)
  class_hat <- max.col(Ppost)

  tibble(
    min_class_prior = min(fit$P),
    min_class_map   = min(tabulate(class_hat, nbins = K)) / N
  )
}

compute_bvr <- function(fit, data, vars) {
  n <- nrow(data)
  K <- get_K(fit)
  results <- c()
  pairs <- combn(vars, 2, simplify = FALSE)

  for (p in pairs) {
    v1 <- p[1]; v2 <- p[2]
    obs_tab <- table(data[[v1]], data[[v2]])
    cat1 <- sort(unique(data[[v1]]))
    cat2 <- sort(unique(data[[v2]]))
    exp_tab <- matrix(0, nrow = length(cat1), ncol = length(cat2))

    idx1 <- which(vars == v1)
    idx2 <- which(vars == v2)

    for (k in 1:K) {
      pk <- fit$P[k]
      probs1 <- fit$probs[[idx1]][k, ]
      probs2 <- fit$probs[[idx2]][k, ]
      exp_tab <- exp_tab + n * pk * outer(probs1, probs2)
    }

    chi2 <- sum((obs_tab - exp_tab)^2 / exp_tab)
    df <- (length(cat1) - 1) * (length(cat2) - 1)
    results[paste(v1, v2, sep = " - ")] <- chi2 / df
  }

  sort(results, decreasing = TRUE)
}

# --- 3. LCA de 1 clase + residuos bivariados (BVR) ---------------------------

fit_1 <- poLCA(f, data = db_lca, nclass = 1,
               maxiter = 3000, verbose = FALSE, graphs = FALSE)

bvr_1class <- compute_bvr(fit_1, db_lca, paste0("a", 1:9))
round(bvr_1class, 2)

# --- 4. Estimación K = 1 a 5 --------------------------------------------------

fit_lca <- function(K, data, formula, nrep = 50, maxiter = 3000, seed = 123) {
  set.seed(seed)
  poLCA(formula, data = data, nclass = K,
        nrep = nrep, maxiter = maxiter,
        verbose = FALSE, graphs = FALSE)
}

Ks <- 1:5
fits <- lapply(Ks, fit_lca, data = db_lca, formula = f, nrep = 30)
names(fits) <- paste0("K", Ks)

# --- 5. Tabla de comparación: ajuste + LRT + entropía + tamaño mínimo de clase --

N_used <- nrow(db_lca)
m_patterns <- n_patterns_from_data(db_lca)

fit_tbl <- tibble(
  K      = map_int(fits, get_K),
  N      = N_used,
  logLik = map_dbl(fits, ~ as.numeric(.x$llik)),
  npar   = map_dbl(fits, get_npar),
  AIC    = map_dbl(fits, ~ as.numeric(.x$aic)),
  BIC    = map_dbl(fits, ~ as.numeric(.x$bic)),
  Gsq    = map_dbl(fits, ~ as.numeric(.x$Gsq)),
  Chisq  = map_dbl(fits, ~ as.numeric(.x$Chisq))
) %>%
  bind_cols(map_dfr(fits, entropy_stats)) %>%
  bind_cols(map_dfr(fits, min_class_share)) %>%
  arrange(K) %>%
  mutate(
    dAIC = AIC - lag(AIC),
    dBIC = BIC - lag(BIC),
    # LRT heurístico K vs K-1: en LCA los modelos no son estrictamente
    # anidados en el sentido regular (problemas de borde en los parámetros),
    # así que este test es orientativo, no una prueba formal válida.
    LRT    = 2 * (logLik - lag(logLik)),
    df_LRT = npar - lag(npar),
    p_LRT  = ifelse(!is.na(LRT) & df_LRT > 0,
                     pchisq(LRT, df = df_LRT, lower.tail = FALSE), NA_real_)
  )

fit_tbl_report <- fit_tbl %>%
  transmute(
    K, N,
    logLik = round(logLik, 1),
    npar,
    AIC = round(AIC, 1), dAIC = round(dAIC, 1),
    BIC = round(BIC, 1), dBIC = round(dBIC, 1),
    LRT = round(LRT, 1), df_LRT, p_LRT = signif(p_LRT, 3),
    entropy = round(entropy, 3),
    mean_max_p = round(mean_max_p, 3),
    min_class_map = round(min_class_map, 3)
  )

fit_tbl_report

# --- 6. Chequeo de dispersión de patrones (sparseness) -----------------------
# Con 9 ítems politómicos, el espacio de patrones posibles es grande frente
# a N=1500; esto debilita Chisq/Gsq como criterios de ajuste absoluto
# (BIC/AIC siguen siendo válidos para comparar modelos).

n_total_possible <- prod(sapply(select(db_lca, a1:a9), function(x) length(unique(x))))

tibble(
  n_obs = N_used,
  n_patterns_observados = m_patterns,
  n_patterns_posibles = n_total_possible,
  pct_cobertura = round(100 * m_patterns / n_total_possible, 1)
)

# --- 7. Visualizaciones -------------------------------------------------------

fit_long <- fit_tbl %>%
  dplyr::select(K, AIC, BIC) %>%
  pivot_longer(-K, names_to = "criterion", values_to = "value")

ggplot(fit_long, aes(x = K, y = value, group = criterion)) +
  geom_line() +
  geom_point() +
  facet_wrap(~criterion, scales = "free_y") +
  labs(x = "Number of classes (K)", y = "Value")

delta_long <- fit_tbl %>%
  dplyr::select(K, dAIC, dBIC) %>%
  pivot_longer(-K, names_to = "delta", values_to = "value")

ggplot(delta_long, aes(x = K, y = value, group = delta)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point() +
  facet_wrap(~delta, scales = "free_y") +
  labs(x = "K", y = "Change vs K-1 (negative = better)")

ggplot(fit_tbl, aes(x = K, y = entropy)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
  geom_line() +
  geom_point() +
  labs(x = "Number of classes (K)", y = "Entropy (relative)",
       title = "Entropía por número de clases")

fit_tbl_report
bvr_1class

# --- 8. Modelo final K=3: asignación de clase y unión con db_proc -----------

fit3 <- fits[["K3"]]

class_assignment <- tibble(
  id_row      = db_lca$id_row,
  class_modal = max.col(fit3$posterior),
  post_max    = apply(fit3$posterior, 1, max)
) %>%
  mutate(class_modal = factor(class_modal,
                               levels = 1:3,
                               labels = c("Estatista", "Promercado", "Pragmático/mixto")))

db_proc <- db_proc %>%
  left_join(class_assignment, by = "id_row")

# Verificación: cuántos casos de db_proc no recibieron clase (deberían ser
# los mismos que se cayeron por NA en el LCA, si los hubo)
sum(is.na(db_proc$class_modal))
summary(db_proc$post_max)
table(db_proc$class_modal)

# 3 clases son el mejor resultado

# ============================================================
# 4.5 Three latent class model
# ============================================================

fit3 <- fits[["K3"]]

# 1) Etiquetas legibles por ítem
item_labels <- c(
  a1 = "idp40l367 (crecimiento vs desigualdad)",
  a2 = "idp40l368 (empleo vs ambiente)",
  a3 = "idp40l369 (empleo vs derechos laborales)",
  a4 = "idp42l375 (pensiones)",
  a5 = "idp42l376 (salud)",
  a6 = "idp42l377 (educación escolar)",
  a7 = "idp42l378 (educación universitaria)",
  a8 = "idp42l379 (cuidado infantil)",
  a9 = "idp42l380 (cuidado personas mayores)"
)

# 2) Extraer probabilidades condicionales por ítem y clase — versión robusta
#    a tamaños distintos de categoría (evita el bug de NA por bind_rows)

probs_long <- purrr::map_dfr(names(fit3$probs), function(item) {
  mat <- fit3$probs[[item]]
  Rj  <- ncol(mat)
  K   <- nrow(mat)
  
  tibble(
    item    = item,
    class   = rep(1:K, times = Rj),
    cat_num = rep(1:Rj, each = K),
    prob    = as.vector(mat)
  )
})

# 3) Calcular valor esperado normalizado por ítem-clase
p_expected <- probs_long %>%
  group_by(item, class) %>%
  summarise(
    Rj      = n(),
    e_score = sum(cat_num * prob),
    e_norm  = (e_score - 1) / (Rj - 1),
    .groups = "drop"
  )

# Verificación: no debería haber NA ahora
sum(is.na(p_expected$e_norm))  # debería dar 0
# Verificación rápida de nombres de columnas (ajustar regex si el formato difiere)
# colnames(fit3$probs[["a1"]])

# 3) Calcular valor esperado normalizado por ítem-clase
#    (sustituye el P(endorsement) del código dicotómico original, que no
#    aplica acá porque los ítems tienen distinto n° de categorías)
p_expected <- probs_long %>%
  group_by(item, class) %>%
  summarise(
    Rj      = n(),
    e_score = sum(cat_num * prob),        # valor esperado en escala original
    e_norm  = (e_score - 1) / (Rj - 1),   # normalizado 0-1
    .groups = "drop"
  )

# 4) Tabla limpia con nombres de variable legibles
p_endorse <- p_expected %>%
  mutate(
    variable = recode(item, !!!item_labels),
    variable = factor(variable, levels = item_labels[c("a1","a2","a3",
                                                         "a4","a5","a6",
                                                         "a7","a8","a9")])
  )

# 5) Gráfico de perfiles (una línea por clase)
ggplot(p_endorse, aes(x = variable, y = e_norm, group = factor(class), colour = factor(class))) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Ítem", 
       y = "Puntaje esperado normalizado\n(0 = Estado/desacuerdo, 1 = Mercado/acuerdo)",
       colour = "Clase") +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))

# 6) Heatmap
ggplot(p_endorse, aes(x = variable, y = factor(class), fill = e_norm)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#A8384A", limits = c(0, 1)) +
  labs(x = "Ítem", y = "Clase", fill = "Puntaje\nnormalizado") +
  theme(legend.position = "right",
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))

# 7) Tabla de reporte (wide): ítems x clases, con prevalencia de clase en el header
class_prev <- tibble(
  class = 1:length(fit3$P),
  pi    = as.numeric(fit3$P)
) %>%
  mutate(pi_pct = 100 * pi)

report_probs <- p_endorse %>%
  dplyr::select(variable, class, e_norm) %>%
  mutate(class = paste0("Class_", class)) %>%
  pivot_wider(names_from = class, values_from = e_norm) %>%
  mutate(across(starts_with("Class_"), ~ round(.x, 3)))

prev_lab <- class_prev %>%
  transmute(class = paste0("Class_", class),
            lab = paste0(class, " (", round(pi_pct, 1), "%)"))

report_probs2 <- report_probs
for (i in seq_len(nrow(prev_lab))) {
  old <- prev_lab$class[i]
  new <- prev_lab$lab[i]
  names(report_probs2)[names(report_probs2) == old] <- new
}

report_probs2 %>%
  kableExtra::kable(format = "html",
                    align = "c",
                    booktabs = TRUE,
                    escape = FALSE,
                    caption = "Puntaje esperado normalizado por ítem y clase (K=3)") %>%
  kableExtra::kable_styling(full_width = TRUE,
                            latex_options = "hold_position",
                            bootstrap_options = c("striped", "bordered", "condensed"),
                            font_size = 23)

fit3$probs


# ============================================================
# 4.6 Perfil de clases LCA — Variables sociodemográficas e ideología
# ============================================================

glimpse(db_proc)

frq(db_proc$idp62clone1)
frq(db_proc$idp60)
frq(db_proc$resp_gender)
frq(db_proc$age_group)
frq(db_proc$idp67clone1)


db_proc <- db_proc %>%
  mutate(
    
    # --- Educación: colapso a 4 niveles CINE (lógica de nivel más alto completado) ---
    educacion_cine = case_when(
      idp60 %in% c(1, 2, 3, 4, 5) ~ "Básica hasta media",
      idp60 %in% c(6, 7)       ~ "Técnica superior",
      idp60 %in% c(8, 9, 10)      ~ "Universitaria o más"
    ),
    educacion_cine = factor(educacion_cine,
                             levels = c("Básica hasta media",
                                        "Técnica superior", "Universitaria o más"),
                             ordered = TRUE),
    
    # --- Ingreso: quintiles ya vienen codificados 12-16, solo se etiquetan ---
    quintil_ingreso = case_when(
      idp62clone1 == 12 ~ "Q1",
      idp62clone1 == 13 ~ "Q2",
      idp62clone1 == 14 ~ "Q3",
      idp62clone1 == 15 ~ "Q4",
      idp62clone1 == 16 ~ "Q5"
    ),
    quintil_ingreso = factor(quintil_ingreso,
                              levels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                              ordered = TRUE),
    
    # --- Sexo: factor nominal simple ---
    sexo = case_when(
      resp_gender == 1 ~ "Hombre",
      resp_gender == 2 ~ "Mujer"
    ),
    sexo = factor(sexo, levels = c("Hombre", "Mujer")),
    
    # --- Tramo etario: ya viene en tramos, solo se etiqueta y ordena ---
    tramo_edad = case_when(
      age_group == 1 ~ "18-29",
      age_group == 2 ~ "30-49",
      age_group == 3 ~ "50-99"
    ),
    tramo_edad = factor(tramo_edad,
                         levels = c("18-29", "30-49", "50-99"),
                         ordered = TRUE),
    
    # --- Ideología política: nominal, sin orden impuesto entre Izq/Centro/Der/Ninguno ---
    ideologia = case_when(
      idp67clone1 == 13 ~ "Izquierda",
      idp67clone1 == 6  ~ "Centro",
      idp67clone1 == 14 ~ "Derecha",
      idp67clone1 == 12 ~ "Ninguno"
    ),
    ideologia = factor(ideologia,
                        levels = c("Izquierda", "Centro", "Derecha", "Ninguno"))
  )


items_redistrib <- db_proc %>% 
  select(idp56l421, idp56l422, idp56l423, idp58l431, idp58l432)

psych::alpha(items_redistrib)  # revisar raw_alpha antes de seguir

db_proc <- db_proc %>%
  mutate(
    pref_redistributiva = rowMeans(
      select(., idp56l421, idp56l422, idp56l423, idp58l431, idp58l432),
      na.rm = TRUE
    )
  )

summary(db_proc$pref_redistributiva)

db_proc <- db_proc %>%
  mutate(
    eval_gobierno = case_when(
      idp36 == 1 ~ "Positiva",
      idp36 == 2 ~ "Negativa",
      idp36 == 3 ~ "Ni positiva ni negativa"
    ),
    eval_gobierno = factor(eval_gobierno,
                            levels = c("Negativa", "Ni positiva ni negativa", "Positiva"))
  )

db_proc$idp48l401 <- as.numeric(db_proc$idp48l401)
db_proc$idp48l402 <- as.numeric(db_proc$idp48l402)

library(survey)

# Diseño muestral con ponderador de calibración (ajustar nombre si es distinto)
design <- svydesign(ids = ~1, weights = ~ponderador, data = db_proc)

# --- Función genérica: crosstab ponderado + test Rao-Scott -------------------

# Categórica: crosstab ponderado + test Rao-Scott (ya definida antes, se reusa)
tab_categorica <- function(var, design, var_label = var) {
  form <- reformulate(c("class_modal", var))
  tab_n   <- svytable(form, design)
  tab_pct <- round(100 * prop.table(tab_n, margin = 1), 1)
  test    <- svychisq(form, design, statistic = "F")
  
  cat("\n===", var_label, "===\n")
  print(tab_pct)
  cat("\nTest Rao-Scott (F):\n")
  print(test)
  
  invisible(list(tabla_pct = tab_pct, tabla_n = tab_n, test = test))
}

# Continua/Likert: medias ponderadas por clase + test F omnibus
tab_continua <- function(var, design, var_label = var) {
  medias <- svyby(reformulate(var), ~class_modal, design, svymean, na.rm = TRUE)
  modelo <- svyglm(reformulate("class_modal", response = var), design)
  test   <- regTermTest(modelo, ~class_modal)
  
  medias_redondeadas <- medias %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  cat("\n===", var_label, "===\n")
  print(medias_redondeadas)
  cat("\nTest F (omnibus, diferencia entre clases):\n")
  print(test)
  
  invisible(list(medias = medias, test = test))
}

# --- Aplicar a las cinco variables -------------------------------------------

vars_perfil <- c(
  ingreso    = "quintil_ingreso",
  educacion  = "educacion_cine",
  sexo       = "sexo",
  edad       = "tramo_edad",
  ideologia  = "ideologia"
)

resultados_perfil <- map2(vars_perfil, names(vars_perfil), 
                           ~ tab_categorica(.x, design, var_label = .y))


# Categórica
res_eval_gobierno <- tab_categorica("eval_gobierno", design, "Evaluación gobierno actual")

# Continuas — percepción de injusticia
res_idp48l401 <- tab_continua("idp48l401", design, "Injusticia: pago en educación")
res_idp48l402 <- tab_continua("idp48l402", design, "Injusticia: pago en salud")

# Continuas — preferencias fiscales
res_idp50l408 <- tab_continua("idp50l408", design, "Prefiero pagar menos impuestos")
res_idp50l409 <- tab_continua("idp50l409", design, "Justo que ricos paguen más impuestos")
res_idp50l410 <- tab_continua("idp50l410", design, "Problema es mal gasto, no bajos impuestos")
res_idp50l411 <- tab_continua("idp50l411", design, "Bajar impuestos a empresas genera empleo")

# Continua — índice de preferencias redistributivas
res_pref_redistributiva <- tab_continua("pref_redistributiva", design, "Índice preferencias redistributivas")

# Continuas — situación económica del hogar
res_idp63 <- tab_continua("idp63", design, "Cómo llegan a fin de mes")
res_idp64 <- tab_continua("idp64", design, "Alza de precios afectó al hogar")
res_idp65 <- tab_continua("idp65", design, "Situación económica vs. hace 6 meses")
res_idp66 <- tab_continua("idp66", design, "Expectativa económica próximos 6 meses")

