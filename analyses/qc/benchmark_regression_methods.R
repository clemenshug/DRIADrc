# library(DRIAD)
library(tidyverse)
library(qs)
library(here)
library(microbenchmark)
library(ordinalForest)
library(ordinalNet)
library(OrdinalLogisticBiplot)
library(ordinalgmifs)
library(glmnetcr)
library(ordinalRidge)
library(furrr)

plan(multisession(workers = 6L))

synapser::synLogin()
syn <- synExtra::synDownloader(here("data"))

traindata <- syn("syn25606985") %>%
  qread()

traindata_all <- syn("syn25607441") %>%
  qread()

train_oforest <- function(data) {
  mdl <- ordinalForest::ordfor(depvar = "Label", data = data)
  preds <- predict(mdl, newdata = dplyr::select(data, -Label))
  as.integer(preds$ypred)
}

train_onet <- function(data) {
  Y <- data[["Label"]]
  XX <- data %>% dplyr::select(-Label) %>% t() %>% cov()

  mdl <- ordinalNet::ordinalNet(
    XX, Y, alpha = 0, threshIn = 1e-4, threshOut = 1e-4
  )
  predict(mdl, type = "response")[, "P[Y=7]"]
}

train_pordlogist <- function(data) {
  Y <- as.integer(data[["Label"]])
  XX <- data %>% dplyr::select(-Label) %>% as.matrix()

  preds <- OrdinalLogisticBiplot::pordlogist(
    Y, XX, penalization = 0.1, tol = 1e-04, maxiter = 200, show = FALSE
  )
  preds$fitted.values[, 7]
}

train_ordinalgmifs <- function(data) {
  mdl <- ordinalgmifs::ordinalgmifs(
    Label~1, x = setdiff(names(data), "Label"), data = data, tol = 1e-04
  )
  preds <- predict(mdl)
  preds$predicted
}

train_glmnetcr <- function(data) {
  XX <- data %>% dplyr::select(-Label) %>% as.matrix()

  mdl <- glmnetcr::glmnetcr(
    x = XX, y = data$Label, alpha = 0.01, thresh = 1e-04
  )
  preds <- predict(mdl)
  probs <- preds$probs[, 7, ]
  probs[, dim(probs)[2]]
}

train_ordinalridge <- function(data) {
  Y <- data[["Label"]]
  XX <- data %>% dplyr::select(-Label) %>% as.matrix()
  X %*% t(X) / ncol(X)
  # K <- XX %*% t(XX)
  K <- cov(t(XX))
  mdl <- ordinalRidge::ordinalRidge(
    K, Y, eps = 1e-04
  )
  predict(mdl, newdata = K)$prob[, "Pr[y >= 6]"]
}

timings <- microbenchmark(
  train_ordinalridge(traindata),
  # train_oforest(traindata),
  train_onet(traindata),
  train_pordlogist(traindata),
  # train_ordinalgmifs(traindata),
  train_glmnetcr(traindata),
  times = 3L, control = list(warmup = 1L)
)

set.seed(42)
res <- map(
  list(
    # oforest = train_oforest,
    onet = train_onet,
    pordlogist= train_pordlogist,
    glmnetcr = train_glmnetcr,
    oridge = train_ordinalridge
  ),
  ~.x(traindata)
)

res_df <- res %>%
  as_tibble() %>%
  mutate(Label = as.integer(traindata$Label))

res_cor <- cor(
  as.matrix(res_df), method = "kendall"
)

# Make random training datasets of 200 genes
all_genes <- setdiff(
  colnames(traindata_all),
  c("Barcode", "ID", "PMI", "AOD", "CDR", "Braak", "Label")
)
set.seed(42)
traindata_random <- map(
  seq_len(50),
  function(...) {
    traindata_all %>%
      dplyr::select(
        Label, all_of(sample(all_genes, size = 200, replace = FALSE))
      )
  }
)

res_all <- list(
  # oforest = train_oforest,
  onet = train_onet,
  pordlogist = train_pordlogist,
  glmnetcr = train_glmnetcr,
  ordinalridge = train_ordinalridge
) %>%
  enframe(name = "method", value = "fun") %>%
  crossing(
    traindata_random %>%
      enframe(name = "gene_set_id", value = "gene_set")
  ) %>%
  mutate(
    res = future_map2(
      fun, gene_set,
      ~possibly(.x, otherwise = NULL)(.y),
      .progress = TRUE,
      .options = furrr_options(
        seed = 42L
      )
    )
  )

cor_all <- res_all %>%
  group_by(method) %>%
  summarize(
    tau = cor(
      res %>%
        as.data.frame() %>%
        as.matrix() %>%
        unname(),
      as.integer(traindata_all$Label),
      method = "kendall"
    ) %>% {
      .[, 1]
    } %>%
      list(),
    rank_score = map_dbl(
      res,
      ~ordinalRidge::evaluate_ranking(
        as.integer(traindata_all$Label),
        .x
      )
    ) %>%
      list(),
    .groups = "drop"
  )

histogram_50_random <- cor_all %>%
  pivot_longer(c(tau, rank_score), names_to = "score_type", values_to = "score") %>%
  unchop(score) %>%
  ggplot(aes(score, color = method)) +
    geom_density(aes(y = stat(ndensity))) +
    theme_minimal() +
    labs(
      x = "Kendall correlation with true labels", y = "Count",
      title = "50 random gene sets (200 genes)"
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~score_type, ncol = 1)

dir.create(here("qc_plots"), showWarnings = FALSE)
ggsave(
  here("qc_plots", "tau_histogram_50_random_gene_sets.pdf"),
  histogram_50_random,
  width = 5, height = 3
)

rank_tau_cor_plot <- cor_all %>%
  unchop(c(tau, rank_score)) %>%
  ggplot(aes(tau, rank_score)) +
    geom_point() +
    facet_wrap(~method)

ggsave(
  here("qc_plots", "rank_score_vs_tau.pdf"),
  rank_tau_cor_plot,
  width = 5, height = 5
)

table_timings_50_random <- gridExtra::tableGrob(
  summary(timings, unit = "s") %>%
    mutate(across(where(is.numeric), signif, digits = 2)),
  rows = NULL,
  theme = gridExtra::ttheme_default()
)

ggsave(
  here("qc_plots", "tau_timings_table_50_random_gene_sets.pdf"),
  table_timings_50_random,
  width = 7, height = 2
)

res_all_hm_data <- res_all %>%
  dplyr::select(-fun, -gene_set, -gene_set_id) %>%
  mutate(
    res = map(
      res,
      ~tibble(prediction = .x, label = as.integer(traindata_all$Label)) %>%
        count(prediction, label)
    )
  ) %>%
  unnest(res) %>%
  group_by(method, prediction, label) %>%
  summarize(across(n, sum), .groups = "drop") %>%
  group_by(method, label) %>%
  mutate(prop_label = n / sum(n)) %>%
  ungroup()

res_all_hm <- res_all_hm_data %>%
  mutate(across(c(label, prediction), ~factor(.x, levels = as.character(1:7)))) %>%
  filter(!method %in% "ordinalridge") %>%
  ggplot(aes(label, prediction, fill = prop_label)) +
    geom_tile() +
    facet_wrap(vars(method)) +
    geom_text(aes(label = signif(prop_label, digits = 2)), color = "white") +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(x = "Label", y = "Prediction", fill = "Proportion of\nlabel (columns)")

ggsave(
  here("qc_plots", "tau_confusion_matrix_50_random.pdf"),
  res_all_hm,
  width = 13, height = 5
)
