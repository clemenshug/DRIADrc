library(DRIAD)
library(tidyverse)
library(qs)
library(here)
library(microbenchmark)
library(ordinalForest)
library(ordinalNet)
library(OrdinalLogisticBiplot)
library(ordinalgmifs)
library(glmnetcr)
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
  predict(mdl, type = "class")
}

train_pordlogist <- function(data) {
  Y <- as.integer(data[["Label"]])
  XX <- data %>% dplyr::select(-Label) %>% as.matrix()

  preds <- OrdinalLogisticBiplot::pordlogist(
    Y, XX, penalization = 0.1, tol = 1e-04, maxiter = 200, show = FALSE
  )
  preds$pred[, 1]
}

train_ordinalgmifs <- function(data) {
  mdl <- ordinalgmifs::ordinalgmifs(
    Label~1, x = setdiff(names(data), "Label"), data = data, tol = 1e-04
  )
  preds <- predict(mdl)
  as.integer(preds$class)
}

train_glmnetcr <- function(data) {
  XX <- data %>% dplyr::select(-Label) %>% as.matrix()

  mdl <- glmnetcr::glmnetcr(
    x = XX, y = data$Label, alpha = 0.01, thresh = 1e-04
  )
  preds <- predict(mdl)
  as.integer(preds$class[, ncol(preds$class)]) + 1L
}

timings <- microbenchmark(
  train_oforest(traindata),
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
    glmnetcr = train_glmnetcr
  ),
  ~.x(traindata)
)

res_df <- as_tibble(res) %>%
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
  pordlogist= train_pordlogist,
  glmnetcr = train_glmnetcr
) %>%
  enframe(name = "method", value = "fun") %>%
  crossing(
    traindata_random %>%
      enframe(name = "gene_set_id", value = "gene_set")
  ) %>%
  mutate(
    res = future_map2(
      fun, gene_set,
      ~.x(.y),
      .progress = TRUE,
      .options = furrr_options(
        seed = 42L
      )
    )
  )

cor_all <- res_all %>%
  group_by(method) %>%
  summarize(
    cor_mat = cor(
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
    .groups = "drop"
  )

histogram_50_random <- cor_all %>%
  unchop(cor_mat) %>%
  ggplot(aes(cor_mat, color = method)) +
    geom_histogram(fill = NA) +
    theme_minimal() +
    labs(
      x = "Kendall correlation with true labels", y = "Count",
      title = "50 random gene sets (200 genes)"
    ) +
    scale_x_continuous(limits = c(0, 1))

dir.create(here("qc_plots"), showWarnings = FALSE)
ggsave(
  here("qc_plots", "tau_histogram_50_random_gene_sets.pdf"),
  histogram_50_random,
  width = 5, height = 3
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
