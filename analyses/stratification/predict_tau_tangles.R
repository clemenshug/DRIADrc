library(DRIAD)
library(tidyverse)
library(here)
library(qs)
library(sva)

synapser::synLogin()
syn <- synExtra::synDownloader(here("data"))

data_dir <- here("analyses", "predict", "oridge")

wd <- here("analyses", "stratification")

genecode <- syn("syn25813317") %>%
  read_csv()

# wrangle_raw_counts <- function() {
data_syn <- tribble(
  ~id, ~dataset,
  "syn8691099", "MSBB",
  "syn8691134", "ROSMAP",
  "syn8690799", "MAYO",
  "syn8690904", "MAYO"
)
data_files <- data_syn %>%
  mutate(data = map(id, ~read_tsv(syn(.x))))

data_long <- data_files %>%
  mutate(data = map(data, pivot_longer, cols = -feature, names_to = "sample", values_to = "count")) %>%
  unnest(data) %>%
  filter(str_starts(feature, "ENSG"))
# }

qsave(
  data_long,
  file.path(wd, "expression_data_long.qs")
)
# data_long <- qread(file.path(wd, "expression_data_long.qs"))

# wrangle_metadata <- function()
rosmap_specimen <- syn("syn21323366") %>%
  read_csv(
    col_types = cols(
      specimenIdSource = col_character(),
      tissueVolume = col_double(),
      samplingAge = col_character(),
      excludeReason = col_character()
    )
  )
rosmap_individual <- syn("syn3191087") %>%
  read_csv()
mayo_specimen <- syn("syn20827192") %>%
  read_csv()
mayo_individual <- syn("syn23277389") %>%
  read_csv()
mayo_additional <- syn("syn14031984") %>%
  read_tsv()
msbb_specimen <- syn("syn21893059") %>%
  read_csv(
    col_types = cols(
      cellType = col_character(),
      specimenIdSource = col_character(),
      sampleStatus = col_character(),
      nucleicAcidSource = col_character()
    )
  )
msbb_individual <- syn("syn6101474") %>%
  read_csv()
# }

meta_batch <- bind_rows(
  transmute(
    rosmap_specimen,
    dataset = "ROSMAP",
    specimen_id = specimenID,
    brain_area = "dorsolateral prefrontal cortex",
    exclude
  ),
  transmute(
    mayo_specimen,
    dataset = "MAYO",
    specimen_id = specimenID,
    brain_area = tissue,
    exclude
  ),
  transmute(
    msbb_specimen,
    dataset = "MSBB",
    specimen_id = specimenID,
    brain_area = tissue,
    exclude
  )
) %>%
  replace_na(list(exclude = FALSE)) %>%
  drop_na() %>%
  mutate(
    sample_id = paste(dataset, specimen_id, sep = "_"),
    batch = factor(paste(dataset, brain_area, sep = "_"))
  ) %>%
  distinct()

expression_mat_raw <- data_long %>%
  transmute(sample = paste(dataset, sample, sep = "_"), feature, count) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  column_to_rownames("feature") %>%
  as.matrix()

qsave(
  expression_mat_raw,
  file.path(wd, "expression_mat_raw.qs")
)

genecode_protein_coding <- genecode %>%
  filter(gene_type == "protein_coding") %>%
  distinct(gene_id, gene_type, gene_name)

meta_batch_include <- meta_batch %>%
  filter(!exclude, sample_id %in% colnames(expression_mat_raw)) %>%
  mutate(across(where(is.factor), fct_drop))

expression_mat_include <- expression_mat_raw[
  intersect(genecode_protein_coding$gene_id, rownames(expression_mat_raw)),
  meta_batch_include$sample_id
]

expression_mat_corrected <- expression_mat_include %>%
  ComBat_seq(batch = meta_batch_include$batch)

qsave(
  expression_mat_corrected,
  file.path(wd, "expression_mat_corrected")
)
# expression_mat_corrected <- qread(file.path(wd, "expression_mat_corrected"))

library(DESeq2)
library(broom)

meta_mayo_tdp <- mayo_additional %>%
  drop_na(`TDP-43`) %>%
  mutate(
    sample_id = paste("MAYO", NETdbID, sep = "_"),
    tau_tangles = factor(`TDP-43`, labels = c("no", "yes")),
    stage = cut(`Braak stage`, breaks = c(-Inf, 2, 4, Inf), labels = c("A", "B", "C"), ordered_result = TRUE),
    `Thal phase` = ordered(`Thal phase`)
  ) %>%
  filter(sample_id %in% colnames(expression_mat_corrected))

ggsave(
  file.path(wd, "mayo_meta_associations.pdf"),
  GGally::ggpairs(select(meta_mayo_tdp, Gender, AgeAtDeath, tau_tangles, stage)), width = 10, height = 9
)

expression_mat_corrected_tdp <- expression_mat_corrected[
  , meta_mayo_tdp$sample_id
]

des <- DESeqDataSetFromMatrix(
  expression_mat_corrected_tdp, mutate(meta_mayo_tdp, across(stage, as.character)),
  design = ~tau_tangles
)

des <- DESeq(des)

des_vst <- varianceStabilizingTransformation(des)

plotPCA(des_vst, intgroup = c("tau_tangles"))

all_genes_pca <- assay(des_vst)[
  order(rowVars(assay(des_vst)), decreasing = TRUE)[1:500],
] %>%
  t() %>%
  prcomp()

all_genes_pca_scores <- all_genes_pca %>%
  tidy(matrix = "scores") %>%
  dplyr::rename(sample_id = row)

all_genes_pca_pcs <- all_genes_pca %>%
  tidy(matrix = "pcs")

all_genes_pca_plot_data <- all_genes_pca_scores %>%
  pivot_wider(id = sample_id, names_from = PC, values_from = value, names_prefix = "PC") %>%
  left_join(
    meta_mayo_tdp
  )

add_margin <- function(p, type = c("density", "scatter"), ratio = 6) {
  library(patchwork)
  margin_impl <- function(p, axis, type) {
    out <- ggplot(p$data, aes(!!p$mapping[[axis]])) +
      switch(
        type,
        density = geom_density(aes(color = !!p$mapping$colour)),
        scatter = list(geom_point(aes(y = !!p$mapping$colour)),
          geom_smooth(aes(y = !!p$mapping$colour), method = "lm"))
      ) +
      guides(color = "none") +
      theme_minimal()
    if (axis == "y")
      out <- out + coord_flip() +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank())
    else
      out <- out + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    out
  }
  patchworkGrob(
    margin_impl(p, "x", type) + plot_spacer() + p + margin_impl(p, "y", type) +
      plot_layout(byrow = TRUE, widths = c(ratio, 1), heights = c(1, ratio))
  )
}

all_genes_pca_plot <- map(
  c("Gender", "AgeAtDeath", "stage", "tau_tangles", "Thal phase"),
  ~{ggplot(
    all_genes_pca_plot_data,
    aes(PC1, PC2, color = !!sym(.x))
  ) +
    geom_point() +
    theme_minimal()} %>%
    add_margin(type = if (is.numeric(all_genes_pca_plot_data[[.x]])) "scatter" else "density", ratio = 5)
    # ggExtra::ggMarginal(type = "density", groupColour = !is.numeric(all_genes_pca_plot_data[[.x]]))
) %>%
  {gridExtra::arrangeGrob(grobs = .)}

ggsave(
  file.path(wd, "pca_mayo_corrected_all_genes_top500.pdf"),
  all_genes_pca_plot, width = 12, height = 10
)

tau_tangle_de <- results(des, name = "tau_tangles_yes_vs_no") %>%
  as.data.frame() %>%
  as_tibble(rownames = "gene_id") %>%
  left_join(distinct(genecode, gene_id, gene_name))

tau_tangle_pca <- assay(des_vst)[
  tau_tangle_de %>%
    arrange(padj) %>%
    pull(gene_id) %>%
    head(n = 500),
] %>%
  t() %>%
  prcomp()

tau_tangle_pca_scores <- tau_tangle_pca %>%
  tidy(matrix = "scores") %>%
  dplyr::rename(sample_id = row)

tau_tangle_pca_pcs <- tau_tangle_pca %>%
  tidy(matrix = "pcs")

tau_tangle_pca_plot_data <- tau_tangle_pca_scores %>%
  pivot_wider(id = sample_id, names_from = PC, values_from = value, names_prefix = "PC") %>%
  left_join(
    meta_mayo_tdp
  )

tau_tangle_pca_plot <- map(
  c("Gender", "AgeAtDeath", "stage", "tau_tangles", "Thal phase"),
  ~{ggplot(
    tau_tangle_pca_plot_data,
    aes(PC1, PC2, color = !!sym(.x))
  ) +
      geom_point() +
      theme_minimal()} %>%
    add_margin(type = if (is.numeric(tau_tangle_pca_plot_data[[.x]])) "scatter" else "density", ratio = 5)
  # ggExtra::ggMarginal(type = "density", groupColour = !is.numeric(all_genes_pca_plot_data[[.x]]))
) %>%
  {gridExtra::arrangeGrob(grobs = .)}

ggsave(
  file.path(wd, "pca_mayo_corrected_tau_de_genes_no_adjustment.pdf"),
  tau_tangle_pca_plot, width = 12, height = 10
)

tau_tangle_all_vs_tau_de_data <- bind_rows(
  top500_var = all_genes_pca_scores,
  tau_tangle_de = tau_tangle_pca_scores,
  .id = "gene_set"
)

tau_tangle_all_vs_tau_de_plot <- tau_tangle_all_vs_tau_de_data %>%
  filter(PC == 1) %>%
  pivot_wider(id_cols = c(PC, sample_id), names_from = gene_set, names_prefix = "PC1_", values_from = value) %>%
  inner_join(
    meta_mayo_tdp
  ) %>%
  ggplot(aes(PC1_top500_var, PC1_tau_tangle_de, color = tau_tangles)) +
    geom_point()

## Modeling

fit_glmnet_model <- function(mat, meta) {
  varstab_count_mayo_long <- mat %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "count")

  set.seed(42)
  tdp_model_cv <- cv.glmnet(
    scale(t(mat)), meta$tau_tangles,
    family = "binomial", type.measure = "class",
    nfolds = 5, alpha = 0.5
  )

  tdp_model_coefs <- coef(tdp_model_cv, s="lambda.1se") %>%
    {
      dplyr::mutate(
        as.data.frame(summary(.)),
        gene_id = dimnames(.)[[1]][i]
      )
    } %>%
    left_join(
      distinct(genecode, gene_id, gene_name)
    ) %>%
    arrange(abs(x)) %>%
    mutate(across(gene_name, fct_inorder))

  top_coef_plot_data <- tdp_model_coefs %>%
    inner_join(
      varstab_count_mayo_long,
      by = "gene_id"
    ) %>%
    inner_join(
      select(meta, sample_id, tau_tangles),
      by = "sample_id"
    )
  top_coef_plot <- ggplot(top_coef_plot_data, aes(gene_name, count, color = tau_tangles)) +
    geom_boxplot(outlier.color = NA) +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.9, alpha = 0.5)

  list(
    model = tdp_model_cv,
    coefs = tdp_model_coefs,
    coef_plot_data = top_coef_plot_data,
    coef_plot = top_coef_plot
  )
}


des <- DESeqDataSetFromMatrix(
  expression_mat_include[
    , meta_mayo_tdp$sample_id
  ], meta_mayo_tdp,
  design = ~Gender + stage
)

des <- estimateSizeFactors(des)

varstab_count_mat_only_mayo <- varianceStabilizingTransformation(des) %>%
  assay() %>% {
    .[
      # Exclude low-expressed genes
      apply(., 1, mean) - min(.) > 0.5,
    ]
  }

des <- DESeqDataSetFromMatrix(expression_mat_corrected_tdp, meta_batch_include, design = ~batch)

des <- estimateSizeFactors(des)

norm_count_mat <- counts(des, normalized = TRUE)

varstab_count_mat <- varianceStabilizingTransformation(des) %>%
  assay()

varstab_count_mat_mayo <- varstab_count_mat[
  # Exclude low-expressed genes
  apply(varstab_count_mat, 1, mean) - min(varstab_count_mat) > 0.5,
  meta_mayo_tdp$sample_id
]




  # geom_text(
  #   aes(label = signif(x, digits = 2), color = NULL), y = 11,
  #   data = tdp_model_coefs
  # )

fit_xgb_model <-

library(xgboost)

set.seed(42)
# tdp_xgb_model_cv <- xgb.cv(
#   data = scale(t(varstab_count_mat_mayo)), label = as.integer(meta_mayo_tdp$tau_tangles) - 1L,
#   objective = "binary:logistic",max_depth = 2,
#   family = "binomial", type.measure = "class", stratified = TRUE,
#   nfold = 5, nrounds = 10, metrics = c("error", "auc", "aucpr"),
#   save_models = TRUE, prediction = TRUE
# )
xgb_params <- list(
  max_depth = 2, nthread = 6, lambda = 0.1, alpha = 0.9
)

tdp_xgb_model_cv <- xgb.cv(
  data = t(varstab_count_mat_mayo), label = as.integer(meta_mayo_tdp$tau_tangles) - 1L,
  params = xgb_params,
  nfold = 5, stratified = TRUE, metrics = c("error", "auc"),
  family = "binomial", type.measure = "class", early_stopping_rounds = 3,
  nrounds = 20, objective = "binary:logistic", prediction = TRUE
)

tdp_xgb_model_best <- xgboost(
  data = t(varstab_count_mat_mayo), label = as.integer(meta_mayo_tdp$tau_tangles) - 1L,
  eval_metric = "auc", nrounds = tdp_xgb_model_cv$best_iteration, objective = "binary:logistic",
  params = xgb_params
)

tdp_xgb_model_importance <- xgb.importance(model = tdp_xgb_model_best)
xgb.plot.importance(importance_matrix = tdp_xgb_model_importance)

top_coef_plot_data <- tdp_xgb_model_importance %>%
  mutate(
    gene_id = Feature
  ) %>%
  left_join(
    distinct(genecode, gene_id, gene_name)
  ) %>%
  inner_join(
    varstab_count_mayo_long,
    by = "gene_id"
  ) %>%
  inner_join(
    select(meta_mayo_tdp, sample_id, tau_tangles),
    by = "sample_id"
  ) %>%
  mutate(across(gene_name, fct_reorder, Importance))

top_coef_plot <- ggplot(top_coef_plot_data, aes(gene_name, count, color = tau_tangles)) +
  geom_boxplot(outlier.color = NA) +
  ggbeeswarm::geom_quasirandom(dodge.width = 0.9, alpha = 0.5)

ggsave(file.path(wd, "xgb_model_best_only_mayo_importance.pdf"), width = 20, height = 6)

## Predict

y <- predict(x, t(varstab_count_mat_mayo), type = "class")

set.seed(123)

varstab_count_mat_mayo_df <- varstab_count_mat_mayo %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  as_tibble() %>%
  inner_join(
    select(meta_mayo_tdp, sample_id, tau_tangles)
  )

library(tidymodels)
mayo_df_split <- initial_split(
  varstab_count_mat_mayo_df %>% select(-sample_id),
  strata = tau_tangles
)

mayo_df_split_train <- training(mayo_df_split)
mayo_df_split_test  <- testing(mayo_df_split)

log_model_spec <- logistic_reg(penalty = 1) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

log_model_recipe <- recipe(head(mayo_df_split_train)) %>%
  update_role(tau_tangles, new_role = "outcome") %>%
  update_role(c(everything(), - tau_tangles), new_role = "predictor")

log_model_wf <- workflow() %>%
  add_model(log_model_spec) %>%
  add_recipe(log_model_recipe)

log_model_fit <- log_model_wf %>%
  fit(data = mayo_df_split_train)

log_model_pred <- log_model_fit %>%
  predict(mayo_df_split_test, type = "prob")


