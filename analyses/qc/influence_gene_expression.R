library(DRIAD)
library(tidyverse)
library(here)
library(qs)
library(ordinalRidge)
library(furrr)

synapser::synLogin()
syn <- synExtra::synDownloader(here("data"))

data_dir <- here("analyses", "predict", "oridge")

wd <- here("analyses", "qc")

# Checking if gene expression in gene set influences predictive performance
# when using ordinalRidge

prediction_tasks_outputs <- file.path(data_dir, "oridge_background_predictions.qs") %>%
  qread()

prediction_tasks_gene_sets <- file.path(data_dir, "prediction_tasks_gene_sets.qs") %>%
  qread()

predictions_tasks <- file.path(data_dir, "prediction_tasks_all.qs") %>%
  qread()

prediction_tasks_gene_sets_long <- prediction_tasks_gene_sets %>%
  mutate(across(background_sets, map, enframe, name = "background_set_id", value = "genes")) %>%
  unnest(background_sets)

expression_data <- with(
  predictions_tasks,
  set_names(
    task %>%
      map(select, where(is.numeric)) %>%
      map(as.matrix),
    brain_region
  )
)

prediction_tasks_gene_sets_expression <- prediction_tasks_gene_sets_long %>%
  mutate(
    avg_expression = pmap_dbl(
      list(genes, brain_region),
      function(genes, brain_region) {
        expression_data[[brain_region]][, genes] %>%
          mean()
      }
    )
  )


prediction_tasks_gene_sets_expression_plot_data <- prediction_tasks_gene_sets_expression %>%
  select(brain_region, dataset, gene_set_size, gene_set_seq_id, background_set_id, avg_expression) %>%
  inner_join(
    select(prediction_tasks_outputs, brain_region, dataset, gene_set_size, background_set_id = Set, gene_set_seq_id, AUC)
  )

prediction_tasks_gene_sets_expression_plot <- prediction_tasks_gene_sets_expression_plot_data %>%
  filter(gene_set_size %in% c(5, 10, 20, 40, 80, 160, 300)) %>%
  ggplot(aes(avg_expression, AUC)) +
    geom_hex(aes(fill = after_stat(ndensity))) +
    geom_smooth(method = "lm") +
    facet_grid(vars(gene_set_size), vars(brain_region)) +
    theme_minimal() +
    labs(x = "Average gene expression", y = "AUC", fill = "Density")

ggsave(
  file.path(wd, "expression_vs_auc_density.pdf"),
  prediction_tasks_gene_sets_expression_plot,
  width = 7, height = 8
)
