library(DRIAD)
library(tidyverse)
library(here)
library(qs)
library(ordinalRidge)

synapser::synLogin()
syn <- synExtra::synDownloader(here("data"))

run_gs_job <- function(task_id, gene_set, directory, ...) {
  library(tidyverse)
  library(furrr)
  library(DRIAD)
  library(qs)

  plan(sequential)

  message("Reading task...", task_id)
  task <- qread(file.path(directory, paste0(task_id, ".qs")))
  out <- tryCatch({
    DRIAD::evalGeneSet(
      intersect(gene_set, colnames(task[["task"]])),
      XY = task[["task"]],
      lP = task[["pairs"]],
      nBK = 0,
      method = "or"
    ) %>%
      dplyr::select(-Feats)
    },
    error = function(e) {
      NULL
  })
  message("Done...")
  out
}

wd <- here("analyses", "predict", "oridge")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)

fnROSMAP <- wrangleROSMAP(here("data", "rosmap"))
# fnMSBB <- Sys.glob(here::here("data-raw", "msbb*.tsv.gz"))
fnMSBB <- wrangleMSBB(here("data", "msbb"))

prediction_tasks_all <- tribble(
  ~brain_region, ~dataset, ~path,
  "Dorsal prefrontal cortex", "ROSMAP", fnROSMAP,
  "BM10", "MSBB", fnMSBB[1],
  "BM22", "MSBB", fnMSBB[2],
  "BM36", "MSBB", fnMSBB[3],
  "BM44", "MSBB", fnMSBB[4],
) %>%
  crossing(
    comparison = c("all")
  ) %>%
  rowwise() %>%
  mutate(
    task = prepareTask(path, comparison) %>%
      list()
  ) %>%
  select(-path) %>%
  ungroup() %>%
  mutate(
    id = paste0("prediction_task_", 1:n())
  )

dir.create(file.path(wd, "data"))
qsave(
  prediction_tasks_all,
  file.path(wd, "data", paste0("prediction_tasks_all", ".qs"))
)

# prediction_tasks_all <- qread(file.path(wd, "data", paste0("prediction_tasks_all", ".qs")))

pwalk(
  prediction_tasks_all,
  function(task, pairs, id, ...) {
    # browser()
    x <- list(task = task, pairs = pairs)
    qsave(
      x,
      file.path(wd, paste0(id, ".qs"))
      # compress = "xz"
    )
  }
)

prediction_tasks <- prediction_tasks_all %>%
  select(-task, -pairs)

dge_gmt <- map(
  here(c("DGE1.gmt", "DGE2.gmt")) %>%
    set_names(c("DGE1", "DGE2")),
  read_gmt
) %>%
  enframe("experiment", "gene_sets") %>%
  mutate(
    gene_sets = map(gene_sets, ~enframe(.x, "drug", "genes"))
  ) %>%
  unnest(gene_sets)

# from or_background.R
background_predictions <- qread(file.path(wd, "oridge_background_predictions_probabilities.qs"))

background_gene_set_sizes <- background_predictions[["gene_set_size"]] %>%
  unique()

library(furrr)
plan(multisession(workers = 20))
# plan(multicore(workers = 8))
# plan(sequential)

dge_tasks <- dge_gmt %>%
  crossing(
    prediction_tasks
  ) %>%
  mutate(
    gene_set_size = map_int(genes, length),
    background_gene_set_size = cut(
      gene_set_size,
      breaks = background_gene_set_sizes,
      labels = FALSE
    ) %>% {
      background_gene_set_sizes[. + 1L]
    },
    genes = map2(
      genes, background_gene_set_size,
      ~head(.x, n = .y)
    )
  )

dge_res <- dge_tasks %>%
  mutate(
    res = future_map2(
      id, genes,
      run_gs_job,
      directory = wd,
      .progress = TRUE,
      .options = furrr_options(
        seed = 42,
        globals = FALSE
      )
    )
  )

qsave(
  dge_res,
  file.path(wd, "dge_res_raw.qs")
)

# dge_res <- qread(file.path(wd, "dge_res_raw.qs"))

dge_res_neat <- dge_res %>%
  dplyr::select(-genes) %>%
  unnest(res)

prediction_task_labels <- prediction_tasks_all %>%
  with(
    set_names(map(task, "Label"), brain_region)
  )

M <- syn("syn11801537") %>% read_csv(col_types=cols()) %>%
  mutate_at( "name", str_to_lower ) %>%
  dplyr::select( LINCSID = lincs_id, URL = link, Drug = name )

dge_res_p <- dge_res_neat  %>%
  inner_join(
    background_predictions %>%
      group_by(brain_region, dataset, gene_set_size) %>%
      summarize(background_scores = list(AUC), .groups = "drop"),
    by = c("brain_region", "dataset", "background_gene_set_size" = "gene_set_size")
  ) %>%
  mutate(
    p = map2_dbl(
      background_scores, AUC,
      ~mean(.x > .y)
    )
  ) %>%
  dplyr::select(-background_scores) %>%
  bind_rows(
    dplyr::select(., experiment, drug, brain_region, p) %>%
      pivot_wider(id_cols = c(experiment, drug), names_from = brain_region, values_from = p) %>%
      rowwise() %>%
      mutate(across(where(is.numeric), pmax, 0.0005)) %>%
      mutate(p = 1 / mean(1 / c_across(where(is.numeric)))) %>%
      ungroup() %>%
      dplyr::transmute(experiment, drug, brain_region = "aggregated", p)
  ) %>%
  left_join(
    dplyr::select(M, LINCSID, drug_name = Drug),
    by = c("drug" = "LINCSID")
  )

write_csv(
  dge_res_p,
  file.path(wd, "dge_p_values.csv.gz")
)

top_p_heatmap <- dge_res_p %>%
  mutate(
    drug_name = paste0(drug_name, " (", str_sub(experiment, start = -1L, end = -1L), ")"),
    p_log = -log10(p + 0.5 * min(p[p > 0])),
    brain_region
  ) %>%
  # arrange(desc(p_log)) %>%
  # arrange(p) %>%
  mutate(
    drug_name = factor(
      drug_name,
      levels = filter(., brain_region == "aggregated") %>% arrange(desc(p_log)) %>% pull(drug_name) %>% rev()
    )
  ) %>%
  filter(as.integer(drug_name) > 30) %>%
  ggplot(aes(brain_region, drug_name, fill = p_log)) +
  geom_tile() +
  # scale_fill_gradientn(
  #   colors = c("white", "yellow", "blue", "red"),
  #   # breaks = scales::rescale(rev(-log10(c(0.001, 0.049, 0.05, 1)))),
  #   breaks = c(0, 0.1, 0.11, 1),
  #   # limits = rev(-log10(c(0.001, 1))),
  #   # breaks = scales::rescale(rev(-log10(c(0.0001, 0.049, 0.05, 1)))),
  #   # breaks = c(0.00001, 0.049, 0.05, 1),
  #   # trans = "log10",
  #   oob = scales::squish
  # )
  scale_fill_viridis_c()

ggsave(
  file.path(wd, "top_p_heatmap.pdf"),
  top_p_heatmap, width = 6, height = 12
)

p_histogram <- dge_res_p %>%
  ggplot(aes(p)) +
  geom_histogram() +
  facet_wrap(vars(brain_region))

ggsave(
  file.path(wd, "p_histogram.pdf"),
  p_histogram, width = 5, height = 3
)

dge_res_w_bk <- dge_res_cor %>%
  mutate(
    background_gene_set_size = cut(
      gene_set_size,
      breaks = background_gene_set_sizes,
      labels = FALSE
    ) %>% {
      background_gene_set_sizes[. + 1L]
    }
  ) %>%
  inner_join(
    background_predictions_tau %>%
      drop_na() %>%
      group_by(brain_region, dataset, gene_set_size) %>%
      summarize(background_taus = list(tau), .groups = "drop"),
    by = c("brain_region", "dataset", "background_gene_set_size" = "gene_set_size")
  ) %>%
  mutate(
    p = map2_dbl(
      background_taus, tau,
      ~mean(.x > .y)
    )
  ) %>%
  left_join(
    dplyr::select(M, LINCSID, drug_name = Drug),
    by = c("drug" = "LINCSID")
  )

drugs = c("ruxolitinib", "nilotinib")

example_enrichment_plot_data <- dge_res_w_bk %>%
  filter(drug_name %in% drugs) %>%
  mutate(
    drug_name = paste0(drug_name, " (", str_sub(experiment, start = -1L, end = -1L), ")"),
    brain_region = as.factor(brain_region)
  )

library(ggridges)

## Generate the ridge plots
example_enrichment_plot <- example_enrichment_plot_data %>%
  mutate(
    drug_name = paste0(drug_name, " (", gene_set_size, ")")
  ) %>%
  dplyr::select(drug_name, brain_region, tau = background_taus) %>%
  unchop(tau) %>%
  ggplot(aes(x=tau, y=brain_region, fill=brain_region) ) +
  facet_wrap( ~drug_name, nrow=1 ) +
  theme_ridges(center_axis_labels=TRUE) +
  geom_vline( xintercept=seq(0.1, 0.9, by=0.1), color="gray90" ) +
  geom_density_ridges2(scale=1.25, size=1, alpha=0.5) +
  geom_segment( aes(x=tau, xend=tau, y=as.numeric(brain_region),
                    yend=as.numeric(brain_region)+0.9),
                data=example_enrichment_plot_data %>%
                  mutate(
                    drug_name = paste0(drug_name, " (", gene_set_size, ")")
                  ) %>%
                  dplyr::select(drug_name, brain_region, tau),

                color="red", lwd=1 ) +
  coord_cartesian(clip="off") + xlim(0, 1) +
  # scale_fill_manual( values=dsPal(), guide=FALSE ) +
  # scale_color_manual( values=c("yes"="red","no"="black"), guide=FALSE ) +
  theme( strip.text.x = element_text(margin=margin(b=5,t=5), face="bold"),
         strip.background = element_blank(), panel.grid.major.x=element_blank() )

ggsave(
  file.path(wd, "example_enrichment_plot.pdf"),
  example_enrichment_plot, width = 14, height = 4
)
