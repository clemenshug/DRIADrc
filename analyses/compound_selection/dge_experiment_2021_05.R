library(tidyverse)
library(here)
library(safejoin)
library(fst)
library(readxl)
library(DRIAD)

synapser::synLogin()
syn <- synExtra::synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("analyses", "compound_selection")

hmsl_compounds_raw <- syn("syn22315918") %>%
  read_rds()

lspci_id_map <- syn("syn23748395") %>%
  read_csv()

drugs_available <- hmsl_compounds_raw %>%
  filter(n_batches > 0) %>%
  inner_join(
    lspci_id_map %>%
      select(lspci_id, id),
    by = c("hms_id" = "id")
  ) %>%
  pull(lspci_id) %>%
  unique()

fnROSMAP <- wrangleROSMAP(here("data", "rosmap"))
# fnMSBB <- Sys.glob(here::here("data-raw", "msbb*.tsv.gz"))
fnMSBB <- wrangleMSBB(here("data", "msbb"))

expression_data_amp_ad <- tribble(
  ~brain_region, ~dataset, ~path,
  "Dorsal prefrontal cortex", "ROSMAP", fnROSMAP,
  "BM10", "MSBB", fnMSBB[1],
  "BM22", "MSBB", fnMSBB[2],
  "BM36", "MSBB", fnMSBB[3],
  "BM44", "MSBB", fnMSBB[4],
) %>%
  mutate(
    data = map(path, read_tsv)
  )

valid_gene_symbols <- expression_data_amp_ad %>%
  pull(data) %>%
  map(colnames) %>%
  reduce(intersect) %>%
  setdiff(
    c("ID", "PMI", "AOD", "CDR",
      "Braak", "Barcode", "Label")
  )

expression_data_amp_ad_agg <- expression_data_amp_ad %>%
  mutate(
    data = map(
      data,
      ~select(.x, all_of(valid_gene_symbols)) %>%
        summarize(
          across(
            everything(),
            .fns = function(x) quantile(x, c(0.33, 0.67), names = FALSE) %>%
                           {.[order(abs(.), decreasing = TRUE)][1]}
          )
        )
    )
  )

expression_data_amp_ad_agg_formatted <- expression_data_amp_ad_agg %>%
  mutate(
    data = map(
      data,
      ~pivot_longer(.x, everything(), names_to = "symbol", values_to = "expression")
    )
  ) %>%
  select(brain_region, data) %>%
  unnest(data) %>%
  pivot_wider(id_cols = symbol, names_from = brain_region, values_from = expression) %>%
  rename_with(paste, "expression", sep = "_", .cols = -symbol)

commercial_info <- syn("syn25173589") %>%
  read_fst()

bbb_info <- read_xlsx(
  file.path(wd, "MCE_ScreeningLibrary_BBB_penetrating_HY-L028_2021-04-14.xlsx"),
  sheet = 2
) %>%
  left_join(
    commercial_info %>%
      distinct(lspci_id, catalog_number),
    by = c("Catalog Number" = "catalog_number")
  )

selectivity <- tibble(class = c("gpcrs", "nuclear_hormone_receptors")) %>%
  rowwise() %>%
  mutate(
    data = read_csv(file.path(wd, paste0(class, ".csv"))) %>% list()
  ) %>%
  unnest(data)

gtex <- read_tsv(
  file.path(wd, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"),
  skip = 2
)

gtex_brain <- gtex %>%
  select(Description, starts_with(fixed("Brain"))) %>%
  group_by(Description) %>%
  summarize(across(everything(), median), .groups = "drop")

gtex_brain_agg <- gtex_brain %>%
  pivot_longer(-Description, names_to = "tissue", values_to = "expression") %>%
  group_by(Description) %>%
  summarize(
    expression = quantile(expression, c(0.33, 0.67), names = FALSE) %>%
      {.[order(abs(.), decreasing = TRUE)][1]},
    .groups = "drop"
  )

master_table <- selectivity %>%
  safe_left_join(
    gtex_brain_agg %>%
      select(symbol = Description, gtex_brain_tpm = expression),
    by = "symbol",
    check = "bcvm"
  ) %>%
  safe_left_join(
    expression_data_amp_ad_agg_formatted,
    by = "symbol",
    check = "bcvm"
  ) %>%
  mutate(
    bbb_library = lspci_id %in% bbb_info$lspci_id,
    available_at_iccb = lspci_id %in% drugs_available
  )

write_csv(
  master_table %>%
    mutate(across(c(where(is.numeric), -lspci_id, -gene_id), signif, digits = 3)),
  file.path(wd, "master_table.csv")
)

drug_research <- read_csv("analyses/compound_selection/prioritized_drugs_manual_research.csv")

prioritized_drugs <- master_table %>%
  rowwise() %>%
  filter(
    gtex_brain_tpm > 0.1 |
      sum(
        c_across(ends_with("_expression")) > 1 | is.na(c_across(ends_with("_expression")))
      ) >= 1
  ) %>%
  ungroup() %>%
  arrange(
    symbol, desc(bbb_library), desc(max_phase), ontarget_ic50_q1
  ) %>%
  mutate(across(c(where(is.numeric), -lspci_id, -gene_id), signif, digits = 3)) %>%
  group_by(symbol) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(
    is_odd = as.integer(as.factor(symbol)) %>%
      magrittr::mod(2)
  ) %>%
  left_join(
    drug_research,
    by = "name"
  )
  # filter(reason != "natural_product")

write_csv(
  prioritized_drugs,
  file.path(wd, "prioritized_drugs_by_target.csv")
)

prioritized_drugs_agg <- prioritized_drugs %>%
  group_by(symbol) %>%
  mutate(rank = seq_len(n())) %>%
  ungroup() %>%
  group_by(
    chembl_id, name, max_phase, bbb_library, available_at_iccb
  ) %>%
  summarize(
    max_rank = min(rank),
    targets = paste(symbol, collapse = " "),
    affinities = paste(
      paste(symbol, ontarget_ic50_q1, sep = ":"),
      collapse = " "
    ), .groups = "drop"
  )

prioritized_drugs_agg_manual <- prioritized_drugs_agg %>%
  left_join(
    drug_research,
    by = "name"
  )
  # filter(reason != "natural_product")

write_csv(
  prioritized_drugs_agg_manual,
  file.path(wd, "prioritized_drugs_by_drug.csv")
)
