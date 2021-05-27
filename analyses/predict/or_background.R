# Run on O2

library(DRIAD)
library(tidyverse)
library(here)
library(batchtools)
library(qs)

wd <- here("analyses", "predict", "oridge")
dir.create(wd)

fnROSMAP <- wrangleROSMAP(file.path(wd, "rosmap"))
# fnMSBB <- Sys.glob(here::here("data-raw", "msbb*.tsv.gz"))
fnMSBB <- wrangleMSBB(file.path(wd, "msbb"))

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
      list(),
    pairs = preparePairs(task) %>%
      list()
  ) %>%
  select(-path) %>%
  ungroup() %>%
  mutate(
    id = paste0("prediction_task_", 1:n())
  )

qsave(
  prediction_tasks_all,
  file.path(wd, paste0("prediction_tasks_all", ".qs"))
)

# prediction_tasks_all <- qread(file.path(wd, paste0("prediction_tasks_all", ".qs")))

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

valid_gene_symbols <- prediction_tasks_all %>%
  pull(task) %>%
  map(colnames) %>%
  reduce(intersect) %>%
  setdiff(
    c("ID", "PMI", "AOD", "CDR",
      "Braak", "Barcode", "Label")
  )

gene_set_sizes <- c(
  5:29,
  seq(30, 300, by = 5)
)

set.seed(42)
prediction_tasks_gene_sets <- prediction_tasks_all %>%
  crossing(
    gene_set_size = gene_set_sizes,
    # Using 100 background gene sets per chunk and we want 1000, so we need 10
    # repeats
    gene_set_seq_id = seq_len(10)
  ) %>%
  rowwise() %>%
  mutate(
    background_sets = DRIAD:::genBK(
      valid_gene_symbols[1:gene_set_size],
      task,
      100
    ) %>%
      set_names(., paste0("BK_", seq_along(.))) %>%
      list()
  ) %>%
  ungroup() %>%
  select(-task, -pairs) %>%
  mutate(
    background_task_id = paste0("background_task_", seq_len(nrow(.)))
  )

qsave(
  prediction_tasks_gene_sets,
  file.path(wd, paste0("prediction_tasks_gene_sets", ".qs"))
)

# prediction_tasks_gene_sets <- qread(here("data", paste0("prediction_tasks_gene_sets", ".qs")))

# Set up jobs
reg <- makeRegistry(
  file.dir = file.path(wd, paste0("registry_", gsub(" ", "_", Sys.time()))),
  seed = 1
)
#reg$cluster.functions <- makeClusterFunctionsSlurm(template = "slurm-simple")

run_bk_job <- function(task_id, background_sets, background_task_id, ...) {
  library(tidyverse)
  library(furrr)
  library(DRIAD)
  library(qs)
  library(here)

  plan(sequential)

  message("Reading task...", task_id)
  task <- qread(here("data", paste0(task_id, ".qs")))
  message("Evaluating gene sets...", background_task_id)
  out <- tryCatch({
      DRIAD::evalGeneSets(
        background_sets,
        XY = task[["task"]],
        lP = task[["pairs"]],
        nBK = 0,
        method = "or"
      ) %>%
        dplyr::select(-Feats, -BK)
    },
    error = function(e) {
      NULL
    }
  )
  message("Done...")
  out
}

batchMap(
  fun = run_bk_job,
  task_id = prediction_tasks_gene_sets[["id"]],
  background_sets = prediction_tasks_gene_sets[["background_sets"]],
  background_task_id = prediction_tasks_gene_sets[["background_task_id"]]
)

# run_bk_job(
#   task_id = prediction_tasks_gene_sets[1:5,][["id"]][[1]],
#   background_sets = prediction_tasks_gene_sets[1:5,][["background_sets"]][[1]],
#   background_task_id = prediction_tasks_gene_sets[1:5,][["background_task_id"]][[1]]
# )

job_table <- findJobs() %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table[findExpired()],
  resources = list(
    memory = "8gb",
    ncpus = 1L,
    partition = "medium",
    walltime = 32*60*60,
    chunks.as.arrayjobs = TRUE,
    # For some reason these nodes fail to execute R because of an "illegal instruction"
    exclude = "compute-f-17-[09-25]"
  )
)

prediction_tasks_outputs_raw <- reduceResultsDataTable(
  findDone()
)
# prediction_tasks_outputs_raw <- reduceResultsDataTable(
#   findDone()
# )

qsave(
  prediction_tasks_outputs_raw,
  file.path(wd, "oridge_background_predictions_raw.qs")
)

prediction_tasks_outputs <- prediction_tasks_outputs_raw %>%
  filter(map_lgl(result, negate(is.null))) %>%
  unnest(result) %>%
  select(job_id = job.id, Set, AUC) %>%
  inner_join(
    prediction_tasks_gene_sets %>%
      mutate(job_id = seq_len(nrow(.))) %>%
      select(job_id, brain_region, dataset, gene_set_size, gene_set_seq_id),
    by = "job_id"
  )

qsave(
  prediction_tasks_outputs,
  file.path(wd, "oridge_background_predictions.qs")
)
