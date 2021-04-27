library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader("~/data/AMP-AD/DGE")

# Wrangle DGE metadata for GEO/SRA submission ----------------------------------
###############################################################################T

# Created in analyses/mechanism/04_dge_differential_expression.R

dge_meta <- here("results", "deseq_meta.csv") %>%
  read_csv()

# Clean up Drug field and split it into drug and treatment column

# Predetermined fields:
# *sample_name	sample_title	bioproject_accession	*organism	*isolate	*age	*biomaterial_provider	*sex	*tissue
# cell_line	cell_subtype	cell_type	culture_collection	dev_stage	disease	disease_stage
# ethnicity	health_state	karyotype	phenotype	population	race	sample_type	treatment	description

# In treatment, putting both drug and stimulation there
# Then have separate columns for drug and stimulation

dge_cleaned <- dge_meta %>%
  # Extract stimulant from Drug column
  extract(
    Drug, "stimulant",
    regex = "\\(([[:alnum:]]+)\\)",
    remove = FALSE
  ) %>%
  mutate(
    # Is the drug column actually a stimulant
    drug_is_stimulant = Drug %in% c("dsrna +lipo", "lps", "naked dsrna", "lipo control", "dsrnami"),
    # If drug is a stimulant, remove
    drug = if_else(drug_is_stimulant, NA_character_, Drug) %>%
      # Remove stimulants that were used in addition to drugs from Drug column
      str_replace("\\s*\\(([[:alnum:]]+)\\)", ""),
    stimulant = if_else(
      drug_is_stimulant,
        Drug,
        stimulant
      ) %>%
      recode(
        lipo = "lipofectamine",
        dsrnami = "dsRNA",
        `dsrna +lipo` = "dsRNA + lipofectamine",
        lps = "lipopolysaccharide",
        `naked dsrna` = "dsRNA",
        `lipo control` = "lipofectamine"
      ),
    # Unite drug and stimulant for treatment column
    treatment = case_when(
      is.na(drug) ~ NA_character_,
      drug == "dmso" ~ "dmso",
      TRUE ~ paste(Concentration, "uM", drug, sep = " ")
    ),
    drug_concentration = if_else(
      drug == "dmso",
      NA_real_,
      Concentration
    )
  ) %>%
  unite(treatment, treatment, stimulant, sep = " + ", na.rm = TRUE, remove = FALSE)

dge_geo <- dge_cleaned %>%
  group_by(experiment, treatment) %>%
  mutate(
    replicate = seq_len(n())
  ) %>%
  ungroup() %>%
  transmute(
    sample_name = Sample,
    sample_title = paste(
      experiment, treatment, "replicate", replicate, sep = " "
    ),
    organism = "Homo sapiens",
    isolate = "not applicable",
    age = "fetal",
    biomaterial_provider = "not applicable",
    sex = "male",
    tissue = "ventral mesencephalon",
    cell_line = "ReNcell VM",
    cell_type = "in-vitro differentiated neural and glial cells",
    sample_type = "mixed cell culture",
    treatment,
    drug_treatment = drug,
    drug_concentration,
    stimulant,
    replicate = paste("biological replicate", replicate),
    batch = experiment
  )

write_tsv(
  dge_geo,
  here("results", "dge_meta_geo.tsv")
)

design_description_text <- "A multiwell cell dispenser (catalog# 5840300, Thermo Scientific, Waltham, MA)
with standard tubing (catalog# 24072670, Thermo Scientific, Waltham, MA) was
used to plate 2500 neural stem cells (ReNcell VM, catalog# SCC008, Millipore, Billerica, MA)
into each well of a 384-well cell culture plate (Perkin Elmer, Waltham, MA).
Neural stem cells where differentiated into mature neural cells for one week and
then treated with compounds or DMSO using a D300 Digital Dispenser
(Hewlett-Packard, Palo Alto, CA). D300 software was used to randomize dispensation
of compounds. After 24 hours, the cells were washed once with PBS using an EL405x
plate washer (BioTek, Winooski, VT) leaving 5-10 ul of PBS behind in each well.
10 ul of 1X TCL lysis buffer (catalog# 1070498, Qiagen, Hilden, Germany) with
1% (v/v) ß-mercaptoethanol was added per well, and the plates were stored at -80°C
until the RNA extraction was performed.

For RNA extraction, the cell lysate plate was thawed and centrifuged for 1 min at
1000 rpm. Using a BRAVO (Agilent, Santa Clara, CA) liquid handler, the lysate
was mixed thoroughly before transferring 10 ul to a 384 well PCR plate. 28 ul
of home-made Serapure SPRI beads (GE Healthcare Life Sciences, Marlborough, MA)
were added directly to the lysate, mixed and incubated for 5 min. The plate was
transferred to a magnetic rack and incubated for 5 min prior to removing the
liquid to aggregate the beads. The beads were washed with 80% ethanol twice,
allowed to dry for 1 min, 20 ul of nuclease free water was added per well, the
plate was removed from the magnetic rack and the beads were thoroughly
resuspended. Following a 5 min incubation, the plate was returned to the
magnetic rack and incubated an additional 5 min before transferring the
supernatant to a fresh PCR plate.

5 ul of the RNA was transferred to a separate
plate containing RT master mix and 3’ and 5’ adapters for reverse transcription
and template switching, and incubated for 90 min at
42°C. The cDNA was pooled and purified with a QIAquick PCR purification kit
according to the manufacturer’s directions with the final elution in 21 ul of
nuclease free water. This was followed by exonuclease I treatment for 30 min at
37°C that was stopped with a 20 min incubation at 80°C. The cDNA was then
amplified using the Advantage 2 PCR Enzyme System (Takara, Fremont, CA) for 6
cycles, and purified using AMPure XP magnetic beads (Beckman Coulter Genomics, Chaska, MN).
Library preparation was performed using a Nextera XT DNA kit (Illumina, San Diego, CA)
on 5 reactions per sample following the manufacturer’s instructions, amplified
12 cycles, and purified with AMPure XP magnetic beads (Beckman Coulter Genomics, Chaska, MN).
The sample was then quantified by qPCR and sequenced on a single Illumina NextSeq
run with 75bp paired end reads at the Harvard University Bauer Core Facility.

The raw fastq files were processed using the bcbio-nextgen single cell/DGE RNA-seq
analysis pipeline (https://bcbio-nextgen.readthedocs.io/). Before submission,
the fastq file was demultiplexed into separate files using the sample specific
barcodes in the first read and unique molecular identifiers were used to filter
duplicate reads of the same original mRNA molecule.
"
  # str_replace_all("(?<!\n)\n(?!\n)", " ") %>%
  # str_replace_all("\n{2}", " ")

well_map <- syn("syn17103737") %>%
  read_tsv() %>%
  transmute(
    well = recode(well, P11_new1 = "P11"),
    barcode
  )

# Split multiplexed fastq file using command
# /n/app/bcbio/dev/anaconda/bin/python /n/app/bcbio/dev/anaconda/bin/umis demultiplex_cells --out_dir demultiplexed --readnumber 1 --prefix dge_1_ work/umis/ad_dge_2.filtered.fq.gz

dge_sra <- dge_geo %>%
  transmute(
    sample_name,
    library_ID = paste(sample_name, "rna_seq", sep = "_"),
    title = paste("RNA-seq of", cell_line, "treated with", treatment),
    library_strategy = "RNA-Seq",
    library_source = "TRANSCRIPTOMIC",
    library_selection = "cDNA_oligo_dT",
    library_layout = "single",
    platform = "ILLUMINA",
    instrument_model = "NextSeq 550",
    design_description = design_description_text,
    filetype = "fastq",
    filename = paste0(library_ID, "_1.fastq.gz"),
    well = str_split_fixed(sample_name, fixed("_"), 2) %>%
      magrittr::extract(, 2)
  ) %>%
  left_join(
    well_map, by = "well"
  ) %>%
  select(-well)

pwalk(
  dge_sra %>%
    filter(str_starts(sample_name, fixed("dge2"))),
  function(barcode, filename, ...) {
    old_fn <- paste0("asdf", barcode, "_R1.fq")
    if (file.exists(old_fn))
      file.rename(
        old_fn,
        str_replace(filename, fixed(".gz"), "")
      )
  }
)

write_tsv(
  dge_sra,
  here("results", "dge_meta_sra.tsv")
)

# Prepare geo submission

# pull SRA metadata

sra_accession <- read_tsv(here("results", "sra_accession_metadata.tsv"))

dge_geo_samples <- dge_geo %>%
  left_join(
    sra_accession %>%
      select(
        sample_name,
        experiment = accession,
        BioSample = biosample_accession
      ),
    by = "sample_name"
  ) %>%
  left_join(
    dge_sra %>%
      select(sample_name, title, barcode),
    by = "sample_name"
  ) %>%
  transmute(
    sample_name,
    title,
    `source name` = cell_type,
    organism,
    `characteristics: cell_type` = cell_type,
    `characteristics: treatment` = treatment,
    `characteristics: cell_line` = cell_line,
    `characteristics: drug_treatment` = drug_treatment,
    `characteristics: drug_concentration` = drug_concentration,
    `characteristics: stimulant` = stimulant,
    molecule = "polyA RNA",
    description = paste(
      "barcode:", barcode, "batch:", batch
    ),
    experiment,
    BioSample
  )



geo_counts_raw <- tribble(
  ~matrix, ~colnames, ~rownames, ~batch,
  "syn18143723", "syn18143725", "syn18143724", "dge1",
  "syn17103650", "syn17103652", "syn17103651", "dge2",
) %>%
  mutate(
    across(-batch, .fns = syn),
    matrix = map(
      matrix,
      ~read_delim(.x, delim = " ", skip = 3, col_names = c("row_id", "col_id", "count"))
    ),
    colnames = map(
      colnames,
      ~read_csv(.x, col_names = "barcode") %>%
        mutate(col_id = seq_len(n()))
    ),
    rownames = map(
      rownames,
      ~read_csv(.x, col_names = "gene_id") %>%
        mutate(row_id = seq_len(n()))
    )
  ) %>%
  mutate(
    sparse_matrix = pmap(
      .,
      function(matrix, colnames, rownames, ...) {
        matrix %>%
          left_join(colnames, by = "col_id") %>%
          left_join(rownames, by = "row_id") %>%
          left_join(well_map, by = "barcode")
      }
    )
  ) %>%
  select(batch, sparse_matrix) %>%
  unnest(sparse_matrix) %>%
  mutate(
    sample_name = paste(batch, well, sep = "_")
  )

geo_counts_long <- geo_counts_raw %>%
  select(sample_name, gene_id, count) %>%
  arrange(sample_name, gene_id)

write_tsv(
  geo_counts_long,
  here("results", "count_matrix_long.csv.gz")
)

geo_counts_wide <- geo_counts_long %>%
  pivot_wider(id_cols = gene_id, names_from = sample_name, values_from = count)

write_tsv(
  geo_counts_wide,
  here("results", "count_matrix.csv.gz")
)
