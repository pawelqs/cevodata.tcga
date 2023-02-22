## code to prepare `tcga` dataset goes here

# Go here -> https://www.cbioportal.org/study/summary?id=pancan_pcawg_2020
# Download the data using the Download Icon next to the dataset name
# Unpack
# Fix the path below

library(tidyverse)
library(cevomod)
set.seed(1234)


## -------------------------------- SNVs ---------------------------------------

# This file was manually created based on the most common classification to impact
# relation in TCGA_BRCA dataset
variant_classification <- openxlsx::read.xlsx("data-raw/Variant_Classification.xlsx") |>
  as_tibble()

mutations <- read_tsv(
  "/mnt/dane/data/cbioportal/pancan_pcawg_2020/data_mutations.txt",
  comment = "#"
)
glimpse(mutations)

snvs_tcga <- mutations %>%
  left_join(variant_classification) %>%
  transmute(
    sample_id = Tumor_Sample_Barcode,
    chrom = Chromosome,
    pos = Start_Position,
    gene_symbol = Hugo_Symbol,
    ref = Reference_Allele,
    alt = Tumor_Seq_Allele2,
    ref_reads = t_ref_count,
    alt_reads = t_alt_count,
    VAF = alt_reads / (alt_reads + ref_reads),
    DP = ref_reads + alt_reads,
    ref_counts_normal = n_ref_count,
    alt_counts_normal = n_alt_count,
    # variant_class = VARIANT_CLASS,
    # variant_type = Variant_Type,
    Variant_Classification,
    # variant_classification,
    impact = IMPACT,
    consequence = Consequence
  )
class(snvs_tcga) <- c("cevo_snvs", class(snvs_tcga))


## --------------------------- Sample metadata --------------------------------

tcga_meta <- mutations |>
  select(
    sample_id = Tumor_Sample_Barcode,
    genome = NCBI_Build
  ) |>
  group_by(sample_id, genome) |>
  summarise(TMB = n())


## --------------------------------- Clinical ---------------------------------

clinical <- read_tsv(
  "/mnt/dane/data/cbioportal/pancan_pcawg_2020/data_clinical_sample.txt",
  comment = "#"
)

colnames(clinical) <- colnames(clinical) |>
  str_to_lower()

glimpse(clinical)

clinical <- clinical |>
  semi_join(tcga_meta) |>
  rename(sample = sample_type) |>
  select(patient_id, sample_id, sample, everything())


## --------------------------------- cevodata ---------------------------------

tcga <- init_cevodata("TCGA PanCancer data") |>
  add_SNV_data(snvs_tcga, name = "TCGA") |>
  add_sample_data(clinical) |>
  add_sample_data(tcga_meta)

usethis::use_data(tcga, overwrite = TRUE)
usethis::use_data(variant_classification, overwrite = TRUE)
