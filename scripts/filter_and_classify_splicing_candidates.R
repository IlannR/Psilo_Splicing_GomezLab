#### Introduction ####

# filter_and_classify_splicing_candidates.R
# Ilann Rouillé
# 04/02/2026
# RStudio Version 2024.12.0+467 

#### library ####

library(dplyr) 
library(stringr) 
library(purrr) 
library(biomaRt)

#### Data Importation ####

# Import alternative splicing data (skipped-exon events).
# The input file is not included in this repository.
# To run this script, place allcelltypes_SE.csv in the local data/ folder.
input_file <- "data/allcelltypes_SE.csv"

if (!file.exists(input_file)) {
  stop("Input file not found. Please place allcelltypes_SE.csv in the data/ folder before running this script.")
}

dataset_SkippedExon <- read.csv(input_file, header = TRUE)

#### Data filtering ####

PV_psilo_SE <- filter(dataset_SkippedExon, Cell_type =="PV", Treatment == "psilo", abs(IncLevelDifference) >= 0.1 ) 

No_01IncLevel <- function(x) { 
  str_split(as.character(x), ",\\s*") |> 
    map_lgl(\(.s) { 
      v <- suppressWarnings(as.numeric(.s))
      !anyNA(v) && all(v > 0 & v < 1) 
    })
}

Filtered_SE <- PV_psilo_SE %>%
  filter(No_01IncLevel(IncLevel1) & No_01IncLevel(IncLevel2)) 

Filtered_genes <-Filtered_SE$GeneID

##### Annotation of genes ####

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(stringr)
})

#############################
## 0) INPUT
#############################

# Put your gene list here (character vector)

genes <- Filtered_genes

# Clean input

genes <- unique(genes)
genes <- genes[!is.na(genes) & genes != ""]

#############################
## 1) BioMart query (mouse)
#############################

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

annot <- getBM(
  attributes = c(
    "external_gene_name",
    "description",
    "gene_biotype",
    "go_id",
    "name_1006",
    "namespace_1003",
    "interpro_short_description",
    "pfam"
  ),
  filters = "external_gene_name",
  values  = genes,
  mart    = mart
)

cat("Input genes:", length(genes), "\n")
cat("Rows returned:", nrow(annot), "\n")
cat("Genes matched:", length(unique(annot$external_gene_name)), "\n\n")

#############################
## 2) Collapse per gene
#############################

gene_tbl <- annot %>%
  filter(!is.na(external_gene_name)) %>%
  group_by(external_gene_name) %>%
  summarise(
    description  = first(description),
    gene_biotype = first(gene_biotype),
    go_terms     = paste(unique(na.omit(name_1006)), collapse = " | "),
    go_ns        = paste(unique(na.omit(namespace_1003)), collapse = " | "),
    interpro     = paste(unique(na.omit(interpro_short_description)), collapse = " | "),
    pfam_ids     = paste(unique(na.omit(pfam)), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    go_terms = ifelse(is.na(go_terms), "", go_terms),
    interpro = ifelse(is.na(interpro), "", interpro),
    pfam_ids = ifelse(is.na(pfam_ids), "", pfam_ids)
  )

#############################
## 3) Robust utility functions
#############################

split_pipe <- function(text) {
  if (length(text) > 1) text <- paste(text, collapse = " | ")
  if (is.na(text) || text == "") return(character(0))
  str_split(text, " \\| ", simplify = FALSE)[[1]]
}

count_hits <- function(text, kws) {
  parts <- split_pipe(text)
  if (length(parts) == 0) return(0L)
  sum(str_detect(parts, regex(paste(kws, collapse="|"), ignore_case=TRUE)))
}

get_evidence <- function(text, kws, n = 3) {
  parts <- split_pipe(text)
  if (length(parts) == 0) return(NA_character_)
  hits <- str_subset(parts, regex(paste(kws, collapse="|"), ignore_case=TRUE))
  hits <- unique(str_trim(hits))
  if (length(hits) == 0) return(NA_character_)
  paste(head(hits, n), collapse="; ")
}

#############################
## 4) Keyword dictionaries
#############################

# Splicing / RNA processing (strict, GO-based)
kw_splicing <- c(
  "splice", "splicing", "spliceosome", "SRE", "ESE", "ISE", "ISS",
  "RNA binding", "mRNA processing", "RNA processing", "isoform","exon","intron",
  "pre-mRNA", "ribonucleoprotein", "nuclear speck", "posttranscrip","RRM",
  "polyadenyl", "RNA editing","RNA-binding", "post-transcrip", "spliceosomal","RBP"
)

kw_splicing_domains <- c(
  "RRM", "RBM", "RBP","RNA recognition motif", "KH domain",
  "DEAD", "DEAH", "helicase",
  "ribonucleoprotein", "snRNP", "spliceosomal"
)

# Transcription / chromatin / epigenetic (GO + InterPro support)
kw_tx_go <- c(
  "transcription", "DNA-binding",
  "chromatin", "histone", "nucleosome", "epigen", "transcription regulator",
  "methyl", "acetyl", "deacetyl", "remodel",
  "coactivator", "corepressor", "repressor", "cis‐regulatory", "enhancers", "silencers"
)

# InterPro motifs that support chromatin/epigenetic regulation (not splicing!)
kw_tx_domains <- c(
  "Bromodomain", "PHD", "Tudor", "MBT",
  "WD40", "WD repeat",
  "SET domain", "SWI", "SNF",
  "Histone", "Chromatin", "WAC"
)

# Ion channel (GO-based; very specific terms)
kw_ion <- c(
  "ion channel", "channel activity", "voltage-gated",
  "potassium channel", "sodium channel", "calcium channel",
  "chloride channel", "cation channel", "conductance"
)

kw_ion_domains <- c(
  "ion channel", "voltage-gated", "potassium channel",
  "sodium channel", "calcium channel", "chloride channel",
  "transmembrane"
)

# Synapse / neurotransmission (GO-based)
kw_synapse <- c(
  "synapse", "synaptic", "presynap", "postsynap", "synap",
  "neurotransmitter receptor", "receptor activity", "neurotransmitter",
  "active zone", "vesicle", "exocytosis", "endocytosis",
  "glutamatergic", "GABA", "AMPA", "NMDA", "kainate", "serotonin"
)

# Neuroplasticity / dendrite / cytoskeleton (GO-based; broad)
kw_plasticity <- c(
  "dendrite", "dendritic", "axon", "neurite",
  "cytoskeleton", "actin", "microtubule",
  "spine", "synapse organization", "axon guidance",
  "cell projection", "neuronal projection", "morphogenesis" 
)

# Signaling pathways (GO-based; broad)
kw_signaling <- c(
  "kinase", "phosphorylation", "GTPase", "Ras", "Rho",
  "second messenger", "cAMP", "MAPK", "PI3K", "mTOR",
  "TGF-beta", "BMP", "Wnt", "Notch",
  "receptor signaling", "signal transduction"
)

kw_signaling_domains <- c(
  "protein kinase", "kinase domain", "GTPase",
  "SH2", "SH3", "PH domain", "Ras", "Rho"
)

#############################
## 5) Count-based scoring + category calls
#############################

# Thresholds (tuneable)
thr_splicing   <- 2  
thr_ion        <- 2 
thr_synapse    <- 2   
thr_plasticity <- 2   
thr_signaling  <- 2   
thr_tx         <- 2   

scored <- gene_tbl %>%
  rowwise() %>%
  mutate(
    # counts
    n_splicing   = count_hits(go_terms, kw_splicing) + count_hits(interpro, kw_splicing_domains),
    n_ion        = count_hits(go_terms, kw_ion) + count_hits(interpro, kw_ion_domains),
    n_synapse    = count_hits(go_terms, kw_synapse),
    n_plasticity = count_hits(go_terms, kw_plasticity),
    n_signaling  = count_hits(go_terms, kw_signaling) + count_hits(interpro, kw_signaling_domains),
    n_tx         = count_hits(go_terms, kw_tx_go) + count_hits(interpro, kw_tx_domains),
    
    # hits using thresholds
    hit_splicing   = n_splicing   >= thr_splicing,
    hit_ion        = n_ion        >= thr_ion,
    hit_synapse    = n_synapse    >= thr_synapse,
    hit_plasticity = n_plasticity >= thr_plasticity,
    hit_signaling  = n_signaling  >= thr_signaling,
    hit_tx         = n_tx         >= thr_tx,
    
    # evidence snippets
    ev_splicing   = ifelse(hit_splicing,   
                           paste(na.omit(c(
                             get_evidence(go_terms, kw_splicing, n=3),
                             get_evidence(interpro, kw_splicing_domains, n=3)
                           )), collapse=" | "),
                           NA),
    ev_ion        = ifelse(hit_ion,
                           paste(na.omit(c(
                             get_evidence(go_terms, kw_ion, n=3),
                             get_evidence(interpro, kw_ion_domains, n=3)
                           )), collapse=" | "),
                           NA),
    ev_synapse    = ifelse(hit_synapse,    get_evidence(go_terms, kw_synapse, n=3), NA),
    ev_plasticity = ifelse(hit_plasticity, get_evidence(go_terms, kw_plasticity, n=3), NA),
    ev_signaling  = ifelse(hit_signaling,
                           paste(na.omit(c(
                             get_evidence(go_terms, kw_signaling, n=3),
                             get_evidence(interpro, kw_signaling_domains, n=3)
                           )), collapse=" | "),
                           NA),
    ev_tx         = ifelse(hit_tx,
                           paste(na.omit(c(
                             get_evidence(go_terms, kw_tx_go, n=3),
                             get_evidence(interpro, kw_tx_domains, n=3)
                           )), collapse=" | "),
                           NA),
    
    # primary category (ONE label)
    primary_category = {
      scores <- c(
        n_splicing,
        n_ion,
        n_synapse,
        n_plasticity,
        n_signaling,
        n_tx
      )
      
      hits <- c(
        hit_splicing,
        hit_ion,
        hit_synapse,
        hit_plasticity,
        hit_signaling,
        hit_tx
      )
      
      labels <- c(
        "Splicing / RNA processing",
        "Ion channel",
        "Synapse / neurotransmission",
        "Neuroplasticity / dendrite / cytoskeleton",
        "Signaling pathways",
        "Transcription / chromatin / epigenetic"
      )
      
      if (!any(hits)) {
        NA_character_
      } else {
        labels[which.max(ifelse(hits, scores, -Inf))]
      }
    },
    
    # Assign the evidence snippets corresponding to the selected primary category      
    primary_evidence = case_when( 
      primary_category == "Splicing / RNA processing" ~ ev_splicing,
      primary_category == "Synapse / neurotransmission" ~ ev_synapse,
      primary_category == "Neuroplasticity / dendrite / cytoskeleton" ~ ev_plasticity,
      primary_category == "Signaling pathways" ~ ev_signaling,
      primary_category == "Transcription / chromatin / epigenetic" ~ ev_tx,
      primary_category == "Ion channel" ~ ev_ion,
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()


#############################
## 6) Output table (drop uncategorized genes)
#############################

out <- scored %>% # Create the final output table from the scored gene table
  filter(!is.na(primary_category)) %>%
  transmute(
    gene = external_gene_name,
    primary_category,
    primary_evidence,
    description,
    # keep counts + flags for QC (you can remove later)
    n_splicing, n_ion, n_synapse, n_plasticity, n_signaling, n_tx,
    hit_splicing, hit_ion, hit_synapse, hit_plasticity, hit_signaling, hit_tx
  ) %>%
  arrange(primary_category, gene) # Arrange the output table by primary category and gene name

# Create output folders if they do not already exist
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/gene_classification")) dir.create("results/gene_classification")

# Export the gene category table
write.csv(
  out,
  "results/gene_classification/gene_category_table.csv",
  row.names = FALSE
)

cat("Categorized genes written:", nrow(out), "\n")
cat("Saved as: gene_category_table.csv in working directory:\n")
cat(getwd(), "\n\n")


#################################################
## Specific_filtering candidates
#################################################

Specific_filtering <- scored %>%
  filter(
    (hit_synapse | hit_plasticity) &
      (hit_signaling | hit_tx)
  ) %>%
  transmute(
    gene = external_gene_name,
    description,
    synapse = hit_synapse,
    plasticity = hit_plasticity,
    ion = hit_ion,
    signaling = hit_signaling,
    transcription = hit_tx
  )

# Export the stricter filtered candidate table
write.csv(
  Specific_filtering,
  "results/gene_classification/specific_filtering_candidates.csv",
  row.names = FALSE
)
cat("Categorized genes written:", nrow(Specific_filtering), "\n")

