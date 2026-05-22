#### Introduction ####

# annotate_splicing_event_transcripts_and_motifs
# Ilann Rouillé
# 28/10/2025
# RStudio Version 2024.12.0+467 

# Note:
# This script requires internet access because it queries Ensembl BioMart.
# If Ensembl BioMart is temporarily unavailable or slow, the script may fail or take longer to run.
# For reproducible execution, restart R before sourcing the script.

#### Script ####

## ============================================================
## Goal: Automatise alternative splicing events analysis
## ============================================================

# Load required packages while suppressing startup messages.
suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(stringr)
})

# Connect to the Ensembl mouse BioMart dataset.
# If a BioMart object is already provided, the function reuses it instead of opening a new connection.
get_mouse_mart <- function(mart = NULL) {
  if (!is.null(mart)) return(mart)
  useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
}

# Convert RefSeq mm39 chromosome accessions to Ensembl-style chromosome names.
# rMATS outputs may use RefSeq accessions such as NC_000083.7,
# whereas Ensembl BioMart generally uses chromosome names such as "17".
normalize_chr <- function(chr) {
  map <- c(
    "NC_000067.7"="1","NC_000068.8"="2","NC_000069.7"="3","NC_000070.7"="4",
    "NC_000071.7"="5","NC_000072.7"="6","NC_000073.7"="7","NC_000074.7"="8",
    "NC_000075.7"="9","NC_000076.7"="10","NC_000077.7"="11","NC_000078.7"="12",
    "NC_000079.7"="13","NC_000080.7"="14","NC_000081.7"="15","NC_000082.7"="16",
    "NC_000083.7"="17","NC_000084.7"="18","NC_000085.7"="19",
    "NC_000086.8"="X","NC_000087.8"="Y","NC_000088.8"="MT"
  )
  x <- chr %>% as.character() %>% str_remove("^chr") %>% trimws()
  if (x %in% names(map)) map[[x]] else x
}

#========================================================
# 1) Annotate one splicing event across all transcripts
#========================================================
annotate_event_all_transcripts <- function(gene, chr, start0, end1,
                                           strand_event = NULL,
                                           mart = NULL) {
  mart     <- get_mouse_mart(mart)
  chr_norm <- normalize_chr(chr)
  start1   <- as.integer(start0) + 1
  
  # Retrieve all transcripts associated with the input gene
  tx_all <- biomaRt::getBM(
    attributes = c("external_gene_name","ensembl_transcript_id",
                   "transcript_is_canonical","transcript_biotype",
                   "chromosome_name","strand"),
    filters = "external_gene_name",
    values  = gene,
    mart    = mart
  )
  if (nrow(tx_all) == 0) {
    return(tibble(
      gene = gene, chr = chr_norm, start1 = start1, end1 = end1,
      ensembl_transcript_id = NA_character_, is_canonical = FALSE,
      transcript_biotype = NA_character_,
      region_type = "unknown_gene",
      exon_rank = NA_character_, intron_rank = NA_character_
    ))
  }
  
  tx_all <- tx_all %>%
    mutate(
      is_canonical    = transcript_is_canonical == 1,
      strand          = ifelse(strand == 1, "+", "-"),
      chromosome_name = as.character(chromosome_name)
    )
  
  # Filter transcripts by chromosome when the chromosome is compatible with Ensembl annotation
  if (chr_norm %in% unique(tx_all$chromosome_name)) {
    tx_all <- dplyr::filter(tx_all, chromosome_name == chr_norm)
    if (nrow(tx_all) == 0) {
      return(tibble(
        gene = gene, chr = chr_norm, start1 = start1, end1 = end1,
        ensembl_transcript_id = NA_character_, is_canonical = FALSE,
        transcript_biotype = NA_character_,
        region_type = "unknown_gene_or_chr",
        exon_rank = NA_character_, intron_rank = NA_character_
      ))
    }
  }
  
  # Retrieve exon coordinates for all transcripts
  ex_all <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id","ensembl_exon_id",
                   "exon_chrom_start","exon_chrom_end","rank",
                   "strand","chromosome_name"),
    filters = "ensembl_transcript_id",
    values  = unique(tx_all$ensembl_transcript_id),
    mart    = mart
  ) %>%
    mutate(
      strand          = ifelse(strand == 1, "+", "-"),
      chromosome_name = as.character(chromosome_name)
    ) %>%
    arrange(ensembl_transcript_id, rank)
  
  # Convert the candidate event coordinates into a GenomicRanges object
  gr_evt <- GRanges(
    seqnames = chr_norm,
    ranges   = IRanges(start1, as.integer(end1)),
    strand   = if (is.null(strand_event)) "*" else strand_event
  )
  
  # For each transcript, determine whether the event overlaps an exon, an intron, or no annotated region
  out <- lapply(split(ex_all, ex_all$ensembl_transcript_id), function(ex_tx) {
    tx_id <- unique(ex_tx$ensembl_transcript_id)[1]
    
    gr_ex <- GRanges(
      seqnames = ex_tx$chromosome_name,
      ranges   = IRanges(ex_tx$exon_chrom_start, ex_tx$exon_chrom_end),
      strand   = ex_tx$strand
    )
    mcols(gr_ex)$exon_rank <- ex_tx$rank
    
    # Test whether the event overlaps annotated exons
    h_ex <- findOverlaps(gr_evt, gr_ex, ignore.strand = TRUE)
    if (length(h_ex) > 0) {
      exon_ranks <- sort(unique(mcols(gr_ex)$exon_rank[subjectHits(h_ex)]))
      return(tibble(
        ensembl_transcript_id = tx_id,
        region_type  = "exon",
        exon_rank    = paste(exon_ranks, collapse = ","),
        intron_rank  = NA_character_
      ))
    }
    
    # If the event does not overlap an exon, build intronic intervals between adjacent exons
    if (nrow(ex_tx) >= 2) {
      introns_tbl <- tibble(
        chromosome_name = ex_tx$chromosome_name[1],
        strand          = ex_tx$strand[1],
        intron_rank  = 1:(nrow(ex_tx) - 1),
        intron_start = ex_tx$exon_chrom_end[1:(nrow(ex_tx)-1)] + 1,
        intron_end   = ex_tx$exon_chrom_start[2:nrow(ex_tx)] - 1
      ) %>% filter(intron_start <= intron_end)
      
      if (nrow(introns_tbl) > 0) {
        gr_in <- GRanges(
          seqnames = introns_tbl$chromosome_name,
          ranges   = IRanges(introns_tbl$intron_start, introns_tbl$intron_end),
          strand   = introns_tbl$strand
        )
        mcols(gr_in)$intron_rank <- introns_tbl$intron_rank
        
        # Test whether the event overlaps annotated introns
        h_in <- findOverlaps(gr_evt, gr_in, ignore.strand = TRUE)
        if (length(h_in) > 0) {
          intron_ranks <- sort(unique(mcols(gr_in)$intron_rank[subjectHits(h_in)]))
          return(tibble(
            ensembl_transcript_id = tx_id,
            region_type  = "intron",
            exon_rank    = NA_character_,
            intron_rank  = paste(intron_ranks, collapse = ",")
          ))
        }
      }
    }
    
    # If the event overlaps neither exons nor introns, classify it as no_overlap
    tibble(
      ensembl_transcript_id = tx_id,
      region_type  = "no_overlap",
      exon_rank    = NA_character_,
      intron_rank  = NA_character_
    )
  }) %>% bind_rows()
  
  # Add transcript metadata and arrange the output table
  out %>%
    left_join(tx_all %>% dplyr::select(ensembl_transcript_id, is_canonical, transcript_biotype),
              by = "ensembl_transcript_id") %>%
    mutate(gene = gene, chr = chr_norm, start1 = start1, end1 = end1,
           is_canonical = ifelse(is.na(is_canonical), FALSE, is_canonical)) %>%
    relocate(gene, chr, start1, end1, ensembl_transcript_id, is_canonical, transcript_biotype) %>%
    arrange(desc(is_canonical), ensembl_transcript_id)
}

#========================================================
# 2) Add protein motif annotations
#    Sources: InterPro, PFAM, SMART, PROSITE, PRINTS
#    Robust to missing columns / NULL values; types are converted safely
#========================================================

add_protein_motifs <- function(base, mart = NULL, only_exonic = TRUE) {
  stopifnot("ensembl_transcript_id" %in% names(base))
  mart <- get_mouse_mart(mart)
  
  # Retrieve protein motif/domain annotations from BioMart-linked resources
  dom_raw <- biomaRt::getBM(
    attributes = c(
      "ensembl_transcript_id","ensembl_peptide_id",
      "interpro","interpro_description","interpro_start","interpro_end",
      "pfam","pfam_start","pfam_end",
      "smart","smart_start","smart_end",
      "scanprosite","scanprosite_start","scanprosite_end",
      "prints","prints_start","prints_end"
    ),
    filters = "ensembl_transcript_id",
    values  = unique(base$ensembl_transcript_id),
    mart    = mart
  )
  
  # Convert one motif source into a standardized table format
  mk <- function(df, id_col, start_col, end_col, src, name_col = NULL) {
    if (!(id_col %in% names(df))) return(tibble())
    out <- df[, c("ensembl_transcript_id","ensembl_peptide_id", id_col), drop = FALSE]
    names(out)[3] <- "motif_id"
    
   
    # Convert fields to character and treat empty strings or "NULL" as missing values
    out$ensembl_transcript_id <- as.character(out$ensembl_transcript_id)
    out$ensembl_peptide_id    <- as.character(out$ensembl_peptide_id)
    out$motif_id              <- as.character(out$motif_id)
    out$motif_id[out$motif_id %in% c("", "NULL")] <- NA_character_
    
    # Add motif amino-acid coordinates when available
    if (!(start_col %in% names(df))) df[[start_col]] <- NA_real_
    if (!(end_col   %in% names(df))) df[[end_col]]   <- NA_real_
    out$motif_aa_start <- suppressWarnings(as.numeric(df[[start_col]]))
    out$motif_aa_end   <- suppressWarnings(as.numeric(df[[end_col]]))
    
    # Add motif names/descriptions when available
    if (!is.null(name_col) && name_col %in% names(df)) {
      out$motif_name <- as.character(df[[name_col]])
      out$motif_name[out$motif_name %in% c("", "NULL")] <- NA_character_
    } else {
      out$motif_name <- NA_character_
    }
    
    out$motif_source <- src
    out <- dplyr::filter(out, !is.na(motif_id))
    out
  }
  
  # Build one combined motif table from the different motif/domain sources
  dom <- bind_rows(
    mk(dom_raw, "interpro",    "interpro_start",    "interpro_end",    "INTERPRO", "interpro_description"),
    mk(dom_raw, "pfam",        "pfam_start",        "pfam_end",        "PFAM"),
    mk(dom_raw, "smart",       "smart_start",       "smart_end",       "SMART"),
    mk(dom_raw, "scanprosite", "scanprosite_start", "scanprosite_end", "PROSITE"),
    mk(dom_raw, "prints",      "prints_start",      "prints_end",      "PRINTS")
  )
  
  # If no motif is found, return the input table with empty motif columns
  if (!nrow(dom)) {
    return(base %>% mutate(
      protein_motif = NA_character_, motif_id = NA_character_,
      motif_source = NA_character_,  motif_name = NA_character_,
      motif_url = NA_character_
    ))
  }
  
  # Add URLs and a human-readable motif label
  dom <- dom %>%
    mutate(
      motif_source = toupper(motif_source),
      motif_url = case_when(
        motif_source == "INTERPRO" ~ sprintf("https://www.ebi.ac.uk/interpro/entry/interpro/%s", motif_id),
        motif_source == "PFAM"     ~ sprintf("https://pfam.xfam.org/family/%s", motif_id),
        motif_source == "SMART"    ~ sprintf("https://smart.embl.de/smart/do_annotation.pl?DOMAIN=%s", motif_id),
        motif_source == "PROSITE"  ~ sprintf("https://prosite.expasy.org/%s", motif_id),
        motif_source == "PRINTS"   ~ sprintf("https://www.ebi.ac.uk/interpro/entry/prints/%s", motif_id),
        TRUE ~ NA_character_
      ),
      protein_motif = ifelse(!is.na(motif_name) & motif_name != "",
                             paste0(motif_name, " (", motif_source, ":", motif_id, ")"),
                             paste0(motif_source, ":", motif_id))
    ) %>%
    mutate(across(c(motif_id, motif_source, motif_name, protein_motif, motif_url),
                  ~as.character(.)))
  
  # Collapse motif annotations to one row per transcript
  dom_agg <- dom %>%
    group_by(ensembl_transcript_id) %>%
    summarise(
      protein_motif = paste(unique(protein_motif), collapse = "; "),
      motif_id      = paste(unique(motif_id), collapse = "; "),
      motif_source  = paste(unique(motif_source), collapse = "; "),
      motif_name    = paste(unique(na.omit(motif_name)), collapse = "; "),
      motif_url     = paste(unique(na.omit(motif_url)), collapse = "; "),
      .groups = "drop"
    )
  
  # Join the motif annotations back to the transcript-level event table
  out <- base %>% left_join(dom_agg, by = "ensembl_transcript_id")
  
  # When requested, keep protein motif annotations only for events overlapping exons
  if (only_exonic && "region_type" %in% names(out)) {
    out <- out %>% mutate(
      protein_motif = ifelse(region_type == "exon", protein_motif, NA_character_),
      motif_id      = ifelse(region_type == "exon", motif_id,      NA_character_),
      motif_source  = ifelse(region_type == "exon", motif_source,  NA_character_),
      motif_name    = ifelse(region_type == "exon", motif_name,    NA_character_),
      motif_url     = ifelse(region_type == "exon", motif_url,     NA_character_)
    )
  }
  out
}


#### Exemple ========================================================== #####

# Open one BioMart connection and reuse it for the example analysis
ensembl_mouse <- get_mouse_mart()

# Example: annotate the Brd4 candidate event across Ensembl transcripts
tab  <- annotate_event_all_transcripts(
  gene   = "Brd4",
  chr    = "NC_000083.7",
  start0 = 32472014,
  end1   = 32472174,
  mart   = ensembl_mouse
)

# Add protein motif/domain annotations to the Brd4 transcript-level table
tabBrd4 <- add_protein_motifs(tab, mart = ensembl_mouse, only_exonic = TRUE)

# Display the resulting table in RStudio
View(tabBrd4)

# Save the exemple output table as a CSV file.
# This file contains transcript-level annotation of the exemple candidate splicing event,
# including whether the event overlaps an exon, intron, or no annotated region,
# and associated protein motif/domain annotations when available.
if (!dir.exists("results")) dir.create("results")

write.csv(
  tabBrd4,
  "results/Brd4_transcript_event_protein_motif_annotation.csv",
  row.names = FALSE
)
