suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

out_dir <- file.path("Results", "01_Preprocessing_maxEE88")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_rank <- function(x) {
  x <- str_replace(as.character(x), "^[dkpcofgs]__", "")
  x[is.na(x) | x %in% c("", "__", "NA", "uncultured", "Unassigned")] <- ""
  x
}

read_raw_taxonomy <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE)
  split <- str_split_fixed(x$Taxon, ";", 8)[, 1:7, drop = FALSE]
  colnames(split) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax <- as.data.frame(split, stringsAsFactors = FALSE)
  tax[] <- lapply(tax, clean_rank)
  tax$FeatureID <- x$`Feature ID`
  tax$Raw_Taxonomy <- x$Taxon
  tax |>
    mutate(
      Taxon_Group = case_when(
        Kingdom == "Fungi" & Phylum == "" ~ "Fungi_k",
        Kingdom == "Fungi" & Phylum == "Ascomycota" & Class == "" ~ "Ascomycota_p",
        Kingdom == "Fungi" ~ "Other_fungal_ASV",
        TRUE ~ "Non_fungal"
      )
    )
}

sequence_entropy <- function(sequence) {
  chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  p <- table(chars) / length(chars)
  -sum(p * log2(p)) / 2
}

analyse_method <- function(method, ps_path, repseq_path, taxonomy_path) {
  ps <- readRDS(ps_path)
  seqs <- getElement(read_qza(repseq_path), "data")
  tax <- read_raw_taxonomy(taxonomy_path)
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  ids <- intersect(intersect(rownames(otu), names(seqs)), tax$FeatureID)
  otu <- otu[ids, , drop = FALSE]
  tax <- tax[match(ids, tax$FeatureID), , drop = FALSE]
  seq_text <- as.character(seqs[ids])
  nextera_fwd_pos <- str_locate(seq_text, fixed("CTGTCTCTTATA"))[, "start"]
  nextera_rev_pos <- str_locate(seq_text, fixed("AGATGTGTATAA"))[, "start"]
  nextera_12_pos <- ifelse(is.na(nextera_fwd_pos), nextera_rev_pos, nextera_fwd_pos)

  asv <- tibble(
    Method = method,
    FeatureID = ids,
    Sequence = seq_text,
    Taxon_Group = tax$Taxon_Group,
    Raw_Taxonomy = tax$Raw_Taxonomy,
    Length_bp = width(seqs[ids]),
    GC_percent = as.numeric(letterFrequency(seqs[ids], letters = "GC", as.prob = TRUE)) * 100,
    Normalized_Base_Entropy = vapply(seq_text, sequence_entropy, numeric(1)),
    Nextera_full_exact = str_detect(seq_text, fixed("CTGTCTCTTATACACATCT")) |
      str_detect(seq_text, fixed("AGATGTGTATAAGAGACAG")),
    Nextera_12mer_exact = str_detect(seq_text, fixed("CTGTCTCTTATA")) |
      str_detect(seq_text, fixed("AGATGTGTATAA")),
    Nextera_12mer_start = nextera_12_pos,
    Nextera_12mer_near_end = !is.na(nextera_12_pos) &
      (nextera_12_pos <= 5 | nextera_12_pos + 11 >= width(seqs[ids]) - 4),
    Total_Reads = as.numeric(rowSums(otu)),
    Prevalence = as.integer(rowSums(otu > 0)),
    Max_Sample_Relative = apply(sweep(otu, 2, colSums(otu), "/"), 1, max, na.rm = TRUE)
  )

  rel <- sweep(otu, 2, colSums(otu), "/")
  groups <- c("Fungi_k", "Ascomycota_p")
  sample_group <- lapply(groups, function(g) {
    ix <- tax$Taxon_Group == g
    value <- if (any(ix)) colSums(rel[ix, , drop = FALSE]) else rep(0, ncol(rel))
    tibble(Method = method, SampleID = colnames(rel), Taxon_Group = g, Relative_Abundance = value)
  }) |>
    bind_rows()

  list(asv = asv, sample_group = sample_group)
}

cutadapt <- analyse_method(
  "Cutadapt",
  "Phyloseq/phy_ITS2_cutadapt_q20_ee88_20261015.rds",
  "benchmark_results/q20_ee88_both_full/dada2/cutadapt/paired/repseq.qza",
  "benchmark_results/q20_ee88_both_full/taxonomy/cutadapt/paired/unite_v10_99_full/taxonomy.tsv"
)
itsxpress <- analyse_method(
  "ITSxpress",
  "Phyloseq/phy_ITS2_itsxpress_q20_ee88_20261015.rds",
  "benchmark_results/q20_ee88_both_full/dada2/itsxpress/paired/repseq.qza",
  "benchmark_results/q20_ee88_both_full/taxonomy/itsxpress/paired/unite_v10_99_full/taxonomy.tsv"
)

asv_metrics <- bind_rows(cutadapt$asv, itsxpress$asv)
sample_group <- bind_rows(cutadapt$sample_group, itsxpress$sample_group)

length_summary <- asv_metrics |>
  filter(Taxon_Group != "Non_fungal") |>
  group_by(Method, Taxon_Group) |>
  summarise(
    N_ASV = n(),
    Median_Length = median(Length_bp),
    Q1_Length = quantile(Length_bp, 0.25),
    Q3_Length = quantile(Length_bp, 0.75),
    Min_Length = min(Length_bp),
    Max_Length = max(Length_bp),
    Median_GC = median(GC_percent),
    Median_Entropy = median(Normalized_Base_Entropy),
    Nextera_full_ASV = sum(Nextera_full_exact),
    Nextera_12mer_ASV = sum(Nextera_12mer_exact),
    Total_Reads = sum(Total_Reads),
    .groups = "drop"
  )

length_tests <- expand_grid(
  Method = c("Cutadapt", "ITSxpress"),
  Taxon_Group = c("Fungi_k", "Ascomycota_p")
) |>
  rowwise() |>
  mutate(
    N_Group = sum(asv_metrics$Method == Method & asv_metrics$Taxon_Group == Taxon_Group),
    P_Wilcoxon_vs_other_fungal = {
      a <- asv_metrics$Length_bp[asv_metrics$Method == Method & asv_metrics$Taxon_Group == Taxon_Group]
      b <- asv_metrics$Length_bp[asv_metrics$Method == Method & asv_metrics$Taxon_Group == "Other_fungal_ASV"]
      if (length(a) > 0 && length(b) > 0) wilcox.test(a, b, exact = FALSE)$p.value else NA_real_
    }
  ) |>
  ungroup() |>
  mutate(P_adjust_BH = p.adjust(P_Wilcoxon_vs_other_fungal, method = "BH"))

cleanup <- read_tsv(
  "benchmark_results/q20_ee88_both_full/cleanup_qc/cleanup_qc_summary.tsv",
  show_col_types = FALSE
) |>
  transmute(
    SampleID = sample_id,
    Nextera_R1_percent = 100 * nextera_r1_matches / input_pairs,
    Nextera_R2_percent = 100 * nextera_r2_matches / input_pairs,
    Nextera_mean_percent = 100 * (nextera_r1_matches + nextera_r2_matches) / (2 * input_pairs),
    Retained_percent = retained_percent
  )

read_fastqc_nextera <- function(path, prefix) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[str_detect(lines, fixed("nextera_transposase_sequence"))]
  parts <- str_split_fixed(lines, "\\t", 2)
  tuples <- str_extract_all(parts[, 2], "\\([0-9]+, [^)]+\\)")
  final_percent <- vapply(tuples, function(x) {
    last <- x[length(x)]
    as.numeric(str_match(last, "\\([0-9]+, ([^)]+)\\)")[, 2])
  }, numeric(1))
  tibble(
    SampleID = str_replace(parts[, 1], "[._]R[12].*$", ""),
    Read = paste0("R", str_match(parts[, 1], "[._]R([12])")[, 2]),
    Percent = final_percent
  ) |>
    pivot_wider(names_from = Read, values_from = Percent) |>
    transmute(
      SampleID,
      !!paste0(prefix, "_R1_percent") := R1,
      !!paste0(prefix, "_R2_percent") := R2,
      !!paste0(prefix, "_mean_percent") := rowMeans(cbind(R1, R2), na.rm = TRUE)
    )
}

raw_fastqc <- read_fastqc_nextera(
  "benchmark_results/q20_ee88_both_full/raw_qc/multiqc/multiqc_data/fastqc_adapter_content_plot.txt",
  "Raw_FastQC_Nextera"
)
post_fastqc <- read_fastqc_nextera(
  "benchmark_results/q20_ee88_both_full/post_cleanup_qc/multiqc/multiqc_data/fastqc_adapter_content_plot.txt",
  "Postcleanup_FastQC_Nextera"
)

metadata <- data.frame(sample_data(readRDS("Phyloseq/phy_ITS2_cutadapt_q20_ee88_20261015.rds")), check.names = FALSE) |>
  tibble::rownames_to_column("SampleID") |>
  select(SampleID, Sample_Type, Subject, Control_Type, PCR_Cycle_Group)

sample_diagnostic <- sample_group |>
  left_join(cleanup, by = "SampleID") |>
  left_join(raw_fastqc, by = "SampleID") |>
  left_join(post_fastqc, by = "SampleID") |>
  left_join(metadata, by = "SampleID")

correlation_table <- sample_diagnostic |>
  mutate(Stratum = "All samples") |>
  bind_rows(
    sample_diagnostic |>
      filter(as.character(Sample_Type) == "True") |>
      mutate(Stratum = "True only"),
    sample_diagnostic |>
      filter(as.character(Sample_Type) == "Control") |>
      mutate(Stratum = "Control only")
  ) |>
  group_by(Method, Taxon_Group, Stratum) |>
  summarise(
    N = sum(complete.cases(Relative_Abundance, Nextera_mean_percent)),
    Spearman_rho = suppressWarnings(cor(Relative_Abundance, Nextera_mean_percent, method = "spearman", use = "complete.obs")),
    P_value = suppressWarnings(cor.test(Relative_Abundance, Nextera_mean_percent, method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) |>
  mutate(P_adjust_BH = p.adjust(P_value, method = "BH"))

fastqc_long <- sample_diagnostic |>
  select(
    Method, SampleID, Taxon_Group, Relative_Abundance, Sample_Type,
    Cleanup_match_percent = Nextera_mean_percent,
    Raw_FastQC_percent = Raw_FastQC_Nextera_mean_percent,
    Postcleanup_FastQC_percent = Postcleanup_FastQC_Nextera_mean_percent
  ) |>
  pivot_longer(
    ends_with("percent"),
    names_to = "Nextera_Metric",
    values_to = "Nextera_percent"
  )

fastqc_correlation_table <- fastqc_long |>
  mutate(Stratum = "All samples") |>
  bind_rows(
    fastqc_long |>
      filter(as.character(Sample_Type) == "True") |>
      mutate(Stratum = "True only")
  ) |>
  group_by(Method, Taxon_Group, Nextera_Metric, Stratum) |>
  summarise(
    N = sum(complete.cases(Relative_Abundance, Nextera_percent)),
    Spearman_rho = suppressWarnings(cor(Relative_Abundance, Nextera_percent, method = "spearman", use = "complete.obs")),
    P_value = suppressWarnings(cor.test(Relative_Abundance, Nextera_percent, method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) |>
  mutate(P_adjust_BH = p.adjust(P_value, method = "BH"))

top_asvs <- asv_metrics |>
  filter(Taxon_Group %in% c("Fungi_k", "Ascomycota_p")) |>
  group_by(Method, Taxon_Group) |>
  slice_max(Total_Reads, n = 15, with_ties = FALSE) |>
  ungroup() |>
  arrange(Method, Taxon_Group, desc(Total_Reads))

cutadapt_top_unresolved <- asv_metrics |>
  filter(Method == "Cutadapt", Taxon_Group %in% c("Fungi_k", "Ascomycota_p"), Length_bp >= 80) |>
  group_by(Taxon_Group) |>
  slice_max(Total_Reads, n = 30, with_ties = FALSE) |>
  ungroup()

itsxpress_candidates <- asv_metrics |>
  filter(Method == "ITSxpress", Length_bp >= 80)

cross_method_matches <- lapply(seq_len(nrow(cutadapt_top_unresolved)), function(i) {
  query <- cutadapt_top_unresolved[i, ]
  hit <- vapply(
    itsxpress_candidates$Sequence,
    function(target) grepl(target, query$Sequence, fixed = TRUE) ||
      grepl(query$Sequence, target, fixed = TRUE),
    logical(1)
  )
  if (!any(hit)) return(NULL)
  itsxpress_candidates[hit, ] |>
    transmute(
      Cutadapt_FeatureID = query$FeatureID,
      Cutadapt_Group = query$Taxon_Group,
      Cutadapt_Length = query$Length_bp,
      Cutadapt_Total_Reads = query$Total_Reads,
      ITSxpress_FeatureID = FeatureID,
      ITSxpress_Group = Taxon_Group,
      ITSxpress_Length = Length_bp,
      ITSxpress_Total_Reads = Total_Reads,
      ITSxpress_Raw_Taxonomy = Raw_Taxonomy
    )
}) |>
  bind_rows()

method_sample_n <- sample_group |>
  distinct(Method, SampleID) |>
  count(Method, name = "N_Samples")

sequence_frequency <- asv_metrics |>
  left_join(method_sample_n, by = "Method") |>
  group_by(Method) |>
  mutate(
    Relative_Total_Frequency = Total_Reads / sum(Total_Reads),
    Prevalence_Fraction = Prevalence / N_Samples,
    Abundance_Rank = rank(-Total_Reads, ties.method = "first")
  ) |>
  ungroup()

cut_all <- sequence_frequency |> filter(Method == "Cutadapt")
its_all <- sequence_frequency |> filter(Method == "ITSxpress")

all_containment_matches <- lapply(seq_len(nrow(its_all)), function(i) {
  target <- its_all[i, ]
  if (target$Length_bp < 100) return(NULL)
  hit_type <- vapply(cut_all$Sequence, function(query) {
    if (identical(query, target$Sequence)) return("Exact_identical")
    if (grepl(target$Sequence, query, fixed = TRUE)) return("ITSxpress_within_Cutadapt")
    if (nchar(query) >= 100 && grepl(query, target$Sequence, fixed = TRUE)) return("Cutadapt_within_ITSxpress")
    ""
  }, character(1))
  hit <- hit_type != ""
  if (!any(hit)) return(NULL)
  cut_all[hit, ] |>
    transmute(
      Cutadapt_FeatureID = FeatureID,
      Cutadapt_Group = Taxon_Group,
      Cutadapt_Length = Length_bp,
      Cutadapt_Total_Reads = Total_Reads,
      Cutadapt_Relative_Total_Frequency = Relative_Total_Frequency,
      Cutadapt_Prevalence = Prevalence,
      ITSxpress_FeatureID = target$FeatureID,
      ITSxpress_Group = target$Taxon_Group,
      ITSxpress_Length = target$Length_bp,
      ITSxpress_Total_Reads = target$Total_Reads,
      ITSxpress_Relative_Total_Frequency = target$Relative_Total_Frequency,
      ITSxpress_Prevalence = target$Prevalence,
      Match_Type = hit_type[hit],
      Match_Length = pmin(Length_bp, target$Length_bp)
    )
}) |>
  bind_rows()

best_containment_matches <- all_containment_matches |>
  arrange(
    ITSxpress_FeatureID,
    desc(Match_Length),
    desc(Cutadapt_Total_Reads)
  ) |>
  group_by(ITSxpress_FeatureID) |>
  slice_head(n = 1) |>
  ungroup()

sequence_match_summary <- tibble(
  Metric = c(
    "Cutadapt ASVs", "ITSxpress ASVs", "Exact-identical ASV pairs",
    "ITSxpress ASVs with >=100-bp exact containment match",
    "Matched ITSxpress ASV read fraction",
    "Matched ITSxpress ASV frequency vs Cutadapt frequency Spearman rho"
  ),
  Value = c(
    nrow(cut_all), nrow(its_all),
    sum(all_containment_matches$Match_Type == "Exact_identical"),
    n_distinct(best_containment_matches$ITSxpress_FeatureID),
    sum(best_containment_matches$ITSxpress_Total_Reads) / sum(its_all$Total_Reads),
    suppressWarnings(cor(
      best_containment_matches$Cutadapt_Relative_Total_Frequency,
      best_containment_matches$ITSxpress_Relative_Total_Frequency,
      method = "spearman"
    ))
  )
)

write_csv(asv_metrics, file.path(out_dir, "05a_unresolved_taxa_ASV_metrics.csv"))
write_csv(length_summary, file.path(out_dir, "05b_unresolved_taxa_length_summary.csv"))
write_csv(length_tests, file.path(out_dir, "05c_unresolved_taxa_length_tests.csv"))
write_csv(sample_diagnostic, file.path(out_dir, "05d_unresolved_taxa_nextera_by_sample.csv"))
write_csv(correlation_table, file.path(out_dir, "05e_unresolved_taxa_nextera_correlations.csv"))
write_csv(fastqc_correlation_table, file.path(out_dir, "05k_unresolved_taxa_all_nextera_metrics_correlations.csv"))
write_csv(top_asvs, file.path(out_dir, "05f_unresolved_taxa_top_ASVs.csv"))
write_csv(cross_method_matches, file.path(out_dir, "05j_cutadapt_unresolved_exact_ITSxpress_matches.csv"))
write_csv(sequence_frequency, file.path(out_dir, "05l_all_ASV_sequence_frequency.csv"))
write_csv(all_containment_matches, file.path(out_dir, "05m_all_exact_containment_pairs.csv"))
write_csv(best_containment_matches, file.path(out_dir, "05n_best_exact_containment_pairs.csv"))
write_csv(sequence_match_summary, file.path(out_dir, "05o_ASV_sequence_match_summary.csv"))
write_csv(
  asv_metrics |>
    filter(Taxon_Group %in% c("Fungi_k", "Ascomycota_p")),
  file.path(out_dir, "05i_unresolved_taxa_all_ASVs_with_sequences.csv")
)

p_length <- asv_metrics |>
  filter(Taxon_Group != "Non_fungal") |>
  ggplot(aes(Taxon_Group, Length_bp, fill = Taxon_Group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.14, outlier.size = 0.7, colour = "grey20") +
  facet_wrap(~Method, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  labs(x = NULL, y = "ASV length (bp)", title = "Sequence length of unresolved fungal ASVs") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(out_dir, "05g_unresolved_taxa_ASV_length.png"), p_length,
       width = 10, height = 5.5, dpi = 200, bg = "white")

p_nextera <- sample_diagnostic |>
  ggplot(aes(Nextera_mean_percent, Relative_Abundance, colour = Taxon_Group, shape = Sample_Type)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_grid(Taxon_Group ~ Method, scales = "free_y") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  labs(
    x = "Nextera matches among input reads (mean of R1/R2, %)",
    y = "Relative abundance",
    title = "Unresolved fungal taxa vs Nextera adapter-match rate"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "05h_unresolved_taxa_vs_nextera.png"), p_nextera,
       width = 11, height = 7.5, dpi = 200, bg = "white")

p_sequence_length <- sequence_frequency |>
  ggplot(aes(Length_bp, colour = Method, fill = Method)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  scale_colour_manual(values = c("Cutadapt" = "#0072B2", "ITSxpress" = "#D55E00")) +
  scale_fill_manual(values = c("Cutadapt" = "#0072B2", "ITSxpress" = "#D55E00")) +
  labs(x = "ASV sequence length (bp)", y = "Density", title = "ASV sequence-length frequency") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "05p_ASV_sequence_length_frequency.png"), p_sequence_length,
       width = 9, height = 5.5, dpi = 200, bg = "white")

p_rank_abundance <- sequence_frequency |>
  mutate(Display_Group = if_else(Taxon_Group %in% c("Fungi_k", "Ascomycota_p"), Taxon_Group, "Other fungal ASV")) |>
  ggplot(aes(Abundance_Rank, Relative_Total_Frequency, colour = Display_Group)) +
  geom_point(alpha = 0.75, size = 1.3) +
  facet_wrap(~Method) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_percent()) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "ASV abundance rank", y = "Fraction of all fungal reads", colour = NULL,
       title = "ASV sequence-frequency rank") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "05q_ASV_rank_abundance.png"), p_rank_abundance,
       width = 10, height = 5.8, dpi = 200, bg = "white")

p_match_frequency <- best_containment_matches |>
  ggplot(aes(
    Cutadapt_Relative_Total_Frequency,
    ITSxpress_Relative_Total_Frequency,
    colour = Match_Type
  )) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_x_log10(labels = scales::label_percent()) +
  scale_y_log10(labels = scales::label_percent()) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    x = "Cutadapt ASV frequency",
    y = "ITSxpress ASV frequency",
    colour = NULL,
    title = ">=100-bp exact sequence-containment matches"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "05r_matched_ASV_sequence_frequency.png"), p_match_frequency,
       width = 8, height = 7, dpi = 200, bg = "white")

cat("\nLENGTH SUMMARY\n")
print(length_summary, n = Inf)
cat("\nLENGTH TESTS\n")
print(length_tests, n = Inf)
cat("\nNEXTERA CORRELATIONS\n")
print(correlation_table, n = Inf)
cat("\nSEQUENCE MATCH SUMMARY\n")
print(sequence_match_summary, n = Inf)
