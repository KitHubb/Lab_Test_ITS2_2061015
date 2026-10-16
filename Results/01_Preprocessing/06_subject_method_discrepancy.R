suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(vegan)
  library(nlme)
})

out_dir <- file.path("Results", "01_Preprocessing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

get_species_count <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  tax <- data.frame(tax_table(ps), check.names = FALSE)
  species_label <- case_when(
    tax$Species != "" ~ str_replace_all(tax$Species, "_", " "),
    tax$Genus != "" ~ paste0(tax$Genus, " sp."),
    tax$Family != "" ~ paste0(tax$Family, " (unclassified)"),
    TRUE ~ "Unclassified fungi"
  )
  rowsum(otu, group = species_label, reorder = FALSE)
}

align_rows <- function(x, rows, columns) {
  out <- matrix(0, nrow = length(rows), ncol = length(columns),
                dimnames = list(rows, columns))
  out[rownames(x), colnames(x)] <- x
  out
}

clr_columns <- function(x, pseudocount = 0.5) {
  lx <- log(x + pseudocount)
  sweep(lx, 2, colMeans(lx), "-")
}

phy_cutadapt <- readRDS("Phyloseq/phy_ITS2_cutadapt_q20_ee22_20261015.rds")
phy_itsxpress <- readRDS("Phyloseq/phy_ITS2_itsxpress_q20_ee22_20261015.rds")

cut <- get_species_count(phy_cutadapt)
its <- get_species_count(phy_itsxpress)
samples <- intersect(colnames(cut), colnames(its))
taxa <- union(rownames(cut), rownames(its))
cut <- align_rows(cut, taxa, samples)
its <- align_rows(its, taxa, samples)

cut_rel <- sweep(cut, 2, colSums(cut), "/")
its_rel <- sweep(its, 2, colSums(its), "/")
cut_clr <- clr_columns(cut)
its_clr <- clr_columns(its)

metadata <- data.frame(sample_data(phy_cutadapt), check.names = FALSE) |>
  rownames_to_column("SampleID") |>
  filter(SampleID %in% samples)

cleanup <- read_tsv(
  "benchmark_results/q20_ee22_both_full/cleanup_qc/cleanup_qc_summary.tsv",
  show_col_types = FALSE
) |>
  transmute(
    SampleID = sample_id,
    Nextera_percent = 100 * (nextera_r1_matches + nextera_r2_matches) / (2 * input_pairs),
    Cleanup_retained_percent = retained_percent
  )

paired_metrics <- tibble(
  SampleID = samples,
  Aitchison_Method_Distance = vapply(samples, function(s) {
    sqrt(sum((cut_clr[, s] - its_clr[, s])^2))
  }, numeric(1)),
  Bray_Method_Distance = vapply(samples, function(s) {
    as.numeric(vegdist(rbind(cut_rel[, s], its_rel[, s]), method = "bray"))
  }, numeric(1)),
  Species_L1_Difference = vapply(samples, function(s) {
    sum(abs(cut_rel[, s] - its_rel[, s])) / 2
  }, numeric(1)),
  Cutadapt_Reads = colSums(cut),
  ITSxpress_Reads = colSums(its),
  Log2_ITSxpress_Cutadapt_Read_Ratio = log2((colSums(its) + 0.5) / (colSums(cut) + 0.5))
) |>
  left_join(metadata, by = "SampleID") |>
  left_join(cleanup, by = "SampleID") |>
  mutate(
    Subject = factor(Subject),
    PCR_Cycle_Group = factor(PCR_Cycle_Group),
    Skin_Site = factor(Skin_Site),
    Biological_Unit = factor(Biological_Unit)
  )

true_metrics <- paired_metrics |>
  filter(as.character(Sample_Type) == "True", Subject %in% c("KSY", "PYH")) |>
  droplevels()

tremulum_cutadapt <- if ("Aureobasidium tremulum" %in% rownames(cut_rel)) {
  cut_rel["Aureobasidium tremulum", ]
} else {
  rep(0, ncol(cut_rel))
}
true_metrics <- true_metrics |>
  mutate(
    Cutadapt_Aureobasidium_tremulum = as.numeric(tremulum_cutadapt[SampleID]),
    ITSxpress_Read_Recovery = ITSxpress_Reads / Cutadapt_Reads
  )

tremulum_correlations <- true_metrics |>
  group_by(Subject) |>
  summarise(
    N = n(),
    Tremulum_vs_Aitchison_rho = cor(Cutadapt_Aureobasidium_tremulum, Aitchison_Method_Distance, method = "spearman"),
    Tremulum_vs_Aitchison_P = cor.test(Cutadapt_Aureobasidium_tremulum, Aitchison_Method_Distance, method = "spearman", exact = FALSE)$p.value,
    Tremulum_vs_ReadRecovery_rho = cor(Cutadapt_Aureobasidium_tremulum, ITSxpress_Read_Recovery, method = "spearman"),
    Tremulum_vs_ReadRecovery_P = cor.test(Cutadapt_Aureobasidium_tremulum, ITSxpress_Read_Recovery, method = "spearman", exact = FALSE)$p.value,
    Tremulum_vs_Qubit_rho = cor(Cutadapt_Aureobasidium_tremulum, Qubit_ng_uL, method = "spearman"),
    Tremulum_vs_Qubit_P = cor.test(Cutadapt_Aureobasidium_tremulum, Qubit_ng_uL, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )

subject_summary <- true_metrics |>
  group_by(Subject) |>
  summarise(
    N = n(),
    Qubit_Median = median(Qubit_ng_uL),
    Qubit_IQR = IQR(Qubit_ng_uL),
    Aitchison_Median = median(Aitchison_Method_Distance),
    Aitchison_IQR = IQR(Aitchison_Method_Distance),
    Bray_Median = median(Bray_Method_Distance),
    Bray_IQR = IQR(Bray_Method_Distance),
    Species_L1_Median = median(Species_L1_Difference),
    ITSxpress_Read_Recovery_Median = median(ITSxpress_Reads / Cutadapt_Reads),
    Nextera_Median = median(Nextera_percent),
    .groups = "drop"
  )

subject_tests <- tibble(
  Metric = c(
    "Qubit_ng_uL", "Aitchison_Method_Distance", "Bray_Method_Distance",
    "Species_L1_Difference", "Log2_ITSxpress_Cutadapt_Read_Ratio", "Nextera_percent"
  ),
  P_Wilcoxon_KSY_vs_PYH = c(
    wilcox.test(Qubit_ng_uL ~ Subject, data = true_metrics, exact = FALSE)$p.value,
    wilcox.test(Aitchison_Method_Distance ~ Subject, data = true_metrics, exact = FALSE)$p.value,
    wilcox.test(Bray_Method_Distance ~ Subject, data = true_metrics, exact = FALSE)$p.value,
    wilcox.test(Species_L1_Difference ~ Subject, data = true_metrics, exact = FALSE)$p.value,
    wilcox.test(Log2_ITSxpress_Cutadapt_Read_Ratio ~ Subject, data = true_metrics, exact = FALSE)$p.value,
    wilcox.test(Nextera_percent ~ Subject, data = true_metrics, exact = FALSE)$p.value
  )
) |>
  mutate(P_adjust_BH = p.adjust(P_Wilcoxon_KSY_vs_PYH, method = "BH"))

qubit_correlations <- true_metrics |>
  mutate(Stratum = "All True samples") |>
  bind_rows(
    true_metrics |> filter(Subject == "KSY") |> mutate(Stratum = "KSY"),
    true_metrics |> filter(Subject == "PYH") |> mutate(Stratum = "PYH")
  ) |>
  group_by(Stratum) |>
  summarise(
    N = n(),
    Qubit_vs_Aitchison_rho = cor(Qubit_ng_uL, Aitchison_Method_Distance, method = "spearman"),
    P_value = cor.test(Qubit_ng_uL, Aitchison_Method_Distance,
                       method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) |>
  mutate(P_adjust_BH = p.adjust(P_value, method = "BH"))

cycle_summary <- true_metrics |>
  group_by(Subject, PCR_Cycle_Group) |>
  summarise(
    N = n(),
    Median_Aitchison = median(Aitchison_Method_Distance),
    Median_Qubit = median(Qubit_ng_uL),
    .groups = "drop"
  )

site_summary <- true_metrics |>
  group_by(Subject, Skin_Site) |>
  summarise(
    N = n(),
    Median_Aitchison = median(Aitchison_Method_Distance),
    Median_Qubit = median(Qubit_ng_uL),
    .groups = "drop"
  )

mixed_model <- lme(
  log(Aitchison_Method_Distance) ~ Subject + scale(log10(Qubit_ng_uL)) +
    PCR_Cycle_Group + scale(Nextera_percent),
  random = ~1 | Biological_Unit,
  data = true_metrics,
  method = "REML",
  na.action = na.omit
)

mixed_model_table <- data.frame(summary(mixed_model)$tTable, check.names = FALSE) |>
  rownames_to_column("Term") |>
  as_tibble() |>
  dplyr::rename(
    Estimate = Value,
    Std_Error = `Std.Error`,
    DF = DF,
    T_value = `t-value`,
    P_value = `p-value`
  ) |>
  mutate(P_adjust_BH = p.adjust(P_value, method = "BH"))

variance_components <- VarCorr(mixed_model)
variance_table <- tibble(
  Component = c("Biological unit (subject/site)", "Residual within biological unit"),
  Variance = c(
    as.numeric(variance_components[1, "Variance"]),
    as.numeric(variance_components[2, "Variance"])
  )
) |>
  mutate(Proportion = Variance / sum(Variance))

abs_diff <- abs(its_rel - cut_rel)
taxon_contribution <- as.data.frame(t(abs_diff), check.names = FALSE) |>
  rownames_to_column("SampleID") |>
  pivot_longer(-SampleID, names_to = "Species", values_to = "Absolute_Relative_Difference") |>
  left_join(true_metrics |> select(SampleID, Subject), by = "SampleID") |>
  filter(!is.na(Subject)) |>
  group_by(Subject, Species) |>
  summarise(Mean_Absolute_Difference = mean(Absolute_Relative_Difference), .groups = "drop") |>
  pivot_wider(names_from = Subject, values_from = Mean_Absolute_Difference, values_fill = 0) |>
  mutate(
    PYH_minus_KSY = PYH - KSY,
    Combined = PYH + KSY
  ) |>
  arrange(desc(PYH_minus_KSY))

top_pyh_contributors <- taxon_contribution |>
  slice_head(n = 20)

taxon_method_direction <- bind_rows(
  as.data.frame(t(cut_rel), check.names = FALSE) |>
    rownames_to_column("SampleID") |>
    pivot_longer(-SampleID, names_to = "Species", values_to = "Relative_Abundance") |>
    mutate(Method = "Cutadapt"),
  as.data.frame(t(its_rel), check.names = FALSE) |>
    rownames_to_column("SampleID") |>
    pivot_longer(-SampleID, names_to = "Species", values_to = "Relative_Abundance") |>
    mutate(Method = "ITSxpress")
) |>
  left_join(true_metrics |> select(SampleID, Subject), by = "SampleID") |>
  filter(!is.na(Subject), Species %in% top_pyh_contributors$Species) |>
  group_by(Subject, Method, Species) |>
  summarise(Mean_Relative_Abundance = mean(Relative_Abundance), .groups = "drop") |>
  pivot_wider(names_from = Method, values_from = Mean_Relative_Abundance, values_fill = 0) |>
  mutate(ITSxpress_minus_Cutadapt = ITSxpress - Cutadapt) |>
  arrange(Subject, desc(abs(ITSxpress_minus_Cutadapt)))

# Underlying biological composition: average the two methods, then average cycles
# within each subject/site biological unit before exploratory PERMANOVA.
consensus_clr <- (cut_clr + its_clr) / 2
unit_meta <- true_metrics |>
  distinct(Biological_Unit, Subject, Skin_Site)
unit_clr <- lapply(levels(true_metrics$Biological_Unit), function(unit) {
  ids <- true_metrics$SampleID[true_metrics$Biological_Unit == unit]
  rowMeans(consensus_clr[, ids, drop = FALSE])
}) |>
  do.call(what = cbind)
colnames(unit_clr) <- levels(true_metrics$Biological_Unit)
unit_meta <- unit_meta[match(colnames(unit_clr), unit_meta$Biological_Unit), , drop = FALSE]
unit_distance <- dist(t(unit_clr))
set.seed(42)
biology_permanova <- adonis2(
  unit_distance ~ Subject + Skin_Site,
  data = unit_meta,
  permutations = 999,
  by = "margin"
)
biology_permanova_table <- data.frame(biology_permanova, check.names = FALSE) |>
  rownames_to_column("Term") |>
  as_tibble()

write_csv(paired_metrics, file.path(out_dir, "06a_paired_method_discrepancy_by_sample.csv"))
write_csv(subject_summary, file.path(out_dir, "06b_subject_discrepancy_summary.csv"))
write_csv(subject_tests, file.path(out_dir, "06c_subject_discrepancy_tests.csv"))
write_csv(qubit_correlations, file.path(out_dir, "06d_qubit_discrepancy_correlations.csv"))
write_csv(cycle_summary, file.path(out_dir, "06e_cycle_discrepancy_summary.csv"))
write_csv(site_summary, file.path(out_dir, "06f_site_discrepancy_summary.csv"))
write_csv(mixed_model_table, file.path(out_dir, "06g_adjusted_mixed_model.csv"))
write_csv(variance_table, file.path(out_dir, "06h_mixed_model_variance.csv"))
write_csv(taxon_contribution, file.path(out_dir, "06i_taxa_method_difference_by_subject.csv"))
write_csv(biology_permanova_table, file.path(out_dir, "06j_subject_site_permanova.csv"))
write_csv(taxon_method_direction, file.path(out_dir, "06o_taxa_method_difference_direction.csv"))
write_csv(tremulum_correlations, file.path(out_dir, "06p_tremulum_correlations.csv"))

p_subject <- ggplot(true_metrics, aes(Subject, Aitchison_Method_Distance, fill = Subject)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(aes(colour = PCR_Cycle_Group), width = 0.13, size = 2, alpha = 0.8) +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  scale_colour_brewer(palette = "Paired") +
  labs(x = NULL, y = "Paired Cutadapt–ITSxpress Aitchison distance",
       colour = "PCR cycle", title = "Method discrepancy by subject") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_qubit <- ggplot(true_metrics, aes(Qubit_ng_uL, Aitchison_Method_Distance, colour = Subject)) +
  geom_point(aes(shape = PCR_Cycle_Group), size = 2.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Qubit DNA concentration (ng/uL)",
       y = "Paired Cutadapt–ITSxpress Aitchison distance",
       colour = "Subject", shape = "PCR cycle",
       title = "Method discrepancy vs DNA concentration") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_site <- ggplot(site_summary, aes(Skin_Site, Median_Aitchison, colour = Subject, group = Subject)) +
  geom_point(size = 2.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = NULL, y = "Median paired Aitchison distance",
       colour = "Subject", title = "Method discrepancy across skin sites") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

p_taxa <- top_pyh_contributors |>
  mutate(Species = reorder(Species, PYH_minus_KSY)) |>
  ggplot(aes(PYH_minus_KSY, Species, fill = PYH_minus_KSY > 0)) +
  geom_col() +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_fill_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#1B9E77"), guide = "none") +
  labs(x = "Mean absolute method difference: PYH minus KSY", y = NULL,
       title = "Taxa contributing more strongly to PYH method discrepancy") +
  theme_bw(base_size = 10)

ggsave(file.path(out_dir, "06k_method_discrepancy_by_subject.png"), p_subject,
       width = 8, height = 6, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "06l_method_discrepancy_vs_qubit.png"), p_qubit,
       width = 9, height = 6.5, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "06m_method_discrepancy_by_site.png"), p_site,
       width = 11, height = 6, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "06n_taxa_driving_PYH_method_difference.png"), p_taxa,
       width = 10, height = 7, dpi = 200, bg = "white")

cat("\nSUBJECT SUMMARY\n")
print(subject_summary, n = Inf)
cat("\nSUBJECT TESTS\n")
print(subject_tests, n = Inf)
cat("\nQUBIT CORRELATIONS\n")
print(qubit_correlations, n = Inf)
cat("\nMIXED MODEL\n")
print(mixed_model_table, n = Inf)
cat("\nVARIANCE\n")
print(variance_table, n = Inf)
cat("\nBIOLOGY PERMANOVA\n")
print(biology_permanova_table, n = Inf)
cat("\nTOP TAXA\n")
print(top_pyh_contributors, n = Inf)
cat("\nTAXA DIRECTION\n")
print(taxon_method_direction, n = Inf)
cat("\nTREMULUM CORRELATIONS\n")
print(tremulum_correlations, n = Inf)
