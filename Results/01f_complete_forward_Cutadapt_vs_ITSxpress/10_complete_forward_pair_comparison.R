suppressPackageStartupMessages({
  library(qiime2R); library(phyloseq); library(dplyr); library(tidyr)
  library(tibble); library(readr); library(stringr); library(purrr)
  library(ggplot2); library(vegan); library(cowplot); library(patchwork)
  library(scales); library(RColorBrewer); library(Biostrings); library(ggrepel)
})
set.seed(42)

out_dir <- file.path("Results", "01f_complete_forward_Cutadapt_vs_ITSxpress")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("Phyloseq", showWarnings = FALSE)

method_levels <- c("Cutadapt Single", "ITSxpress Single")
method_colors <- c("Cutadapt Single" = "#1F78B4", "ITSxpress Single" = "#E31A1C")
comparison_spec <- tibble::tribble(
  ~Comparison, ~Method_1, ~Method_2,
  "Complete forward: Cutadapt vs ITSxpress", "Cutadapt Single", "ITSxpress Single"
)

metadata <- data.frame(
  phyloseq::sample_data(readRDS("Phyloseq/phy_ITS2_cutadapt_single_complete_q20_ee2_20261016.rds")),
  check.names = FALSE
) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::mutate(
    Sample_Type = dplyr::case_when(
      stringr::str_detect(Subject, stringr::regex("mock|zymo", ignore_case = TRUE)) ~ "Mock",
      stringr::str_detect(Subject, stringr::regex("^(ASwab|CSwab|ETube)", ignore_case = TRUE)) ~ "Control",
      TRUE ~ "True"
    ),
    Control_Type = dplyr::case_when(
      stringr::str_detect(Subject, stringr::regex("^ASwab", ignore_case = TRUE)) ~ "Air swab",
      stringr::str_detect(Subject, stringr::regex("^CSwab", ignore_case = TRUE)) ~ "Clean bench swab",
      stringr::str_detect(Subject, stringr::regex("^ETube", ignore_case = TRUE)) ~ "Empty tube",
      Sample_Type == "Mock" ~ "Mock", TRUE ~ "Skin sample"
    ),
    Skin_Site = dplyr::case_when(
      Sample_Type == "True" & (is.na(Direction) | Direction %in% c("", "-")) ~ as.character(Region),
      Sample_Type == "True" ~ paste(Region, Direction), TRUE ~ as.character(Control_Type)
    ),
    Shape_Group = dplyr::case_when(
      Sample_Type == "True" ~ as.character(Subject),
      Sample_Type == "Control" ~ as.character(Control_Type), TRUE ~ "Mock"
    )
  )

clean_tax_text <- function(x) {
  x <- stringr::str_replace(as.character(x), "^[dkpcofgs]__", "")
  x[is.na(x) | x %in% c("", "__", "NA", "uncultured", "Unassigned")] <- ""
  x
}

tax_clean <- function(ps) {
  tax <- data.frame(phyloseq::tax_table(ps), check.names = FALSE)
  tax[] <- lapply(tax, clean_tax_text)
  for (i in seq_len(nrow(tax))) {
    if (tax$Species[i] != "") tax$Species[i] <- stringr::str_replace_all(tax$Species[i], " ", "_")
    else if (tax$Phylum[i] == "") tax[i, 2:7] <- paste0(tax$Kingdom[i], "_k")
    else if (tax$Class[i] == "") tax[i, 3:7] <- paste0(tax$Phylum[i], "_p")
    else if (tax$Order[i] == "") tax[i, 4:7] <- paste0(tax$Class[i], "_c")
    else if (tax$Family[i] == "") tax[i, 5:7] <- paste0(tax$Order[i], "_o")
    else if (tax$Genus[i] == "") tax[i, 6:7] <- paste0(tax$Family[i], "_f")
    else tax$Species[i] <- paste0(tax$Genus[i], "_g")
  }
  phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tax)); ps
}

make_phyloseq <- function(method, table_path, taxonomy_path) {
  count <- qiime2R::read_qza(table_path)[["data"]]
  tax_in <- readr::read_tsv(taxonomy_path, show_col_types = FALSE)
  tax <- stringr::str_split_fixed(tax_in$Taxon, ";", 8)[, 1:7, drop = FALSE]
  colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(tax) <- tax_in$`Feature ID`
  missing <- setdiff(rownames(count), rownames(tax))
  if (length(missing)) tax <- rbind(tax, matrix("", length(missing), 7, dimnames = list(missing, colnames(tax))))
  md <- metadata |> dplyr::filter(SampleID %in% colnames(count)) |>
    dplyr::arrange(match(SampleID, colnames(count))) |> tibble::column_to_rownames("SampleID") |>
    data.frame(check.names = FALSE)
  count <- count[, rownames(md), drop = FALSE]; tax <- tax[rownames(count), , drop = FALSE]
  raw <- phyloseq::phyloseq(phyloseq::otu_table(count, taxa_are_rows = TRUE),
                            phyloseq::tax_table(tax), phyloseq::sample_data(md))
  cleaned_tax <- tax_clean(raw)
  fungi <- phyloseq::prune_taxa(as.character(phyloseq::tax_table(cleaned_tax)[, "Kingdom"]) == "Fungi", cleaned_tax)
  final <- phyloseq::prune_samples(phyloseq::sample_sums(fungi) > 0, fungi)
  final <- phyloseq::prune_taxa(phyloseq::taxa_sums(final) > 0, final)
  audit <- tibble::tibble(
    Method = method,
    Stage = c("Imported ASV table", "Taxonomy cleaned", "Kingdom == Fungi retained", "Zero-sum samples/taxa removed"),
    Samples = c(nsamples(raw), nsamples(cleaned_tax), nsamples(fungi), nsamples(final)),
    Taxa = c(ntaxa(raw), ntaxa(cleaned_tax), ntaxa(fungi), ntaxa(final)),
    Total_Reads = c(sum(sample_sums(raw)), sum(sample_sums(cleaned_tax)), sum(sample_sums(fungi)), sum(sample_sums(final)))
  )
  list(raw = raw, cleaned = final, audit = audit)
}

root <- "benchmark_results/q20_ee22_both_full_cutadapt_forward"
paths <- tibble::tribble(
  ~Method, ~table_path, ~taxonomy_path, ~repseq_path,
  "Cutadapt Single", file.path(root, "dada2/cutadapt/single/table.qza"), file.path(root, "taxonomy/cutadapt/single/unite_v10_99_full/taxonomy.tsv"), file.path(root, "dada2/cutadapt/single/repseq.qza"),
  "ITSxpress Single", file.path(root, "dada2/itsxpress/single/table.qza"), file.path(root, "taxonomy/itsxpress/single/unite_v10_99_full/taxonomy.tsv"), file.path(root, "dada2/itsxpress/single/repseq.qza")
)

results <- purrr::pmap(
  paths[, c("Method", "table_path", "taxonomy_path")],
  function(Method, table_path, taxonomy_path) make_phyloseq(Method, table_path, taxonomy_path)
)
names(results) <- paths$Method
ps_list <- purrr::map(results, "cleaned")
audit <- purrr::map_dfr(results, "audit")

rds_names <- c(
  "Cutadapt Single" = "phy_ITS2_cutadapt_single_complete_q20_ee2_20261016.rds",
  "ITSxpress Single" = "phy_ITS2_itsxpress_single_complete_q20_ee2_20261016.rds"
)
purrr::iwalk(ps_list, ~saveRDS(.x, file.path("Phyloseq", rds_names[[.y]])))

get_species_count <- function(ps) {
  otu <- as(phyloseq::otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  tax <- data.frame(tax_table(ps), check.names = FALSE)
  label <- dplyr::case_when(
    tax$Species != "" ~ stringr::str_replace_all(tax$Species, "_", " "),
    tax$Genus != "" ~ paste0(tax$Genus, " sp."),
    tax$Family != "" ~ paste0(tax$Family, " (unclassified)"), TRUE ~ "Unclassified fungi"
  )
  rowsum(otu, group = label, reorder = FALSE)
}
species_list <- purrr::map(ps_list, get_species_count)
common_samples <- Reduce(intersect, purrr::map(species_list, colnames))
all_species <- Reduce(union, purrr::map(species_list, rownames))
align_matrix <- function(x) {
  z <- matrix(0, length(all_species), length(common_samples), dimnames = list(all_species, common_samples))
  z[rownames(x), colnames(x)] <- x; z
}
aligned <- purrr::map(species_list, align_matrix)

rank_summary <- purrr::imap_dfr(ps_list, function(ps, method) {
  tax <- data.frame(tax_table(ps), check.names = FALSE)
  tibble::tibble(Method = method, ASV = ntaxa(ps), Kingdom = n_distinct(tax$Kingdom),
                 Phylum = n_distinct(tax$Phylum), Class = n_distinct(tax$Class),
                 Order = n_distinct(tax$Order), Family = n_distinct(tax$Family),
                 Genus = n_distinct(tax$Genus), Species = n_distinct(tax$Species))
})

sample_qc <- purrr::imap_dfr(ps_list, function(ps, method) {
  otu <- as(phyloseq::otu_table(ps), "matrix"); if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  tibble::tibble(SampleID = sample_names(ps), Method = method,
                 Read_Depth = as.numeric(sample_sums(ps)),
                 Observed_ASV = as.integer(rowSums(otu > 0)))
}) |> dplyr::left_join(metadata, by = "SampleID") |>
  dplyr::mutate(Method = factor(Method, levels = method_levels))

depth_summary <- sample_qc |> dplyr::group_by(Method, Sample_Type) |>
  dplyr::summarise(N = dplyr::n(), Min = min(Read_Depth), Q1 = quantile(Read_Depth, .25),
                   Median = median(Read_Depth), Q3 = quantile(Read_Depth, .75), Max = max(Read_Depth), .groups = "drop")

qc_wide <- sample_qc |> dplyr::select(SampleID, Method, Read_Depth, Observed_ASV) |>
  tidyr::pivot_wider(names_from = Method, values_from = c(Read_Depth, Observed_ASV))
paired_qc_tests <- purrr::map_dfr(c("Read_Depth", "Observed_ASV"), function(metric) {
  x <- qc_wide[[paste0(metric, "_Cutadapt Single")]]
  y <- qc_wide[[paste0(metric, "_ITSxpress Single")]]
  tibble::tibble(
    Comparison = comparison_spec$Comparison, Metric = metric, N_Pairs = sum(complete.cases(x, y)),
    Median_Method_1 = median(x, na.rm = TRUE), Median_Method_2 = median(y, na.rm = TRUE),
    Median_Ratio_Method_2_to_1 = median(y / x, na.rm = TRUE),
    P_value = wilcox.test(y, x, paired = TRUE, exact = FALSE)$p.value
  )
}) |> dplyr::mutate(P_adjust_BH = p.adjust(P_value, "BH"))

qc_pair_plot_data <- sample_qc |> dplyr::select(SampleID, Method, Read_Depth, Observed_ASV) |>
  tidyr::pivot_longer(c(Read_Depth, Observed_ASV), names_to = "Metric", values_to = "Value") |>
  dplyr::mutate(Comparison = comparison_spec$Comparison)
p_qc <- ggplot(qc_pair_plot_data, aes(Method, Value, fill = Method)) +
  geom_boxplot(width = .58, outlier.shape = NA, alpha = .72) +
  geom_line(aes(group = SampleID), colour = "grey70", alpha = .18, linewidth = .25) +
  geom_point(aes(group = SampleID), position = position_jitter(width = .06), size = .55, alpha = .45) +
  facet_grid(Metric ~ Comparison, scales = "free_y") + scale_y_log10(labels = label_comma()) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(x = NULL, y = "Fungal count (log10)", title = "Paired sample QC across complete-forward pipelines") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 28, hjust = 1))
ggsave(file.path(out_dir, "10a_paired_sample_QC.png"), p_qc, width = 10, height = 7.2, dpi = 210, bg = "white")

clr <- function(x) { z <- log(t(x) + .5); z - rowMeans(z) }
clr_list <- purrr::map(aligned, clr)
shape_values <- c("Air swab" = 0, "Clean bench swab" = 1, "Empty tube" = 2, "KSY" = 5, "PYH" = 6, "Mock" = 4)

analyse_beta_pair <- function(Method_1, Method_2) {
  x <- clr_list[[Method_1]]; y <- clr_list[[Method_2]]
  rownames(x) <- paste(Method_1, rownames(x), sep = "__")
  rownames(y) <- paste(Method_2, rownames(y), sep = "__")
  combined <- rbind(x, y); dist_all <- dist(combined)
  meta_pair <- tibble::tibble(RowID = rownames(combined)) |>
    tidyr::separate(RowID, c("Method", "SampleID"), sep = "__", extra = "merge") |>
    dplyr::left_join(metadata |> dplyr::select(SampleID, Shape_Group), by = "SampleID") |>
    dplyr::mutate(Method = factor(Method, levels = c(Method_1, Method_2))) |> data.frame()
  rownames(meta_pair) <- rownames(combined)
  set.seed(42); ad <- vegan::adonis2(dist_all ~ Method, data = meta_pair, permutations = 999, strata = meta_pair$SampleID)
  bd <- vegan::betadisper(dist_all, meta_pair$Method); set.seed(42); bdp <- vegan::permutest(bd, 999)
  ord1 <- cmdscale(dist(x), k = 2, eig = TRUE, add = TRUE)
  ord2 <- cmdscale(dist(y), k = 2, eig = TRUE, add = TRUE)
  set.seed(42); pro <- vegan::protest(ord1$points, ord2$points, permutations = 999)
  pc <- cmdscale(dist_all, k = 2, eig = TRUE, add = TRUE)
  variance <- pc$eig[pc$eig > 0] / sum(pc$eig[pc$eig > 0]) * 100
  colnames(pc$points) <- c("Axis1", "Axis2")
  score <- as_tibble(pc$points, rownames = "RowID") |>
    dplyr::left_join(tibble::rownames_to_column(meta_pair, "RowID"), by = "RowID")
  main <- ggplot(score, aes(Axis1, Axis2, colour = Method, shape = Shape_Group)) +
    geom_hline(yintercept = 0, colour = "grey88") + geom_vline(xintercept = 0, colour = "grey88") +
    geom_line(aes(group = SampleID), colour = "grey76", alpha = .5, linewidth = .28) +
    geom_point(size = 2.1, alpha = .82) +
    scale_colour_manual(values = method_colors) + scale_shape_manual(values = shape_values, drop = TRUE) +
    labs(title = comparison_spec$Comparison, x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
         y = paste0("PCoA2 (", round(variance[2], 1), "%)"), colour = NULL, shape = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  xd <- cowplot::axis_canvas(main, axis = "x") +
    geom_density(data = score, aes(Axis1, fill = Method, colour = Method), alpha = .2, linewidth = .35) +
    scale_fill_manual(values = method_colors) + scale_colour_manual(values = method_colors)
  yd <- cowplot::axis_canvas(main, axis = "y", coord_flip = TRUE) +
    geom_density(data = score, aes(Axis2, fill = Method, colour = Method), alpha = .2, linewidth = .35) +
    scale_fill_manual(values = method_colors) + scale_colour_manual(values = method_colors) + coord_flip()
  plot <- main |> cowplot::insert_xaxis_grob(xd, grid::unit(.58, "in"), position = "top") |>
    cowplot::insert_yaxis_grob(yd, grid::unit(.58, "in"), position = "right") |> cowplot::ggdraw()
  stat <- tibble::tibble(
    Comparison = comparison_spec$Comparison, Method_1, Method_2, N_Pairs = length(common_samples),
    PERMANOVA_R2 = ad$R2[1], PERMANOVA_P = ad$`Pr(>F)`[1],
    Betadisper_P = bdp$tab$`Pr(>F)`[1], Procrustes_r = pro$t0, PROTEST_P = pro$signif
  ) |> dplyr::mutate(PERMANOVA_P_BH = PERMANOVA_P, Betadisper_P_BH = Betadisper_P)
  list(stat = stat, plot = plot)
}
beta_res <- analyse_beta_pair("Cutadapt Single", "ITSxpress Single")
beta_stats <- beta_res$stat
ggsave(file.path(out_dir, "10b_aitchison_complete_forward.png"), beta_res$plot,
       width = 10.5, height = 7.2, dpi = 210, bg = "white")

relative_long <- purrr::imap_dfr(aligned, function(mat, method) {
  rel <- sweep(mat, 2, colSums(mat), "/")
  data.frame(t(rel), check.names = FALSE) |> rownames_to_column("SampleID") |>
    pivot_longer(-SampleID, names_to = "Species", values_to = "Relative_Abundance") |>
    mutate(Method = method)
}) |> left_join(metadata, by = "SampleID")

Method_1 <- "Cutadapt Single"; Method_2 <- "ITSxpress Single"
means <- relative_long |> group_by(Method, Species) |> summarise(Mean = mean(Relative_Abundance), .groups = "drop") |>
  pivot_wider(names_from = Method, values_from = Mean, values_fill = 0)
means$Combined_Mean <- (means[[Method_1]] + means[[Method_2]]) / 2
means <- arrange(means, desc(Combined_Mean)); top_taxa <- head(means$Species, 14)
paired <- relative_long |> filter(Species %in% top_taxa) |> select(SampleID, Species, Method, Relative_Abundance) |>
  pivot_wider(names_from = Method, values_from = Relative_Abundance, values_fill = 0)
taxa_tests <- paired |> group_by(Species) |> summarise(
  Comparison = comparison_spec$Comparison, Method_1, Method_2, N_Pairs = n(),
  Mean_Method_1 = mean(.data[[Method_1]]), Mean_Method_2 = mean(.data[[Method_2]]),
  Mean_Difference_2_minus_1 = mean(.data[[Method_2]] - .data[[Method_1]]),
  P_value = if (all(.data[[Method_1]] == .data[[Method_2]])) 1 else wilcox.test(.data[[Method_2]], .data[[Method_1]], paired = TRUE, exact = FALSE)$p.value,
  .groups = "drop") |> mutate(P_adjust_BH = p.adjust(P_value, "BH")) |> arrange(P_adjust_BH)
detected1 <- rowSums(aligned[[Method_1]]) > 0; detected2 <- rowSums(aligned[[Method_2]]) > 0
taxa_overlap <- tibble::tibble(
  Comparison = comparison_spec$Comparison, Method_1, Method_2,
  Detected_Method_1 = sum(detected1), Detected_Method_2 = sum(detected2),
  Shared = sum(detected1 & detected2), Method_1_only = sum(detected1 & !detected2),
  Method_2_only = sum(!detected1 & detected2)
)
plot_data <- means |> mutate(Label = if_else(Species %in% top_taxa, Species, NA_character_))
p_taxa <- ggplot(plot_data, aes(.data[[Method_1]] + 1e-6, .data[[Method_2]] + 1e-6)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(aes(colour = !is.na(Label)), alpha = .75, size = 1.7) +
  ggrepel::geom_text_repel(data = filter(plot_data, !is.na(Label)), aes(label = Label), seed = 42,
                           size = 2.6, max.overlaps = Inf, min.segment.length = 0) +
  scale_x_log10(labels = label_percent()) + scale_y_log10(labels = label_percent()) +
  scale_colour_manual(values = c("TRUE" = "#1B9E77", "FALSE" = "grey75"), guide = "none") +
  coord_equal() + labs(title = comparison_spec$Comparison,
    subtitle = paste0("Species-level means; Spearman rho = ", round(cor(means[[Method_1]], means[[Method_2]], method = "spearman"), 3)),
    x = Method_1, y = Method_2) + theme_bw(base_size = 10)
ggsave(file.path(out_dir, "10c_species_pairwise_abundance.png"), p_taxa,
       width = 8, height = 7.2, dpi = 210, bg = "white")

overall_means <- relative_long |> group_by(Species) |> summarise(Mean = mean(Relative_Abundance), .groups = "drop") |> arrange(desc(Mean))
top_species <- head(overall_means$Species, 14)
taxon_levels <- c(top_species, "Others")
pal <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2")); taxon_colors <- setNames(pal[seq_along(taxon_levels)], taxon_levels)
sample_order <- metadata |> mutate(Subject_Profile = if_else(Sample_Type == "True", Subject, Control_Type)) |>
  arrange(Subject_Profile, PCR_Cycle_Group, Skin_Site, SampleID) |> pull(SampleID)
composition <- relative_long |> mutate(
  Subject_Profile = case_when(Sample_Type == "True" ~ Subject, Sample_Type == "Control" ~ Control_Type, TRUE ~ "Mock"),
  Species = if_else(Species %in% top_species, Species, "Others"), Species = factor(Species, levels = rev(taxon_levels)),
  SampleID = factor(SampleID, levels = sample_order)) |>
  group_by(Method, SampleID, Subject_Profile, PCR_Cycle_Group, Species) |>
  summarise(Relative_Abundance = sum(Relative_Abundance), .groups = "drop")
draw_composition <- function(method) ggplot(filter(composition, Method == method), aes(SampleID, Relative_Abundance, fill = Species)) +
  geom_col(width = .9) + facet_wrap(vars(Subject_Profile, PCR_Cycle_Group), scales = "free_x", nrow = 1, drop = TRUE) +
  scale_y_continuous(labels = label_percent(), expand = c(0, 0)) + scale_fill_manual(values = taxon_colors, drop = FALSE) +
  labs(title = paste0(method, " — sample-level Species composition"), x = "SampleID", y = "Relative abundance", fill = "Species") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), legend.position = "right")
composition_plots <- purrr::map(method_levels, draw_composition); names(composition_plots) <- method_levels
for (i in seq_along(composition_plots)) {
  ggsave(file.path(out_dir, paste0("10d", i, "_composition_", gsub(" ", "_", method_levels[i]), ".png")),
         composition_plots[[i]], width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
}

sequence_metrics <- purrr::map2_dfr(paths$Method, paths$repseq_path, function(method, path) {
  seqs <- qiime2R::read_qza(path)[["data"]]; ps <- ps_list[[method]]
  otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  ids <- intersect(names(seqs), rownames(otu))
  tibble::tibble(Method = method, FeatureID = ids, Sequence = as.character(seqs[ids]), Length_bp = width(seqs[ids]),
                 Total_Reads = rowSums(otu[ids, , drop = FALSE]), Prevalence = rowSums(otu[ids, , drop = FALSE] > 0)) |>
    group_by(Method) |> mutate(Read_Fraction = Total_Reads / sum(Total_Reads), Abundance_Rank = rank(-Total_Reads, ties.method = "first")) |> ungroup()
})
sequence_summary <- sequence_metrics |> group_by(Method) |> summarise(
  ASVs = n(), Median_Length = median(Length_bp), Q1_Length = quantile(Length_bp, .25), Q3_Length = quantile(Length_bp, .75),
  Min_Length = min(Length_bp), Max_Length = max(Length_bp), Total_Reads = sum(Total_Reads), .groups = "drop")
p_len <- ggplot(sequence_metrics, aes(Length_bp, colour = Method, fill = Method, weight = Total_Reads)) +
  geom_density(alpha = .15, linewidth = .75) + facet_wrap(~Method, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
  labs(x = "ASV length (bp)", y = "Read-weighted density", title = "ASV sequence-length distribution") +
  theme_bw(base_size = 10) + theme(legend.position = "none")
p_rank <- ggplot(sequence_metrics, aes(Abundance_Rank, Read_Fraction, colour = Method)) + geom_point(alpha = .55, size = .9) +
  facet_wrap(~Method, ncol = 2) + scale_x_log10() + scale_y_log10(labels = label_percent()) +
  scale_colour_manual(values = method_colors) + labs(x = "ASV abundance rank", y = "Fraction of fungal reads", title = "ASV rank-abundance") +
  theme_bw(base_size = 10) + theme(legend.position = "none")
p_seq <- p_len + p_rank
ggsave(file.path(out_dir, "10e_ASV_sequence_length_frequency.png"), p_seq,
       width = 14, height = 10, dpi = 210, bg = "white")

raphanus_id <- "5c8f111a2f4a8e9e41a1cb58ebc985f6"
cut_raph <- sequence_metrics |> filter(Method == "Cutadapt Single", FeatureID == raphanus_id)
its_seq <- sequence_metrics |> filter(Method == "ITSxpress Single")
raphanus_matches <- if (nrow(cut_raph)) {
  its_seq |> filter(vapply(Sequence, function(x) grepl(x, cut_raph$Sequence[1], fixed = TRUE), logical(1)))
} else its_seq[0, ]
raphanus_summary <- tibble::tibble(
  Item = c("Cutadapt ASV total reads", "Cutadapt ASV length", "ITSxpress contained ASVs", "ITSxpress matched reads"),
  Value = c(ifelse(nrow(cut_raph), cut_raph$Total_Reads[1], 0), ifelse(nrow(cut_raph), cut_raph$Length_bp[1], NA),
            nrow(raphanus_matches), sum(raphanus_matches$Total_Reads))
)

readr::write_csv(audit, file.path(out_dir, "10a_preprocessing_audit.csv"))
readr::write_csv(depth_summary, file.path(out_dir, "10b_depth_summary.csv"))
readr::write_csv(paired_qc_tests, file.path(out_dir, "10c_paired_qc_tests.csv"))
readr::write_csv(beta_stats, file.path(out_dir, "10d_beta_statistics.csv"))
readr::write_csv(taxa_overlap, file.path(out_dir, "10e_species_overlap.csv"))
readr::write_csv(taxa_tests, file.path(out_dir, "10f_species_paired_tests.csv"))
readr::write_csv(sequence_summary, file.path(out_dir, "10g_ASV_sequence_summary.csv"))
readr::write_csv(sequence_metrics, file.path(out_dir, "10h_ASV_sequence_metrics.csv"))
readr::write_csv(raphanus_summary, file.path(out_dir, "10i_Raphanus_ASV_check.csv"))
saveRDS(list(ps_list = ps_list, audit = audit, sample_qc = sample_qc, beta_stats = beta_stats,
             taxa_overlap = taxa_overlap, taxa_tests = taxa_tests, sequence_metrics = sequence_metrics),
        file.path(out_dir, "10_complete_forward_pair_objects.rds"))
