suppressPackageStartupMessages({
  library(qiime2R); library(phyloseq); library(dplyr); library(tidyr)
  library(tibble); library(readr); library(stringr); library(purrr)
  library(ggplot2); library(vegan); library(cowplot); library(patchwork)
  library(scales); library(RColorBrewer); library(Biostrings)
})
set.seed(42)

out_dir <- file.path("Results", "01_Preprocessing_2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("Phyloseq", showWarnings = FALSE)

method_levels <- c("Cutadapt Single", "Cutadapt Paired", "ITSxpress Single", "ITSxpress Paired")
method_colors <- c(
  "Cutadapt Single" = "#1F78B4", "Cutadapt Paired" = "#A6CEE3",
  "ITSxpress Single" = "#E31A1C", "ITSxpress Paired" = "#FB9A99"
)
comparison_spec <- tibble::tribble(
  ~Comparison, ~Method_1, ~Method_2,
  "Cutadapt: Single vs Paired", "Cutadapt Single", "Cutadapt Paired",
  "ITSxpress: Single vs Paired", "ITSxpress Single", "ITSxpress Paired",
  "Single-end: Cutadapt vs ITSxpress", "Cutadapt Single", "ITSxpress Single",
  "Paired-end: Cutadapt vs ITSxpress", "Cutadapt Paired", "ITSxpress Paired"
)

metadata <- data.frame(
  phyloseq::sample_data(readRDS("Phyloseq/phy_ITS2_cutadapt_q20_ee22_20261015.rds")),
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

paths <- tibble::tribble(
  ~Method, ~table_path, ~taxonomy_path, ~repseq_path,
  "Cutadapt Single", "benchmark_results/q20_ee22_both_full/dada2/cutadapt/single/table.qza", "benchmark_results/q20_ee22_both_full/taxonomy/cutadapt/single/unite_v10_99_full/taxonomy.tsv", "benchmark_results/q20_ee22_both_full/dada2/cutadapt/single/repseq.qza",
  "Cutadapt Paired", "benchmark_results/q20_ee22_both_full/dada2/cutadapt/paired/table.qza", "benchmark_results/q20_ee22_both_full/taxonomy/cutadapt/paired/unite_v10_99_full/taxonomy.tsv", "benchmark_results/q20_ee22_both_full/dada2/cutadapt/paired/repseq.qza",
  "ITSxpress Single", "benchmark_results/q20_ee22_both_full/dada2/itsxpress/single/table.qza", "benchmark_results/q20_ee22_both_full/taxonomy/itsxpress/single/unite_v10_99_full/taxonomy.tsv", "benchmark_results/q20_ee22_both_full/dada2/itsxpress/single/repseq.qza",
  "ITSxpress Paired", "benchmark_results/q20_ee22_both_full/dada2/itsxpress/paired/table.qza", "benchmark_results/q20_ee22_both_full/taxonomy/itsxpress/paired/unite_v10_99_full/taxonomy.tsv", "benchmark_results/q20_ee22_both_full/dada2/itsxpress/paired/repseq.qza"
)

results <- purrr::pmap(
  paths[, c("Method", "table_path", "taxonomy_path")],
  function(Method, table_path, taxonomy_path) make_phyloseq(Method, table_path, taxonomy_path)
)
names(results) <- paths$Method
ps_list <- purrr::map(results, "cleaned")
audit <- purrr::map_dfr(results, "audit")

rds_names <- c(
  "Cutadapt Single" = "phy_ITS2_cutadapt_single_q20_ee22_20261016.rds",
  "Cutadapt Paired" = "phy_ITS2_cutadapt_paired_q20_ee22_20261016.rds",
  "ITSxpress Single" = "phy_ITS2_itsxpress_single_q20_ee22_20261016.rds",
  "ITSxpress Paired" = "phy_ITS2_itsxpress_paired_q20_ee22_20261016.rds"
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

# Preprocessing/rank/read-depth summaries
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

paired_qc_tests <- purrr::pmap_dfr(comparison_spec, function(Comparison, Method_1, Method_2) {
  w <- sample_qc |> dplyr::filter(Method %in% c(Method_1, Method_2)) |>
    dplyr::select(SampleID, Method, Read_Depth, Observed_ASV) |>
    tidyr::pivot_wider(names_from = Method, values_from = c(Read_Depth, Observed_ASV))
  purrr::map_dfr(c("Read_Depth", "Observed_ASV"), function(metric) {
    x <- w[[paste0(metric, "_", Method_1)]]; y <- w[[paste0(metric, "_", Method_2)]]
    tibble::tibble(Comparison, Metric = metric, N_Pairs = sum(complete.cases(x, y)),
                   Median_Method_1 = median(x, na.rm = TRUE), Median_Method_2 = median(y, na.rm = TRUE),
                   Median_Ratio_Method_2_to_1 = median(y / x, na.rm = TRUE),
                   P_value = wilcox.test(y, x, paired = TRUE, exact = FALSE)$p.value)
  })
}) |> dplyr::group_by(Metric) |> dplyr::mutate(P_adjust_BH = p.adjust(P_value, "BH")) |> dplyr::ungroup()

sample_qc_long <- sample_qc |> dplyr::select(SampleID, Method, Read_Depth, Observed_ASV)
qc_pair_plot_data <- purrr::pmap_dfr(comparison_spec, function(Comparison, Method_1, Method_2) {
  sample_qc_long |> dplyr::filter(Method %in% c(Method_1, Method_2)) |>
    tidyr::pivot_longer(c(Read_Depth, Observed_ASV), names_to = "Metric", values_to = "Value") |>
    dplyr::mutate(Comparison = Comparison, Method = factor(Method, levels = c(Method_1, Method_2)))
})
p_qc <- ggplot(qc_pair_plot_data, aes(Method, Value, fill = Method)) +
  geom_boxplot(width = .58, outlier.shape = NA, alpha = .72) +
  geom_line(aes(group = SampleID), colour = "grey70", alpha = .18, linewidth = .25) +
  geom_point(aes(group = SampleID), position = position_jitter(width = .06), size = .55, alpha = .45) +
  facet_grid(Metric ~ Comparison, scales = "free_y") + scale_y_log10(labels = label_comma()) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(x = NULL, y = "Fungal count (log10)", title = "Paired sample QC across four pipelines") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 28, hjust = 1))
ggsave(file.path(out_dir, "08a_paired_sample_QC.png"), p_qc, width = 16, height = 7.2, dpi = 210, bg = "white")

# Aitchison pairwise comparisons and density PCoA
clr <- function(x) { z <- log(t(x) + .5); z - rowMeans(z) }
clr_list <- purrr::map(aligned, clr)
shape_values <- c("Air swab" = 0, "Clean bench swab" = 1, "Empty tube" = 2, "KSY" = 5, "PYH" = 6, "Mock" = 4)

analyse_beta_pair <- function(Comparison, Method_1, Method_2) {
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
  score <- as_tibble(pc$points, rownames = "RowID") |> dplyr::rename(Axis1 = V1, Axis2 = V2) |>
    dplyr::left_join(tibble::rownames_to_column(meta_pair, "RowID"), by = "RowID")
  main <- ggplot(score, aes(Axis1, Axis2, colour = Method, shape = Shape_Group)) +
    geom_hline(yintercept = 0, colour = "grey88") + geom_vline(xintercept = 0, colour = "grey88") +
    geom_line(aes(group = SampleID), colour = "grey76", alpha = .5, linewidth = .28) +
    geom_point(size = 2.1, alpha = .82) +
    scale_colour_manual(values = method_colors) + scale_shape_manual(values = shape_values, drop = TRUE) +
    labs(title = Comparison, x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
         y = paste0("PCoA2 (", round(variance[2], 1), "%)"), colour = NULL, shape = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  xd <- cowplot::axis_canvas(main, axis = "x") + geom_density(data = score, aes(Axis1, fill = Method, colour = Method), alpha = .2, linewidth = .35) +
    scale_fill_manual(values = method_colors) + scale_colour_manual(values = method_colors)
  yd <- cowplot::axis_canvas(main, axis = "y", coord_flip = TRUE) + geom_density(data = score, aes(Axis2, fill = Method, colour = Method), alpha = .2, linewidth = .35) +
    scale_fill_manual(values = method_colors) + scale_colour_manual(values = method_colors) + coord_flip()
  plot <- main |> cowplot::insert_xaxis_grob(xd, grid::unit(.58, "in"), position = "top") |>
    cowplot::insert_yaxis_grob(yd, grid::unit(.58, "in"), position = "right") |> cowplot::ggdraw()
  stat <- tibble::tibble(Comparison, Method_1, Method_2, N_Pairs = length(common_samples),
                         PERMANOVA_R2 = ad$R2[1], PERMANOVA_P = ad$`Pr(>F)`[1],
                         Betadisper_P = bdp$tab$`Pr(>F)`[1], Procrustes_r = pro$t0, PROTEST_P = pro$signif)
  list(stat = stat, plot = plot)
}
beta_res <- purrr::pmap(comparison_spec, analyse_beta_pair)
beta_stats <- purrr::map_dfr(beta_res, "stat") |>
  dplyr::mutate(PERMANOVA_P_BH = p.adjust(PERMANOVA_P, "BH"), Betadisper_P_BH = p.adjust(Betadisper_P, "BH"))
for (i in seq_along(beta_res)) ggsave(file.path(out_dir, paste0("08b", i, "_aitchison_", i, ".png")), beta_res[[i]]$plot, width = 10.5, height = 7.2, dpi = 210, bg = "white")

global_clr <- purrr::imap_dfr(clr_list, function(x, method) {
  data.frame(x, check.names = FALSE) |> tibble::rownames_to_column("SampleID") |> dplyr::mutate(Method = method)
})
global_meta <- global_clr |> dplyr::select(SampleID, Method) |>
  dplyr::left_join(metadata |> dplyr::select(SampleID, Shape_Group), by = "SampleID")
global_matrix <- global_clr |> dplyr::select(-SampleID, -Method) |> as.matrix()
global_ord <- cmdscale(dist(global_matrix), k = 2, eig = TRUE, add = TRUE)
global_var <- global_ord$eig[global_ord$eig > 0] / sum(global_ord$eig[global_ord$eig > 0]) * 100
colnames(global_ord$points) <- c("Axis1", "Axis2")
global_score <- dplyr::bind_cols(global_meta, tibble::as_tibble(global_ord$points)) |>
  dplyr::mutate(Method = factor(Method, levels = method_levels))
p_global <- ggplot(global_score, aes(Axis1, Axis2, colour = Method, shape = Shape_Group)) +
  geom_hline(yintercept = 0, colour = "grey88") + geom_vline(xintercept = 0, colour = "grey88") +
  geom_point(size = 2.1, alpha = .75) + stat_ellipse(aes(group = Method), linetype = 2, linewidth = .55) +
  scale_colour_manual(values = method_colors) + scale_shape_manual(values = shape_values, drop = TRUE) +
  labs(title = "Four-pipeline overview — Aitchison PCoA", x = paste0("PCoA1 (", round(global_var[1], 1), "%)"),
       y = paste0("PCoA2 (", round(global_var[2], 1), "%)"), colour = NULL, shape = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "08b0_aitchison_four_pipeline_overview.png"), p_global, width = 10.5, height = 7.2, dpi = 210, bg = "white")

# Supervised two-factor contrast plane. These axes are method contrasts rather
# than unconstrained PCoA axes: x = Paired -> Single, y = Cutadapt -> ITSxpress
# after removing the component parallel to the read-mode contrast.
factor_meta <- global_meta |>
  dplyr::mutate(
    Read_Mode = if_else(stringr::str_detect(Method, "Single$"), "Single", "Paired"),
    Trimming = if_else(stringr::str_detect(Method, "^ITSxpress"), "ITSxpress", "Cutadapt"),
    Pipeline = factor(Method, levels = method_levels)
  )
centered_global <- sweep(global_matrix, 2, colMeans(global_matrix), "-")
read_vector <- colMeans(global_matrix[factor_meta$Read_Mode == "Single", , drop = FALSE]) -
  colMeans(global_matrix[factor_meta$Read_Mode == "Paired", , drop = FALSE])
trim_vector <- colMeans(global_matrix[factor_meta$Trimming == "ITSxpress", , drop = FALSE]) -
  colMeans(global_matrix[factor_meta$Trimming == "Cutadapt", , drop = FALSE])
read_axis <- read_vector / sqrt(sum(read_vector^2))
trim_residual <- trim_vector - sum(trim_vector * read_axis) * read_axis
trim_axis <- trim_residual / sqrt(sum(trim_residual^2))
contrast_score <- dplyr::bind_cols(
  factor_meta,
  tibble::tibble(
    Paired_to_Single = as.numeric(centered_global %*% read_axis),
    Cutadapt_to_ITSxpress = as.numeric(centered_global %*% trim_axis)
  )
)
if (with(contrast_score, mean(Paired_to_Single[Read_Mode == "Single"]) < mean(Paired_to_Single[Read_Mode == "Paired"]))) {
  contrast_score$Paired_to_Single <- -contrast_score$Paired_to_Single
}
if (with(contrast_score, mean(Cutadapt_to_ITSxpress[Trimming == "ITSxpress"]) < mean(Cutadapt_to_ITSxpress[Trimming == "Cutadapt"]))) {
  contrast_score$Cutadapt_to_ITSxpress <- -contrast_score$Cutadapt_to_ITSxpress
}
contrast_centroid <- contrast_score |>
  dplyr::group_by(Pipeline, Read_Mode, Trimming) |>
  dplyr::summarise(
    Paired_to_Single = mean(Paired_to_Single),
    Cutadapt_to_ITSxpress = mean(Cutadapt_to_ITSxpress), .groups = "drop"
  )

p_factor_plane <- ggplot(
  contrast_score,
  aes(Paired_to_Single, Cutadapt_to_ITSxpress, colour = Pipeline)
) +
  geom_hline(yintercept = 0, colour = "grey78", linewidth = .45) +
  geom_vline(xintercept = 0, colour = "grey78", linewidth = .45) +
  geom_line(aes(group = interaction(SampleID, Trimming)), colour = "grey68", alpha = .10, linewidth = .22) +
  geom_line(aes(group = interaction(SampleID, Read_Mode)), colour = "grey68", alpha = .10, linewidth = .22) +
  geom_point(aes(shape = Shape_Group), size = 1.65, alpha = .48) +
  stat_ellipse(aes(group = Pipeline), linetype = 2, linewidth = .7, level = .80) +
  geom_point(data = contrast_centroid, shape = 21, fill = "white", size = 5, stroke = 1.4) +
  ggrepel::geom_label_repel(
    data = contrast_centroid, aes(label = Pipeline), seed = 42,
    size = 3.1, fontface = "bold", fill = alpha("white", .88), show.legend = FALSE
  ) +
  annotate("text", x = Inf, y = 0, label = "Single  →", hjust = 1.05, vjust = -0.6, colour = "grey35", fontface = "bold") +
  annotate("text", x = -Inf, y = 0, label = "←  Paired", hjust = -0.05, vjust = -0.6, colour = "grey35", fontface = "bold") +
  annotate("text", x = 0, y = Inf, label = "↑  ITSxpress", hjust = -0.08, vjust = 1.3, colour = "grey35", fontface = "bold") +
  annotate("text", x = 0, y = -Inf, label = "↓  Cutadapt", hjust = -0.08, vjust = -0.3, colour = "grey35", fontface = "bold") +
  scale_colour_manual(values = method_colors) + scale_shape_manual(values = shape_values, drop = TRUE) +
  labs(
    title = "Four-pipeline factorial contrast plane",
    subtitle = "Axes are supervised CLR method contrasts; lines connect the same SampleID",
    x = "Read-mode contrast score: Paired  →  Single",
    y = "Trimming contrast score: Cutadapt  →  ITSxpress",
    colour = NULL, shape = NULL
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "08b5_four_pipeline_factorial_contrast_plane.png"), p_factor_plane,
       width = 11.5, height = 8.4, dpi = 220, bg = "white")

factorial_meta <- factor_meta |>
  dplyr::mutate(
    Read_Mode = factor(Read_Mode, levels = c("Paired", "Single")),
    Trimming = factor(Trimming, levels = c("Cutadapt", "ITSxpress"))
  )
set.seed(42)
factorial_adonis <- vegan::adonis2(
  dist(global_matrix) ~ Trimming * Read_Mode,
  data = factorial_meta,
  permutations = 999,
  strata = factorial_meta$SampleID,
  by = "terms"
)
factorial_permanova <- tibble::as_tibble(factorial_adonis, rownames = "Term") |>
  dplyr::filter(!Term %in% c("Residual", "Total")) |>
  dplyr::rename(R2 = R2, P_value = `Pr(>F)`) |>
  dplyr::mutate(P_adjust_BH = p.adjust(P_value, method = "BH"))
four_group_dispersion <- vegan::betadisper(dist(global_matrix), factorial_meta$Pipeline)
set.seed(42)
four_group_dispersion_test <- vegan::permutest(four_group_dispersion, permutations = 999)
factorial_dispersion <- tibble::tibble(
  Test = "Four-pipeline betadisper",
  F_value = four_group_dispersion_test$tab$F[1],
  P_value = four_group_dispersion_test$tab$`Pr(>F)`[1]
)

# Species overlap, abundance tests and scatter plots for allowed pairs
relative_long <- purrr::imap_dfr(aligned, function(mat, method) {
  rel <- sweep(mat, 2, colSums(mat), "/")
  data.frame(t(rel), check.names = FALSE) |> rownames_to_column("SampleID") |>
    pivot_longer(-SampleID, names_to = "Species", values_to = "Relative_Abundance") |>
    mutate(Method = method)
}) |> left_join(metadata, by = "SampleID")

analyse_taxa_pair <- function(Comparison, Method_1, Method_2) {
  sub <- relative_long |> filter(Method %in% c(Method_1, Method_2))
  means <- sub |> group_by(Method, Species) |> summarise(Mean = mean(Relative_Abundance), .groups = "drop") |>
    pivot_wider(names_from = Method, values_from = Mean, values_fill = 0)
  means$Combined_Mean <- (means[[Method_1]] + means[[Method_2]]) / 2
  means <- arrange(means, desc(Combined_Mean)); top <- head(means$Species, 14)
  paired <- sub |> filter(Species %in% top) |> select(SampleID, Species, Method, Relative_Abundance) |>
    pivot_wider(names_from = Method, values_from = Relative_Abundance, values_fill = 0)
  tests <- paired |> group_by(Species) |> summarise(
    Comparison = Comparison, Method_1 = Method_1, Method_2 = Method_2, N_Pairs = n(),
    Mean_Method_1 = mean(.data[[Method_1]]), Mean_Method_2 = mean(.data[[Method_2]]),
    Mean_Difference_2_minus_1 = mean(.data[[Method_2]] - .data[[Method_1]]),
    P_value = if (all(.data[[Method_1]] == .data[[Method_2]])) 1 else wilcox.test(.data[[Method_2]], .data[[Method_1]], paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop") |> mutate(P_adjust_BH = p.adjust(P_value, "BH")) |> arrange(P_adjust_BH)
  detected1 <- rowSums(aligned[[Method_1]]) > 0; detected2 <- rowSums(aligned[[Method_2]]) > 0
  overlap <- tibble::tibble(Comparison, Method_1, Method_2,
                            Detected_Method_1 = sum(detected1), Detected_Method_2 = sum(detected2),
                            Shared = sum(detected1 & detected2), Method_1_only = sum(detected1 & !detected2),
                            Method_2_only = sum(!detected1 & detected2))
  plot_data <- means |> mutate(Label = if_else(Species %in% top, Species, NA_character_))
  p <- ggplot(plot_data, aes(.data[[Method_1]] + 1e-6, .data[[Method_2]] + 1e-6)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_point(aes(colour = !is.na(Label)), alpha = .75, size = 1.7) +
    ggrepel::geom_text_repel(data = filter(plot_data, !is.na(Label)), aes(label = Label), seed = 42,
                             size = 2.6, max.overlaps = Inf, min.segment.length = 0) +
    scale_x_log10(labels = label_percent()) + scale_y_log10(labels = label_percent()) +
    scale_colour_manual(values = c("TRUE" = "#1B9E77", "FALSE" = "grey75"), guide = "none") +
    coord_equal() + labs(title = Comparison,
      subtitle = paste0("Species-level means; Spearman rho = ", round(cor(means[[Method_1]], means[[Method_2]], method = "spearman"), 3)),
      x = Method_1, y = Method_2) + theme_bw(base_size = 10)
  list(tests = tests, overlap = overlap, plot = p, top = top)
}
taxa_res <- purrr::pmap(comparison_spec, analyse_taxa_pair)
taxa_tests <- purrr::map_dfr(taxa_res, "tests")
taxa_overlap <- purrr::map_dfr(taxa_res, "overlap")
p_taxa <- patchwork::wrap_plots(purrr::map(taxa_res, "plot"), ncol = 2) + patchwork::plot_annotation(title = "Allowed pairwise Species abundance comparisons")
ggsave(file.path(out_dir, "08c_species_pairwise_abundance.png"), p_taxa, width = 15, height = 13, dpi = 210, bg = "white")

# Four sample-level composition plots with one shared taxon palette
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
for (i in seq_along(composition_plots)) ggsave(file.path(out_dir, paste0("08d", i, "_composition_", gsub(" ", "_", method_levels[i]), ".png")), composition_plots[[i]], width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)

# ASV sequence length/frequency from the four repseq artifacts
sequence_metrics <- purrr::map2_dfr(paths$Method, paths$repseq_path, function(method, path) {
  seqs <- qiime2R::read_qza(path)[["data"]]; ps <- ps_list[[method]]
  otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  ids <- intersect(names(seqs), rownames(otu))
  tibble::tibble(Method = method, FeatureID = ids, Length_bp = width(seqs[ids]),
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
ggsave(file.path(out_dir, "08e_ASV_sequence_length_frequency.png"), p_seq, width = 14, height = 10, dpi = 210, bg = "white")

# Cross-pipeline ASV sequence families. Raw feature IDs cannot be compared
# directly because trimming/read mode changes sequence boundaries. Families are
# connected components of exact identity or >=100-bp exact containment edges
# across the four pre-specified method comparisons.
asv_nodes <- purrr::map2_dfr(paths$Method, paths$repseq_path, function(method, path) {
  seqs <- qiime2R::read_qza(path)[["data"]]
  ps <- ps_list[[method]]
  otu <- as(phyloseq::otu_table(ps), "matrix"); if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  tax <- data.frame(phyloseq::tax_table(ps), check.names = FALSE)
  ids <- intersect(intersect(names(seqs), rownames(otu)), rownames(tax))
  tibble::tibble(
    NodeID = paste(method, ids, sep = "__"), Method = method, FeatureID = ids,
    Sequence = as.character(seqs[ids]), Length_bp = Biostrings::width(seqs[ids]),
    Total_Reads = as.numeric(rowSums(otu[ids, , drop = FALSE])),
    Taxon_Label = stringr::str_replace_all(as.character(tax[ids, "Species"]), "_", " ")
  )
})

find_containment_edges <- function(method_1, method_2, minimum_bp = 100L) {
  a <- dplyr::filter(asv_nodes, Method == method_1)
  b <- dplyr::filter(asv_nodes, Method == method_2)
  a_set <- Biostrings::DNAStringSet(a$Sequence); b_set <- Biostrings::DNAStringSet(b$Sequence)
  edge_ab <- lapply(which(a$Length_bp >= minimum_bp), function(i) {
    candidate <- which(b$Length_bp >= a$Length_bp[i])
    if (!length(candidate)) return(NULL)
    hit <- candidate[Biostrings::vcountPattern(a_set[[i]], b_set[candidate], fixed = TRUE) > 0]
    if (!length(hit)) return(NULL)
    tibble::tibble(Node_1 = a$NodeID[i], Node_2 = b$NodeID[hit], Match_bp = a$Length_bp[i])
  }) |> dplyr::bind_rows()
  edge_ba <- lapply(which(b$Length_bp >= minimum_bp), function(i) {
    candidate <- which(a$Length_bp > b$Length_bp[i])
    if (!length(candidate)) return(NULL)
    hit <- candidate[Biostrings::vcountPattern(b_set[[i]], a_set[candidate], fixed = TRUE) > 0]
    if (!length(hit)) return(NULL)
    tibble::tibble(Node_1 = a$NodeID[hit], Node_2 = b$NodeID[i], Match_bp = b$Length_bp[i])
  }) |> dplyr::bind_rows()
  dplyr::bind_rows(edge_ab, edge_ba) |> dplyr::distinct(Node_1, Node_2, .keep_all = TRUE) |>
    dplyr::mutate(Comparison = paste(method_1, "vs", method_2))
}

asv_edges <- purrr::pmap_dfr(comparison_spec, function(Comparison, Method_1, Method_2) {
  find_containment_edges(Method_1, Method_2) |> dplyr::mutate(Comparison = Comparison)
})

node_index <- stats::setNames(seq_len(nrow(asv_nodes)), asv_nodes$NodeID)
parent <- seq_len(nrow(asv_nodes))
find_root <- function(x) {
  while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
  x
}
union_nodes <- function(a, b) {
  ra <- find_root(a); rb <- find_root(b)
  if (ra != rb) parent[rb] <<- ra
}
if (nrow(asv_edges)) {
  for (i in seq_len(nrow(asv_edges))) union_nodes(node_index[[asv_edges$Node_1[i]]], node_index[[asv_edges$Node_2[i]]])
}
roots <- vapply(seq_len(nrow(asv_nodes)), find_root, integer(1))
family_lookup <- stats::setNames(sprintf("ASVfam_%04d", match(roots, unique(roots))), asv_nodes$NodeID)
asv_nodes$ASV_Family <- unname(family_lookup[asv_nodes$NodeID])

family_membership <- asv_nodes |>
  dplyr::group_by(ASV_Family) |>
  dplyr::mutate(Family_Members = dplyr::n(), Methods_in_Family = dplyr::n_distinct(Method),
                Max_ASVs_Per_Method = max(table(Method))) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    # Prefer the most resolved available taxonomy as the family label. Choosing
    # only the most abundant member can hide a species-level Cutadapt call behind
    # an abundant ITSxpress member labelled merely Fungi_k.
    Taxonomy_Resolution = dplyr::case_when(
      stringr::str_detect(Taxon_Label, " k$") ~ 1L,
      stringr::str_detect(Taxon_Label, " p$") ~ 2L,
      stringr::str_detect(Taxon_Label, " c$") ~ 3L,
      stringr::str_detect(Taxon_Label, " o$") ~ 4L,
      stringr::str_detect(Taxon_Label, " f$") ~ 5L,
      stringr::str_detect(Taxon_Label, " g$") ~ 6L,
      stringr::str_detect(Taxon_Label, " sp$") ~ 6L,
      TRUE ~ 7L
    )
  )
family_summary <- family_membership |>
  dplyr::group_by(ASV_Family) |>
  dplyr::arrange(dplyr::desc(Taxonomy_Resolution), dplyr::desc(Total_Reads), .by_group = TRUE) |>
  dplyr::slice_head(n = 1) |>
  dplyr::transmute(
    ASV_Family, Representative_FeatureID = FeatureID,
    Representative_Taxon = Taxon_Label, Representative_Length_bp = Length_bp,
    Family_Members, Methods_in_Family, Max_ASVs_Per_Method,
    Ambiguous_Family = Max_ASVs_Per_Method > 1
  ) |>
  dplyr::ungroup()

family_count_list <- purrr::imap(ps_list, function(ps, method) {
  otu <- as(phyloseq::otu_table(ps), "matrix"); if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  ids <- intersect(rownames(otu), family_membership$FeatureID[family_membership$Method == method])
  fam <- family_membership$ASV_Family[match(paste(method, ids, sep = "__"), family_membership$NodeID)]
  rowsum(otu[ids, common_samples, drop = FALSE], group = fam, reorder = FALSE)
})
all_families <- unique(family_membership$ASV_Family)
family_total <- matrix(0, nrow = length(all_families), ncol = length(method_levels),
                       dimnames = list(all_families, method_levels))
for (method in method_levels) {
  x <- family_count_list[[method]]
  family_total[rownames(x), method] <- rowSums(x)
}
method_totals <- colSums(family_total)
family_frequency <- sweep(family_total + 0.5, 2, method_totals + 0.5 * nrow(family_total), "/")
family_effect <- tibble::as_tibble(family_frequency, rownames = "ASV_Family") |>
  dplyr::left_join(family_summary, by = "ASV_Family") |>
  dplyr::mutate(
    Read_Mode_Effect = 0.5 * (log2(`Cutadapt Single`) + log2(`ITSxpress Single`)) -
      0.5 * (log2(`Cutadapt Paired`) + log2(`ITSxpress Paired`)),
    Trimming_Effect = 0.5 * (log2(`ITSxpress Single`) + log2(`ITSxpress Paired`)) -
      0.5 * (log2(`Cutadapt Single`) + log2(`Cutadapt Paired`)),
    Interaction = (log2(`ITSxpress Single`) - log2(`Cutadapt Single`)) -
      (log2(`ITSxpress Paired`) - log2(`Cutadapt Paired`)),
    Mean_Frequency = rowMeans(dplyr::pick(dplyr::all_of(method_levels))),
    Max_Frequency = apply(dplyr::pick(dplyr::all_of(method_levels)), 1, max),
    Raw_Total_Reads = rowSums(family_total[ASV_Family, , drop = FALSE]),
    Quadrant = dplyr::case_when(
      Read_Mode_Effect >= 0 & Trimming_Effect >= 0 ~ "ITSxpress / Single",
      Read_Mode_Effect < 0 & Trimming_Effect >= 0 ~ "ITSxpress / Paired",
      Read_Mode_Effect >= 0 & Trimming_Effect < 0 ~ "Cutadapt / Single",
      TRUE ~ "Cutadapt / Paired"
    )
  )

asv_family_audit <- tibble::tibble(
  Metric = c("Raw ASVs across four pipelines", "Exact-containment edges", "ASV sequence families",
             "Families represented in >=2 pipelines", "Single-pipeline families", "Ambiguous families"),
  N = c(nrow(asv_nodes), nrow(asv_edges), nrow(family_summary),
        sum(family_summary$Methods_in_Family >= 2), sum(family_summary$Methods_in_Family == 1),
        sum(family_summary$Ambiguous_Family))
)
# ASV-family differential abundance: CLR within each sample/pipeline followed by
# paired 2x2 factorial contrasts across the same 83 SampleIDs. No abundance or
# prevalence filter is applied to hypothesis testing.
family_clr_list <- purrr::imap(family_count_list, function(x, method) {
  aligned_count <- matrix(0, nrow = length(all_families), ncol = length(common_samples),
                          dimnames = list(all_families, common_samples))
  aligned_count[rownames(x), colnames(x)] <- x
  logged <- log2(t(aligned_count) + 0.5)
  logged - rowMeans(logged)
})
read_contrast <- 0.5 * (family_clr_list[["Cutadapt Single"]] + family_clr_list[["ITSxpress Single"]]) -
  0.5 * (family_clr_list[["Cutadapt Paired"]] + family_clr_list[["ITSxpress Paired"]])
trimming_contrast <- 0.5 * (family_clr_list[["ITSxpress Single"]] + family_clr_list[["ITSxpress Paired"]]) -
  0.5 * (family_clr_list[["Cutadapt Single"]] + family_clr_list[["Cutadapt Paired"]])
interaction_contrast <- (family_clr_list[["ITSxpress Single"]] - family_clr_list[["Cutadapt Single"]]) -
  (family_clr_list[["ITSxpress Paired"]] - family_clr_list[["Cutadapt Paired"]])

test_paired_contrast <- function(mat, contrast_name) {
  tibble::tibble(
    ASV_Family = colnames(mat),
    Contrast = contrast_name,
    Mean_Effect = colMeans(mat),
    Median_Effect = apply(mat, 2, median),
    P_value = apply(mat, 2, function(z) {
      if (all(z == 0)) 1 else suppressWarnings(stats::wilcox.test(z, mu = 0, exact = FALSE)$p.value)
    })
  ) |>
    dplyr::mutate(P_adjust_BH = p.adjust(P_value, method = "BH"))
}
asv_da_long <- dplyr::bind_rows(
  test_paired_contrast(read_contrast, "Read_Mode"),
  test_paired_contrast(trimming_contrast, "Trimming"),
  test_paired_contrast(interaction_contrast, "Interaction")
)
asv_da <- asv_da_long |>
  tidyr::pivot_wider(
    names_from = Contrast,
    values_from = c(Mean_Effect, Median_Effect, P_value, P_adjust_BH),
    names_glue = "{.value}_{Contrast}"
  ) |>
  dplyr::left_join(family_effect, by = "ASV_Family") |>
  dplyr::mutate(
    Main_Effect_Significance = dplyr::case_when(
      P_adjust_BH_Read_Mode < 0.05 & P_adjust_BH_Trimming < 0.05 ~ "Both main effects",
      P_adjust_BH_Read_Mode < 0.05 ~ "Read mode only",
      P_adjust_BH_Trimming < 0.05 ~ "Trimming only",
      TRUE ~ "Not significant"
    ),
    Interaction_Significant = P_adjust_BH_Interaction < 0.05,
    Effect_Quadrant = dplyr::case_when(
      Mean_Effect_Read_Mode >= 0 & Mean_Effect_Trimming >= 0 ~ "ITSxpress / Single",
      Mean_Effect_Read_Mode < 0 & Mean_Effect_Trimming >= 0 ~ "ITSxpress / Paired",
      Mean_Effect_Read_Mode >= 0 & Mean_Effect_Trimming < 0 ~ "Cutadapt / Single",
      TRUE ~ "Cutadapt / Paired"
    )
  )

asv_da_summary <- asv_da |>
  dplyr::count(Main_Effect_Significance, Interaction_Significant, name = "N_ASV_Families") |>
  dplyr::arrange(Main_Effect_Significance, Interaction_Significant)

asv_effect_range <- dplyr::bind_rows(
  tibble::tibble(Contrast = "Read mode", Value = asv_da$Mean_Effect_Read_Mode),
  tibble::tibble(Contrast = "Trimming", Value = asv_da$Mean_Effect_Trimming),
  tibble::tibble(Contrast = "Interaction", Value = asv_da$Mean_Effect_Interaction)
) |>
  dplyr::group_by(Contrast) |>
  dplyr::summarise(
    Minimum = min(Value), P01 = quantile(Value, .01), Median = median(Value),
    P99 = quantile(Value, .99), Maximum = max(Value), .groups = "drop"
  )

aureobasidium_families <- asv_da |>
  dplyr::filter(stringr::str_detect(Representative_Taxon, regex("Aureobasidium", ignore_case = TRUE))) |>
  dplyr::arrange(dplyr::desc(Mean_Frequency))

label_families <- dplyr::bind_rows(
  asv_da |> dplyr::filter(Main_Effect_Significance != "Not significant") |>
    dplyr::slice_max(Max_Frequency, n = 15, with_ties = FALSE),
  asv_da |> dplyr::filter(Interaction_Significant) |>
    dplyr::slice_max(abs(Mean_Effect_Interaction), n = 8, with_ties = FALSE),
  # Always label the four-axis extrema and Aureobasidium families so that an
  # extreme or biologically important family cannot disappear by label ranking.
  asv_da |> dplyr::slice_min(Mean_Effect_Read_Mode, n = 3, with_ties = FALSE),
  asv_da |> dplyr::slice_max(Mean_Effect_Read_Mode, n = 3, with_ties = FALSE),
  asv_da |> dplyr::slice_min(Mean_Effect_Trimming, n = 3, with_ties = FALSE),
  asv_da |> dplyr::slice_max(Mean_Effect_Trimming, n = 3, with_ties = FALSE),
  aureobasidium_families
) |>
  dplyr::distinct(ASV_Family, .keep_all = TRUE) |>
  dplyr::mutate(Label = paste0(ASV_Family, " | ", Representative_Taxon))
da_colors <- c(
  "Both main effects" = "#6A3D9A", "Read mode only" = "#1F78B4",
  "Trimming only" = "#E31A1C", "Not significant" = "grey75"
)
p_asv_quadrant <- ggplot(asv_da, aes(Mean_Effect_Read_Mode, Mean_Effect_Trimming)) +
  geom_hline(yintercept = 0, colour = "grey55") + geom_vline(xintercept = 0, colour = "grey55") +
  geom_point(aes(colour = Main_Effect_Significance, size = Mean_Frequency,
                 shape = Interaction_Significant), alpha = .68) +
  ggrepel::geom_text_repel(
    data = label_families, aes(label = Label, colour = Main_Effect_Significance),
    seed = 42, size = 2.65, max.overlaps = Inf, min.segment.length = 0,
    box.padding = .5, point.padding = .18, show.legend = FALSE
  ) +
  annotate("label", x = Inf, y = Inf, label = "ITSxpress / Single", hjust = 1.05, vjust = 1.2, fill = alpha("white", .85)) +
  annotate("label", x = -Inf, y = Inf, label = "ITSxpress / Paired", hjust = -.05, vjust = 1.2, fill = alpha("white", .85)) +
  annotate("label", x = Inf, y = -Inf, label = "Cutadapt / Single", hjust = 1.05, vjust = -.2, fill = alpha("white", .85)) +
  annotate("label", x = -Inf, y = -Inf, label = "Cutadapt / Paired", hjust = -.05, vjust = -.2, fill = alpha("white", .85)) +
  scale_colour_manual(values = da_colors) +
  scale_size_continuous(range = c(.45, 6), labels = scales::label_percent()) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
  labs(
    title = "ASV-level 2×2 differential-abundance effect map",
    subtitle = sprintf(
      "All values shown (x %.2f to %.2f; y %.2f to %.2f); paired CLR contrasts across 83 SampleIDs",
      min(asv_da$Mean_Effect_Read_Mode), max(asv_da$Mean_Effect_Read_Mode),
      min(asv_da$Mean_Effect_Trimming), max(asv_da$Mean_Effect_Trimming)
    ),
    x = "Mean CLR effect: Paired  ←  0  →  Single",
    y = "Mean CLR effect: Cutadapt  ←  0  →  ITSxpress",
    colour = "BH-significant main effect", size = "Mean read fraction",
    shape = "Interaction BH < 0.05"
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "08f_ASV_family_four_quadrant.png"), p_asv_quadrant,
       width = 13.5, height = 9.8, dpi = 220, bg = "white")

write_csv(audit, file.path(out_dir, "08a_preprocessing_audit.csv")); write_csv(rank_summary, file.path(out_dir, "08b_taxonomy_rank_summary.csv"))
write_csv(depth_summary, file.path(out_dir, "08c_depth_summary.csv")); write_csv(paired_qc_tests, file.path(out_dir, "08d_paired_QC_tests.csv"))
write_csv(beta_stats, file.path(out_dir, "08e_pairwise_beta_statistics.csv")); write_csv(taxa_overlap, file.path(out_dir, "08f_species_overlap.csv"))
write_csv(taxa_tests, file.path(out_dir, "08g_species_paired_tests.csv")); write_csv(sequence_summary, file.path(out_dir, "08h_ASV_sequence_summary.csv"))
write_csv(factorial_permanova, file.path(out_dir, "08i_factorial_PERMANOVA.csv"))
write_csv(factorial_dispersion, file.path(out_dir, "08j_four_pipeline_betadisper.csv"))
write_csv(contrast_centroid, file.path(out_dir, "08k_factorial_contrast_centroids.csv"))
write_csv(asv_family_audit, file.path(out_dir, "08l_ASV_family_mapping_audit.csv"))
write_csv(family_membership, file.path(out_dir, "08m_ASV_family_membership.csv"))
write_csv(asv_da, file.path(out_dir, "08n_ASV_family_factorial_DA.csv"))
write_csv(asv_da_summary, file.path(out_dir, "08o_ASV_factorial_DA_summary.csv"))
write_csv(asv_effect_range, file.path(out_dir, "08p_ASV_effect_ranges.csv"))
write_csv(aureobasidium_families, file.path(out_dir, "08q_Aureobasidium_ASV_families.csv"))
saveRDS(list(ps_list = ps_list, aligned_species = aligned, comparison_spec = comparison_spec,
             sample_qc = sample_qc, beta_stats = beta_stats, taxa_tests = taxa_tests,
             factorial_permanova = factorial_permanova,
             factorial_dispersion = factorial_dispersion,
             contrast_score = contrast_score,
             asv_family_membership = family_membership,
             asv_da = asv_da,
             composition_plots = composition_plots), file.path(out_dir, "08_four_pipeline_objects.rds"))

cat("COMMON_SAMPLES=", length(common_samples), "\n", sep = "")
print(audit); print(beta_stats); print(factorial_permanova); print(factorial_dispersion)
print(asv_family_audit); print(asv_da_summary); print(sequence_summary)
