suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(phangorn)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
})

set.seed(42)

out_dir <- file.path("Results", "01_Preprocessing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

parse_taxonomy <- function(x) {
  parts <- stringr::str_split_fixed(x, ";", 8)[, seq_along(rank_names), drop = FALSE]
  colnames(parts) <- rank_names
  parts <- as.data.frame(parts, stringsAsFactors = FALSE)
  parts[] <- lapply(parts, function(z) {
    z <- stringr::str_replace(z, "^[dkpcofgs]__", "")
    z[is.na(z) | z %in% c("", "NA", "Unassigned", "uncultured")] <- NA_character_
    z
  })
  parts
}

deepest_rank <- function(taxonomy) {
  parsed <- parse_taxonomy(taxonomy)
  apply(parsed, 1, function(z) {
    placeholder <- stringr::str_detect(
      dplyr::coalesce(z, ""),
      stringr::regex("Incertae_sedis|_sp$|^(Fungi|Ascomycota|Basidiomycota)_(phy|cls|ord|fam|gen)_", ignore_case = TRUE)
    )
    observed <- which(!is.na(z) & z != "" & !placeholder)
    if (length(observed) == 0) 0L else max(observed)
  })
}

deepest_label <- function(taxonomy) {
  parsed <- parse_taxonomy(taxonomy)
  apply(parsed, 1, function(z) {
    placeholder <- stringr::str_detect(
      dplyr::coalesce(z, ""),
      stringr::regex("Incertae_sedis|_sp$|^(Fungi|Ascomycota|Basidiomycota)_(phy|cls|ord|fam|gen)_", ignore_case = TRUE)
    )
    observed <- which(!is.na(z) & z != "" & !placeholder)
    if (length(observed) == 0) "Unassigned" else z[max(observed)]
  })
}

asv <- readr::read_csv(
  file.path(out_dir, "05a_unresolved_taxa_ASV_metrics.csv"),
  show_col_types = FALSE
)
matches <- readr::read_csv(
  file.path(out_dir, "05m_all_exact_containment_pairs.csv"),
  show_col_types = FALSE
)

cut_info <- asv |>
  dplyr::filter(Method == "Cutadapt") |>
  dplyr::select(
    Cutadapt_FeatureID = FeatureID,
    Cutadapt_Sequence = Sequence,
    Cutadapt_Raw_Taxonomy = Raw_Taxonomy
  )

its_info <- asv |>
  dplyr::filter(Method == "ITSxpress") |>
  dplyr::select(
    ITSxpress_FeatureID = FeatureID,
    ITSxpress_Sequence = Sequence,
    ITSxpress_Raw_Taxonomy = Raw_Taxonomy
  )

boundary_pairs <- matches |>
  dplyr::filter(
    Cutadapt_Group %in% c("Fungi_k", "Ascomycota_p"),
    Match_Type == "ITSxpress_within_Cutadapt"
  ) |>
  dplyr::left_join(cut_info, by = "Cutadapt_FeatureID") |>
  dplyr::left_join(its_info, by = "ITSxpress_FeatureID") |>
  dplyr::mutate(
    ITS_Start = stringr::str_locate(Cutadapt_Sequence, stringr::fixed(ITSxpress_Sequence))[, "start"],
    ITS_End = ITS_Start + ITSxpress_Length - 1L,
    Left_Flank_bp = ITS_Start - 1L,
    Right_Flank_bp = Cutadapt_Length - ITS_End,
    Total_Flank_bp = Left_Flank_bp + Right_Flank_bp,
    ITS_Fraction_of_Cutadapt = ITSxpress_Length / Cutadapt_Length,
    Cutadapt_Deepest_Rank_N = deepest_rank(Cutadapt_Raw_Taxonomy),
    ITSxpress_Deepest_Rank_N = deepest_rank(ITSxpress_Raw_Taxonomy),
    Cutadapt_Deepest_Rank = dplyr::if_else(
      Cutadapt_Deepest_Rank_N == 0L, "Unassigned", rank_names[Cutadapt_Deepest_Rank_N]
    ),
    ITSxpress_Deepest_Rank = dplyr::if_else(
      ITSxpress_Deepest_Rank_N == 0L, "Unassigned", rank_names[ITSxpress_Deepest_Rank_N]
    ),
    ITSxpress_Deepest_Label = deepest_label(ITSxpress_Raw_Taxonomy),
    Rank_Improvement = ITSxpress_Deepest_Rank_N - Cutadapt_Deepest_Rank_N,
    Classification_Improved = Rank_Improvement > 0,
    Read_Ratio_ITSxpress_to_Cutadapt = ITSxpress_Total_Reads / Cutadapt_Total_Reads
  ) |>
  dplyr::arrange(dplyr::desc(Cutadapt_Total_Reads))

boundary_summary <- boundary_pairs |>
  dplyr::group_by(Cutadapt_Group) |>
  dplyr::summarise(
    Matched_Cutadapt_ASVs = dplyr::n_distinct(Cutadapt_FeatureID),
    Classification_Improved_ASVs = sum(Classification_Improved),
    Improvement_Percent = 100 * mean(Classification_Improved),
    Median_Cutadapt_Length_bp = stats::median(Cutadapt_Length),
    Median_ITSxpress_Length_bp = stats::median(ITSxpress_Length),
    Median_Total_Flank_bp = stats::median(Total_Flank_bp),
    Median_ITS_Fraction_Percent = 100 * stats::median(ITS_Fraction_of_Cutadapt),
    Cutadapt_Reads_in_Matched_ASVs = sum(Cutadapt_Total_Reads),
    Improved_ASV_Cutadapt_Reads = sum(Cutadapt_Total_Reads[Classification_Improved]),
    .groups = "drop"
  )

rank_transition <- boundary_pairs |>
  dplyr::count(
    Cutadapt_Group,
    ITSxpress_Deepest_Rank,
    ITSxpress_Deepest_Label,
    wt = Cutadapt_Total_Reads,
    name = "Cutadapt_Reads"
  ) |>
  dplyr::arrange(Cutadapt_Group, dplyr::desc(Cutadapt_Reads))

boundary_table <- boundary_pairs |>
  dplyr::transmute(
    Cutadapt_FeatureID,
    Cutadapt_Group,
    Cutadapt_Length_bp = Cutadapt_Length,
    ITSxpress_FeatureID,
    ITSxpress_Length_bp = ITSxpress_Length,
    Left_Flank_bp,
    Right_Flank_bp,
    ITS_Fraction_Percent = 100 * ITS_Fraction_of_Cutadapt,
    Cutadapt_Deepest_Rank,
    ITSxpress_Deepest_Rank,
    ITSxpress_Deepest_Label,
    Classification_Improved,
    Cutadapt_Total_Reads,
    ITSxpress_Total_Reads,
    Read_Ratio_ITSxpress_to_Cutadapt
  )

readr::write_csv(boundary_table, file.path(out_dir, "07a_ITS2_boundary_containment_pairs.csv"))
readr::write_csv(boundary_summary, file.path(out_dir, "07b_ITS2_boundary_summary.csv"))
readr::write_csv(rank_transition, file.path(out_dir, "07c_ITS2_taxonomy_rank_transitions.csv"))

plot_pairs <- boundary_pairs |>
  dplyr::filter(Classification_Improved) |>
  dplyr::slice_max(Cutadapt_Total_Reads, n = 12, with_ties = FALSE) |>
  dplyr::mutate(
    Pair_Label = paste0(
      Cutadapt_Group, " → ", ITSxpress_Deepest_Label,
      "  (", scales::comma(Cutadapt_Total_Reads), " reads)"
    ),
    Pair_Label = factor(Pair_Label, levels = rev(Pair_Label))
  )

taxa_levels <- unique(as.character(plot_pairs$ITSxpress_Deepest_Label))
taxa_palette <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(taxa_levels)),
  taxa_levels
)

p_boundary <- ggplot2::ggplot(plot_pairs, ggplot2::aes(y = Pair_Label)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 1, xend = Cutadapt_Length, yend = Pair_Label),
    linewidth = 5, colour = "grey78", lineend = "butt"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = ITS_Start, xend = ITS_End, yend = Pair_Label,
      colour = ITSxpress_Deepest_Label
    ),
    linewidth = 5, lineend = "butt"
  ) +
  ggplot2::geom_point(ggplot2::aes(x = ITS_Start), shape = 124, size = 3) +
  ggplot2::geom_point(ggplot2::aes(x = ITS_End), shape = 124, size = 3) +
  ggplot2::scale_colour_manual(values = taxa_palette) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.04))) +
  ggplot2::labs(
    x = "Position in Cutadapt ASV (bp)", y = NULL, colour = "ITSxpress taxonomy",
    title = "ITSxpress sequence is an exact internal segment of the longer Cutadapt ASV",
    subtitle = "Grey: complete Cutadapt ASV; colour: exact ITSxpress-contained interval"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(out_dir, "07d_ITS2_boundary_exact_containment.png"),
  p_boundary, width = 12, height = 7.2, dpi = 220, bg = "white"
)

p_rank <- boundary_pairs |>
  dplyr::mutate(
    ITSxpress_Deepest_Rank = factor(ITSxpress_Deepest_Rank, levels = rank_names)
  ) |>
  ggplot2::ggplot(ggplot2::aes(Cutadapt_Group, weight = Cutadapt_Total_Reads, fill = ITSxpress_Deepest_Rank)) +
  ggplot2::geom_bar(position = "fill", colour = "white", linewidth = 0.25, width = 0.68) +
  ggplot2::scale_y_continuous(labels = scales::label_percent()) +
  ggplot2::scale_fill_brewer(palette = "Dark2", drop = FALSE) +
  ggplot2::labs(
    x = NULL, y = "Fraction of matched Cutadapt reads",
    fill = "Deepest ITSxpress rank", title = "Taxonomic resolution of exact-contained ITSxpress ASVs"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  file.path(out_dir, "07e_ITS2_taxonomy_resolution_improvement.png"),
  p_rank, width = 9.5, height = 5.6, dpi = 220, bg = "white"
)

# A compact representative alignment/tree is a visual cross-check only. Exact
# containment coordinates above are the direct evidence for boundary trimming.
tree_pairs <- boundary_pairs |>
  dplyr::filter(Classification_Improved) |>
  dplyr::slice_max(Cutadapt_Total_Reads, n = 9, with_ties = FALSE)

tree_sequences <- c(tree_pairs$Cutadapt_Sequence, tree_pairs$ITSxpress_Sequence)
tree_names <- c(
  paste0("Cut|", tree_pairs$Cutadapt_Group, "|", substr(tree_pairs$Cutadapt_FeatureID, 1, 8)),
  paste0("ITS|", tree_pairs$ITSxpress_Deepest_Label, "|", substr(tree_pairs$ITSxpress_FeatureID, 1, 8))
)
names(tree_sequences) <- make.unique(tree_names)
seqs <- Biostrings::DNAStringSet(tree_sequences)
Biostrings::writeXStringSet(seqs, file.path(out_dir, "07f_ITS2_ASV_internal_reference.fasta"))

alignment <- DECIPHER::AlignSeqs(seqs, verbose = FALSE)
Biostrings::writeXStringSet(alignment, file.path(out_dir, "07g_ITS2_ASV_alignment.fasta"))
phy <- phangorn::phyDat(as.matrix(alignment), type = "DNA")
dm <- phangorn::dist.ml(phy)
tree0 <- phangorn::NJ(dm)
fit0 <- phangorn::pml(tree0, data = phy, k = 4)
fit <- suppressWarnings(phangorn::optim.pml(
  fit0,
  model = "GTR",
  optGamma = TRUE,
  rearrangement = "stochastic",
  control = phangorn::pml.control(trace = 0)
))
ape::write.tree(fit$tree, file = file.path(out_dir, "07h_ITS2_representative_ML_tree.nwk"))

grDevices::png(
  file.path(out_dir, "07i_ITS2_representative_ML_tree.png"),
  width = 2200, height = 1500, res = 220, bg = "white"
)
graphics::par(mar = c(2.5, 1.5, 5, 1.5))
method_tip <- ifelse(stringr::str_starts(fit$tree$tip.label, "Cut\\|"), "Cutadapt", "ITSxpress")
tip_col <- c("Cutadapt" = "#0072B2", "ITSxpress" = "#D55E00")[method_tip]
ape::plot.phylo(
  fit$tree,
  type = "phylogram",
  cex = 0.72,
  tip.color = tip_col,
  no.margin = FALSE,
  main = "Representative exact-contained ASVs: GTR+Gamma ML tree"
)
ape::add.scale.bar(cex = 0.7)
graphics::mtext(
  "Blue = Cutadapt full ASV; orange = ITSxpress internal ASV",
  side = 3, line = 0.4, cex = 0.8
)
grDevices::dev.off()

tree_fit_summary <- tibble::tibble(
  Representative_Pairs = nrow(tree_pairs),
  Tree_Tips = length(fit$tree$tip.label),
  Alignment_Columns = Biostrings::width(alignment)[1],
  Model = "GTR+Gamma",
  Log_Likelihood = as.numeric(stats::logLik(fit))
)
readr::write_csv(tree_fit_summary, file.path(out_dir, "07j_ITS2_tree_fit_summary.csv"))

cat("\nBOUNDARY SUMMARY\n")
print(boundary_summary, n = Inf)
cat("\nTOP IMPROVED PAIRS\n")
print(boundary_table |> dplyr::filter(Classification_Improved) |> dplyr::slice_head(n = 12), n = Inf)
cat("\nTREE SUMMARY\n")
print(tree_fit_summary)
