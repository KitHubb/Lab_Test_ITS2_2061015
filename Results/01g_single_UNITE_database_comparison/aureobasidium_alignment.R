# Coordinate-preserving alignment of the same ASV before/after ITSx.
au_dir <- file.path('Results', '01g_single_UNITE_database_comparison', 'Aureobasidium_alignment')
dir.create(au_dir, recursive = TRUE, showWarnings = FALSE)
au_audit <- readr::read_csv(file.path('Results', '01g_single_UNITE_database_comparison',
  'ITSx_ASV_audit', 'ASV_audit.csv'), show_col_types = FALSE)
au_tx <- readr::read_csv(file.path('Results', '01g_single_UNITE_database_comparison',
  'ITSx_ASV_audit', 'taxonomy_transitions.csv'), show_col_types = FALSE) |>
  dplyr::filter(Before_Genus == 'Aureobasidium')
au_ids <- unique(au_tx$FeatureID)
au_seq <- au_audit |> dplyr::filter(FeatureID %in% au_ids) |>
  dplyr::arrange(dplyr::desc(Total_Reads))
stopifnot(nrow(au_seq) > 0, !anyDuplicated(au_seq$FeatureID))
au_pairs <- purrr::map_dfr(seq_len(nrow(au_seq)), function(i) {
  x <- au_seq[i, ]
  purrr::map_dfr(c('Default', 'Partial50'), function(branch) {
    seq <- x[[paste0(branch, '_Sequence')]]
    stopifnot(!is.na(seq))
    hits <- Biostrings::matchPattern(seq, x$Original_Sequence, fixed = TRUE)
    stopifnot(length(hits) == 1L)
    left <- Biostrings::start(hits)[1] - 1L
    right <- x$Original_Length - Biostrings::end(hits)[1]
    tibble::tibble(FeatureID = x$FeatureID, Stage = branch, Reads = x$Total_Reads,
      Original_bp = x$Original_Length, Extracted_bp = nchar(seq),
      Start_1based = left + 1L, End_1based = left + nchar(seq),
      Removed_5prime_bp = left, Removed_3prime_bp = right,
      Retained_region_identity_pct = 100, Original_coverage_pct = 100*nchar(seq)/x$Original_Length,
      Original = x$Original_Sequence, Extracted = seq,
      Aligned_extracted = paste0(strrep('-', left), seq, strrep('-', right)))
  })
})
au_labels <- au_tx |> dplyr::select(FeatureID, DB, Stage, Before_Label, Display_Taxon, Confidence, Total_Reads)
readr::write_csv(au_pairs, file.path(au_dir, 'same_ASV_alignment.csv'))
readr::write_csv(au_labels, file.path(au_dir, 'taxonomy_labels.csv'))
# All candidate ASVs have the same length and extraction coordinates: direct
# coordinate alignment counts substitutions without an arbitrary gap penalty.
au_default <- au_pairs |> dplyr::filter(Stage == 'Default')
stopifnot(dplyr::n_distinct(au_default$Original_bp) == 1,
  dplyr::n_distinct(au_default$Start_1based) == 1,
  dplyr::n_distinct(au_default$End_1based) == 1)
au_cross <- purrr::map_dfr(seq_len(nrow(au_default)), function(i) {
  purrr::map_dfr(seq_len(nrow(au_default)), function(j) {
    purrr::map_dfr(c('Original', 'Removed 5prime', 'Retained ITS2'), function(region) {
      a <- au_default$Original[i]; b <- au_default$Original[j]
      if (region == 'Removed 5prime') {
        a <- substr(a, 1, au_default$Removed_5prime_bp[i])
        b <- substr(b, 1, au_default$Removed_5prime_bp[j])
      }
      if (region == 'Retained ITS2') {a <- au_default$Extracted[i]; b <- au_default$Extracted[j]}
      mismatches <- sum(strsplit(a, '')[[1]] != strsplit(b, '')[[1]])
      tibble::tibble(ASV1 = au_default$FeatureID[i], ASV2 = au_default$FeatureID[j],
        Region = region, Aligned_bp = nchar(a), Mismatches = mismatches,
        Identity_pct = 100*(1-mismatches/nchar(a)))
    })
  })
})
readr::write_csv(au_cross, file.path(au_dir, 'between_ASV_coordinate_similarity.csv'))
au_order <- au_default$FeatureID
au_heat <- au_cross |>
  dplyr::mutate(ASV1 = factor(ASV1, levels = au_order), ASV2 = factor(ASV2, levels = rev(au_order)),
    Region = factor(Region, levels = c('Original', 'Removed 5prime', 'Retained ITS2'))) |>
  ggplot2::ggplot(ggplot2::aes(ASV1, ASV2, fill = Identity_pct)) +
  ggplot2::geom_tile(colour = 'white') + ggplot2::geom_text(ggplot2::aes(label = Mismatches), size = 2.5) +
  ggplot2::facet_wrap(~Region, nrow = 1) +
  ggplot2::scale_x_discrete(labels = function(x) substr(x, 1, 10)) +
  ggplot2::scale_y_discrete(labels = function(x) substr(x, 1, 10)) +
  ggplot2::scale_fill_distiller(palette = 'Blues', direction = 1) + soyeon_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) +
  ggplot2::labs(x = 'Original ASV ID (first 10 characters)', y = NULL, fill = 'Identity (%)',
    subtitle = 'Numbers in tiles: nucleotide mismatches at aligned coordinates')
ggplot2::ggsave(file.path(au_dir, 'between_ASV_similarity.png'), au_heat, width = 15, height = 6, dpi = 160)
au_colour <- function(s, faded = FALSE) {
  bases <- strsplit(s, '')[[1]]
  cols <- c(A = '#d9f0d3', C = '#c6dbef', G = '#fee391', T = '#fcbba1', '-' = '#eeeeee')
  paste0(vapply(bases, function(b) sprintf('<span style="background:%s;%s">%s</span>',
    ifelse(b %in% names(cols), cols[b], '#ffffff'),
    ifelse(faded, 'opacity:0.6;', ''), b), character(1)), collapse = '')
}
au_html <- '<style>.au-seq{overflow-x:auto;white-space:pre;font:13px/1.8 Consolas,monospace;padding:12px;background:#fafafa;border:1px solid #ddd}.au-card{margin:20px 0}</style>'
for (i in seq_len(nrow(au_default))) {
  x <- au_default[i, ]
  tx <- au_labels |> dplyr::filter(FeatureID == x$FeatureID, Stage == 'ITSx default')
  label <- paste(paste0(tx$DB, ': ', tx$Before_Label, ' → ', tx$Display_Taxon), collapse = '<br>')
  block <- ''
  for (start in seq(1, x$Original_bp, by = 70)) {
    end <- min(start+69, x$Original_bp)
    marks <- paste0(ifelse(strsplit(substr(x$Aligned_extracted, start, end), '')[[1]] == '-', ' ', '|'), collapse = '')
    block <- paste0(block, sprintf('%d–%d\n', start, end),
      'Cutadapt  ', au_colour(substr(x$Original, start, end)), '\n',
      'Match     ', marks, '\n',
      'ITSx      ', au_colour(substr(x$Aligned_extracted, start, end)), '\n\n')
  }
  au_html <- paste0(au_html, '<div class="au-card"><strong>', x$FeatureID, ' — ',
    format(x$Reads, big.mark = ','), ' reads</strong><p>', label,
    '</p><div class="au-seq">', block, '</div></div>')
}
writeLines(au_html, file.path(au_dir, 'nucleotide_alignments.html'), useBytes = TRUE)

# Same original ASV IDs only: no ITSxpress matching in this comparison.
au_conf_before <- dplyr::bind_rows(
  parse_taxonomy(file.path(root, 'taxonomy/cutadapt/single/unite_v10_99_full/taxonomy.tsv'), db_levels[1]),
  parse_taxonomy(file.path(root, 'taxonomy/cutadapt/single/unite_v10_dynamic_all_full/taxonomy.tsv'), db_levels[2])) |>
  dplyr::filter(FeatureID %in% au_ids) |>
  dplyr::select(FeatureID, DB, Cutadapt_label = Display_Taxon,
    Cutadapt_confidence = Confidence, Cutadapt_raw_taxonomy = Raw_Taxon)
au_confidence <- au_tx |>
  dplyr::filter(Stage == 'ITSx default') |>
  dplyr::select(FeatureID, DB, ITSx_label = Display_Taxon,
    ITSx_confidence = Confidence, ITSx_raw_taxonomy = Raw_Taxon) |>
  dplyr::left_join(au_conf_before, by = c('FeatureID', 'DB')) |>
  dplyr::left_join(au_seq |> dplyr::select(FeatureID, Total_Reads, Original_Length, Default_Length), by = 'FeatureID') |>
  dplyr::select(FeatureID, Total_Reads, DB, Original_Length, Default_Length,
    Cutadapt_label, Cutadapt_confidence, ITSx_label, ITSx_confidence,
    Cutadapt_raw_taxonomy, ITSx_raw_taxonomy) |>
  dplyr::arrange(dplyr::desc(Total_Reads), DB)
stopifnot(nrow(au_confidence) == 2*nrow(au_seq))
au_partial_check <- au_tx |>
  dplyr::select(FeatureID, DB, Stage, Raw_Taxon, Confidence) |>
  tidyr::pivot_wider(names_from = Stage, values_from = c(Raw_Taxon, Confidence))
stopifnot(all(au_partial_check[['Raw_Taxon_ITSx default']] == au_partial_check[['Raw_Taxon_ITSx partial50']]),
  all(au_partial_check[['Confidence_ITSx default']] == au_partial_check[['Confidence_ITSx partial50']]))
au_blast_saved <- readr::read_csv('Results/01_Preprocessing/Aureobasidium_tremulum_BLAST_summary.csv', show_col_types = FALSE) |>
  dplyr::filter(grepl('^Cutadapt_', Query))
au_blast_fasta <- Biostrings::readDNAStringSet('Results/01_Preprocessing/Aureobasidium_tremulum_BLAST_queries.fasta')
au_blast_evidence <- purrr::map_dfr(seq_len(nrow(au_seq)), function(i) {
  a <- au_seq[i, ]
  purrr::map_dfr(unique(au_blast_saved$Query), function(query) {
    fasta_name <- paste0(query, '_', au_blast_saved$Query_length_bp[match(query, au_blast_saved$Query)], 'bp')
    stopifnot(fasta_name %in% names(au_blast_fasta))
    bseq <- as.character(au_blast_fasta[[fasta_name]])
    exact <- identical(a$Original_Sequence, bseq)
    contained <- grepl(a$Original_Sequence, bseq, fixed = TRUE)
    tibble::tibble(FeatureID = a$FeatureID, Saved_query = query,
      ASV_length = a$Original_Length, Saved_query_length = nchar(bseq),
      Sequence_relation = if(exact) 'Exact full-query match' else if(contained)
        'ASV is exact substring of longer saved query; NOT direct ASV BLAST' else 'No exact query/containment match',
      Exact_query_match = exact, Contained_in_query = contained)
  })
})
au_confidence <- au_confidence |>
  dplyr::left_join(au_blast_evidence |> dplyr::select(FeatureID, BLAST_evidence = Sequence_relation), by = 'FeatureID')
au_blast_related <- au_blast_evidence |> dplyr::filter(Exact_query_match | Contained_in_query) |>
  dplyr::left_join(au_blast_saved, by = c('Saved_query' = 'Query')) |>
  dplyr::mutate(Bit_score = NA_real_, Alignment_length = NA_integer_,
    Missing_fields = 'Bit score and alignment length not recorded in saved summary')
readr::write_csv(au_confidence, file.path(au_dir, 'Cutadapt_ITSx_confidence_by_DB.csv'))
readr::write_csv(au_blast_evidence, file.path(au_dir, 'BLAST_query_sequence_verification.csv'))
readr::write_csv(au_blast_related, file.path(au_dir, 'BLAST_related_query_hits_NOT_direct_ASV_hits.csv'))

# Batch-BLAST evidence available in blast_results. The external ITSx 173-bp
# sequences are not present there, so they remain explicitly missing.
au_read_blast <- function(path) readr::read_tsv(path, show_col_types = FALSE,
  col_types = readr::cols(.default = readr::col_character())) |>
  dplyr::mutate(dplyr::across(dplyr::any_of(c('QIIME_Confidence','QIIME_ConfidenceNumeric',
    'BLAST_Top1_Pident','BLAST_Top1_Qcovus','BLAST_Top1_Bitscore',
    'BLAST_Top1_AlignmentLength')), as.numeric))
au_blast_cut <- au_read_blast(file.path('blast_results',
  'q20_ee22_cutadapt_forward_single_species_rescue', 'taxonomy_blast_evidence.tsv')) |>
  dplyr::filter(ASV %in% au_ids)
stopifnot(nrow(au_blast_cut) == length(au_ids))
au_blast_ixp_all <- au_read_blast(file.path('blast_results',
  'q20_ee22_itsxpress_forward_single_species_rescue', 'taxonomy_blast_evidence.tsv'))
au_blast_ixp <- au_blast_ixp_all |>
  dplyr::filter(grepl('Aureobasidium', dplyr::coalesce(QIIME_Genus, '')) |
    grepl('Aureobasidium', dplyr::coalesce(BLAST_Top1_Genus, '')))
au_ixp_ids <- au_blast_ixp$ASV
au_ixp_reads <- as.numeric(rowSums(res[['ITSxpress Single']]$counts)[au_ixp_ids])
stopifnot(!anyNA(au_ixp_reads))
au_blast_status <- function(x) dplyr::case_when(
  grepl('blast_species_rescued', x) ~ 'BLAST species rescued',
  grepl('resolved_species', x) ~ 'QIIME species retained',
  grepl('genus_mismatch', x) ~ 'Genus discordant',
  grepl('ambiguous|cutoff_fail|no_blast', x) ~ 'BLAST inconclusive',
  TRUE ~ 'Other')
au_rank <- function(x) dplyr::case_when(
  grepl('^[A-Za-z]+ [a-z-]+$', x) & !grepl(' sp\\.$', x) ~ 'Species label',
  grepl('Aureobasidium', x) ~ 'Genus', TRUE ~ 'Below genus')
au_blast_metrics <- function(x, processing, reads) x |>
  dplyr::transmute(FeatureID = ASV, Processing = processing, Total_Reads = reads,
    BLAST_Identity = BLAST_Top1_Pident, BLAST_Query_Coverage = BLAST_Top1_Qcovus,
    BLAST_Composite = BLAST_Top1_Pident * BLAST_Top1_Qcovus / 100,
    BLAST_Bitscore = BLAST_Top1_Bitscore, BLAST_Alignment_Length = BLAST_Top1_AlignmentLength,
    BLAST_Top1_Species, BLAST_Ambiguous = BLAST_AmbiguousTopHit,
    BLAST_Genus_Match = BLAST_GenusMatch, Replacement_Status,
    Evidence = au_blast_status(Replacement_Status))
au_cut_metrics <- au_blast_metrics(au_blast_cut, 'Cutadapt original',
  au_seq$Total_Reads[match(au_blast_cut$ASV, au_seq$FeatureID)])
au_ixp_metrics <- au_blast_metrics(au_blast_ixp, 'ITSxpress', au_ixp_reads)
au_classifier_cut <- au_confidence |>
  dplyr::transmute(FeatureID, DB, Processing = 'Cutadapt original',
    QIIME_Label = Cutadapt_label, QIIME_Confidence = Cutadapt_confidence)
au_classifier_itsx <- au_confidence |>
  dplyr::transmute(FeatureID, DB, Processing = 'Cutadapt -> ITSx',
    QIIME_Label = ITSx_label, QIIME_Confidence = ITSx_confidence)
au_classifier_ixp <- taxonomy_long |>
  dplyr::filter(Method == 'ITSxpress Single', FeatureID %in% au_ixp_ids) |>
  dplyr::transmute(FeatureID, DB, Processing = 'ITSxpress',
    QIIME_Label = Display_Taxon, QIIME_Confidence = Confidence)
au_blast_plot_data <- dplyr::bind_rows(au_classifier_cut, au_classifier_itsx, au_classifier_ixp) |>
  dplyr::left_join(dplyr::bind_rows(au_cut_metrics, au_ixp_metrics),
    by = c('FeatureID','Processing')) |>
  dplyr::mutate(Total_Reads = dplyr::case_when(
      Processing == 'Cutadapt -> ITSx' ~ au_seq$Total_Reads[match(FeatureID, au_seq$FeatureID)],
      TRUE ~ Total_Reads),
    QIIME_Rank = au_rank(QIIME_Label),
    Processing = factor(Processing,
      levels = c('Cutadapt original','Cutadapt -> ITSx','ITSxpress')),
    DB = factor(DB, levels = db_levels),
    Evidence = factor(Evidence, levels = c('QIIME species retained','BLAST species rescued',
      'BLAST inconclusive','Genus discordant','Other')))
stopifnot(nrow(au_blast_plot_data) == 2*(length(au_ids)+length(au_ids)+length(au_ixp_ids)))
readr::write_csv(au_blast_plot_data, file.path(au_dir,
  'Aureobasidium_BLAST_classifier_scatter_data.csv'))
au_itsx_fasta <- Biostrings::DNAStringSet(au_seq$Default_Sequence)
names(au_itsx_fasta) <- paste0(au_seq$FeatureID, '|ITSx_default|', au_seq$Default_Length, 'bp')
Biostrings::writeXStringSet(au_itsx_fasta, file.path(au_dir,
  'Aureobasidium_ITSx_default_10_queries.fasta'))
au_blast_plot_ready <- au_blast_plot_data |>
  dplyr::filter(!is.na(BLAST_Composite), Processing != 'Cutadapt -> ITSx') |>
  dplyr::mutate(Processing = droplevels(Processing))
au_scatter <- au_blast_plot_ready |>
  ggplot2::ggplot(ggplot2::aes(BLAST_Composite, QIIME_Confidence,
    colour = Evidence, shape = QIIME_Rank, size = log10(Total_Reads + 1))) +
  ggplot2::geom_point(alpha = .8, stroke = .9) +
  ggrepel::geom_text_repel(data = function(x) x |>
      dplyr::filter(Total_Reads >= 100 | Evidence %in% c('BLAST species rescued','Genus discordant')),
    ggplot2::aes(label = substr(FeatureID, 1, 10)), size = 3.5,
    min.segment.length = 0, max.overlaps = Inf, direction = 'y',
    xlim = c(96.8, 100.05), ylim = c(.71, .995), show.legend = FALSE) +
  ggplot2::facet_grid(DB ~ Processing, drop = FALSE) +
  ggplot2::scale_colour_brewer(palette = 'Dark2', drop = TRUE) +
  ggplot2::scale_shape_manual(values = c('Species label' = 16, 'Genus' = 15, 'Below genus' = 17)) +
  ggplot2::scale_size_continuous(range = c(2.5, 9),
    breaks = log10(c(10,100,1000,10000,100000)+1),
    labels = scales::comma(c(10,100,1000,10000,100000)), name = 'Original reads') +
  ggplot2::coord_cartesian(xlim = c(96.7,100.1), ylim = c(.7,1.005)) +
  soyeon_theme +
  ggplot2::guides(colour = ggplot2::guide_legend(order = 1, nrow = 1),
    shape = ggplot2::guide_legend(order = 2, nrow = 1),
    size = ggplot2::guide_legend(order = 3, nrow = 1)) +
  ggplot2::theme(text = ggplot2::element_text(size = 15),
    axis.title = ggplot2::element_text(size = 16),
    axis.text = ggplot2::element_text(size = 13),
    strip.text = ggplot2::element_text(size = 14, face = 'bold'),
    legend.position = 'bottom', legend.box = 'vertical',
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 13),
    legend.key.height = grid::unit(.55, 'cm')) +
  ggplot2::labs(x = 'BLAST coverage-weighted identity = identity x query coverage / 100 (%)',
    y = 'QIIME classifier confidence', colour = 'QIIME-BLAST outcome', shape = 'QIIME rank',
    subtitle = 'Aureobasidium-related ASVs; shared axes are zoomed to the observed high-score range')
ggplot2::ggsave(file.path(au_dir, 'Aureobasidium_BLAST_classifier_scatter.png'),
  au_scatter, width = 12, height = 9, dpi = 180)
