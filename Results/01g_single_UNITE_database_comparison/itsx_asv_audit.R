# ASV-level ITSx audit; sourced after the existing 01g analysis.
ix_dir <- file.path(root, 'ITSx_res')
ix_out <- file.path(out_dir, 'ITSx_ASV_audit')
dir.create(ix_out, recursive = TRUE, showWarnings = FALSE)
ix_read <- function(path) {
  x <- Biostrings::readDNAStringSet(path)
  ids <- sub(' .*$', '', names(x))
  stopifnot(!anyDuplicated(ids))
  tibble(FeatureID = ids, Sequence = as.character(x), Length = nchar(as.character(x)))
}
ix_input <- ix_read(file.path(ix_dir, 'input/dna-sequences.fasta'))
ix_default <- ix_read(file.path(ix_dir, 'default/itsx.ITS2.fasta'))
ix_partial <- ix_read(file.path(ix_dir, 'partial50/itsx.ITS2.full_and_partial.fasta'))
ix_original <- qiime2R::read_qza(paths$RepSeq[1])$data
stopifnot(setequal(ix_input$FeatureID, names(ix_original)),
          all(ix_input$Sequence == as.character(ix_original[ix_input$FeatureID])),
          all(ix_default$FeatureID %in% ix_input$FeatureID),
          all(ix_partial$FeatureID %in% ix_input$FeatureID))
ix_counts <- res[['Cutadapt Single']]$counts
# QIIME imports can retain one-dimensional array attributes on dimnames.
dimnames(ix_counts) <- lapply(dimnames(ix_counts), as.character)
stopifnot(setequal(rownames(ix_counts), ix_input$FeatureID),
          all(colnames(ix_counts) %in% metadata$SampleID))
ix_positions <- readr::read_tsv(file.path(ix_dir, 'default/itsx.positions.txt'),
  col_names = c('FeatureID','Input_length','SSU','ITS1','S58','ITS2','LSU','Note'),
  col_types = cols(.default = col_character()), show_col_types = FALSE)
ix_audit <- ix_input |>
  rename(Original_Sequence = Sequence, Original_Length = Length) |>
  left_join(ix_default |> rename(Default_Sequence = Sequence, Default_Length = Length), by = 'FeatureID') |>
  left_join(ix_partial |> rename(Partial50_Sequence = Sequence, Partial50_Length = Length), by = 'FeatureID') |>
  left_join(ix_positions, by = 'FeatureID') |>
  mutate(Total_Reads = as.numeric(rowSums(ix_counts)[FeatureID]),
    Status = case_when(!is.na(Default_Sequence) & !is.na(Partial50_Sequence) ~ 'Both',
      !is.na(Partial50_Sequence) ~ 'Partial50 only', !is.na(Default_Sequence) ~ 'Default only',
      TRUE ~ 'No ITS2 output'),
    Same_Sequence = !is.na(Default_Sequence) & !is.na(Partial50_Sequence) & Default_Sequence == Partial50_Sequence,
    Boundary = case_when(is.na(Default_Sequence) ~ 'No default ITS2 output',
      str_detect(coalesce(LSU, ''), 'Not found|No start|No end') ~ 'ITS2 output; incomplete LSU anchor',
      TRUE ~ 'ITS2 output; inspect positions'),
    Removed_bp = Original_Length - Default_Length)
ix_check <- tibble(Check = c('Input ASV IDs and sequences match Cutadapt repseq',
  'Default and partial50 ID sets identical','All shared extracted sequences identical'),
  Passed = c(TRUE, setequal(ix_default$FeatureID, ix_partial$FeatureID),
    all(ix_audit$Same_Sequence[ix_audit$Status == 'Both'])))
ix_stage_ids <- list('Cutadapt input ASVs' = ix_input$FeatureID,
  'ITSx default' = ix_default$FeatureID, 'ITSx partial50' = ix_partial$FeatureID)
ix_summary <- imap_dfr(ix_stage_ids, function(ids, stage) {
  seqs <- if(stage == 'Cutadapt input ASVs') ix_input else if(stage == 'ITSx default') ix_default else ix_partial
  tibble(Stage = stage, Original_ASVs = length(ids), Unique_sequences = n_distinct(seqs$Sequence),
    Linked_reads = sum(ix_counts[ids, , drop = FALSE]), ASV_retention = length(ids)/nrow(ix_counts),
    Read_retention = sum(ix_counts[ids, , drop = FALSE])/sum(ix_counts))
})
ix_sample <- imap_dfr(ix_stage_ids, function(ids, stage) {
  tibble(SampleID = colnames(ix_counts), Stage = stage,
    Reads = colSums(ix_counts[ids, , drop = FALSE]), Input_reads = colSums(ix_counts))
}) |> left_join(metadata, by = 'SampleID') |>
  mutate(Retention = Reads/Input_reads, Stage = factor(Stage, levels = names(ix_stage_ids)))
ix_boundary_summary <- ix_audit |> group_by(Status, Boundary) |>
  summarise(ASVs = n(), Linked_reads = sum(Total_Reads), .groups = 'drop')
ix_tax <- imap_dfr(c(default = 'ITSx default', partial50 = 'ITSx partial50'), function(stage, branch) {
  bind_rows(parse_taxonomy(file.path(ix_dir, branch, 'taxonomy/unite_v10_99_full/taxonomy.tsv'), db_levels[1]),
    parse_taxonomy(file.path(ix_dir, branch, 'taxonomy/unite_v10_dynamic_all_full/taxonomy.tsv'), db_levels[2])) |>
    mutate(Stage = stage)
})
stopifnot(!anyDuplicated(ix_tax[c('Stage','DB','FeatureID')]))
for (stage in c('ITSx default', 'ITSx partial50')) for (db in db_levels) {
  expected_ids <- ix_stage_ids[[stage]]
  observed_ids <- ix_tax$FeatureID[ix_tax$Stage == stage & ix_tax$DB == db]
  stopifnot(setequal(expected_ids, observed_ids))
}
ix_before <- taxonomy_long |> filter(Method == 'Cutadapt Single') |>
  select(-Method) |> mutate(Stage = 'Cutadapt input ASVs')
ix_tax_all <- bind_rows(ix_before, ix_tax)
ix_transition <- ix_tax |>
  left_join(ix_before |> select(FeatureID, DB, Before_Genus = Genus, Before_Species = Species,
    Before_Kingdom = Kingdom_Group, Before_Label = Display_Taxon), by = c('FeatureID','DB')) |>
  left_join(ix_audit |> select(FeatureID, Total_Reads, Original_Length, Default_Length, Boundary), by = 'FeatureID') |>
  mutate(Genus_changed = coalesce(Before_Genus, 'Unassigned') != coalesce(Genus, 'Unassigned'),
    Species_changed = coalesce(Before_Species, 'Unassigned') != coalesce(Species, 'Unassigned'))
ix_tax_summary <- ix_transition |> group_by(Stage, DB) |>
  summarise(ASVs = n(), Genus_changed_ASVs = sum(Genus_changed),
    Species_changed_ASVs = sum(Species_changed), Fungal_ASVs = sum(Kingdom_Group == 'Fungi'),
    Fungal_linked_reads = sum(Total_Reads[Kingdom_Group == 'Fungi']), .groups = 'drop')
ix_focal_ids <- ix_before |> filter(str_detect(Raw_Taxon, regex('Malassezia|Aureobasidium|Raphanus', TRUE))) |>
  pull(FeatureID) |> unique()
ix_focal <- ix_tax_all |> filter(FeatureID %in% ix_focal_ids) |>
  left_join(ix_audit |> select(FeatureID, Total_Reads, Original_Length, Default_Length, Status), by = 'FeatureID') |>
  select(Stage, DB, FeatureID, Display_Taxon, Kingdom_Group, Total_Reads, Original_Length, Default_Length, Status) |>
  arrange(desc(Total_Reads), FeatureID, DB, Stage)
ix_tax_compare <- ix_tax |> select(FeatureID, DB, Stage, Raw_Taxon) |>
  pivot_wider(names_from = Stage, values_from = Raw_Taxon) |>
  mutate(Same_taxonomy = `ITSx default` == `ITSx partial50`)
ix_check <- bind_rows(ix_check, tibble(Check = 'Shared ASVs: taxonomy labels identical in both DBs',
  Passed = all(ix_tax_compare$Same_taxonomy[ix_tax_compare$FeatureID %in% intersect(ix_default$FeatureID, ix_partial$FeatureID)])))
ix_retention_plot <- ix_sample |> filter(Stage != 'Cutadapt input ASVs') |>
  ggplot(aes(SampleID, Retention, fill = Stage)) + geom_col(position = 'dodge', width = .75) +
  facet_wrap(vars(Sample_Type, Subject_Group), scales = 'free_x', nrow = 2) +
  scale_fill_manual(values = c('ITSx default' = '#33A02C', 'ITSx partial50' = '#6A3D9A')) +
  scale_y_continuous(labels = percent, limits = c(0,1)) + soyeon_theme +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1)) +
  labs(x = 'SampleID', y = 'Retained original-ASV reads / Cutadapt ASV reads', fill = NULL)
ix_length_plot <- ix_audit |> filter(!is.na(Default_Length)) |>
  ggplot(aes(Original_Length, Default_Length, colour = Boundary, size = Total_Reads)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey60') +
  geom_point(alpha = .5) + scale_colour_brewer(palette = 'Dark2') + soyeon_theme +
  labs(x = 'Original ASV length (bp)', y = 'ITSx-extracted ITS2 length (bp)', size = 'Original reads')
# Shared original sample denominator shows extraction loss explicitly.
ix_count_long <- as.data.frame(ix_counts) |> rownames_to_column('FeatureID') |>
  pivot_longer(-FeatureID, names_to = 'SampleID', values_to = 'Reads') |> filter(Reads > 0)
ix_comp <- ix_tax_all |> select(FeatureID, DB, Stage, Kingdom_Group, Genus) |>
  inner_join(ix_count_long, by = 'FeatureID', relationship = 'many-to-many') |>
  mutate(Taxon = case_when(is.na(Kingdom_Group) | Kingdom_Group == 'Unassigned' ~ 'Unassigned kingdom',
    Kingdom_Group != 'Fungi' ~ paste0('Non-fungal: ', Kingdom_Group),
    is.na(Genus) ~ 'Unclassified fungi', TRUE ~ Genus)) |>
  group_by(Stage, DB, SampleID, Taxon) |> summarise(Reads = sum(Reads), .groups = 'drop')
ix_top <- ix_comp |> group_by(Taxon) |> summarise(Reads = sum(Reads), .groups = 'drop') |>
  slice_max(Reads, n = 14, with_ties = FALSE) |> pull(Taxon)
ix_comp <- ix_comp |> mutate(Taxon = if_else(Taxon %in% ix_top, Taxon, 'Others')) |>
  group_by(Stage, DB, SampleID, Taxon) |> summarise(Reads = sum(Reads), .groups = 'drop')
ix_missing <- expand_grid(Stage = names(ix_stage_ids), DB = db_levels, SampleID = colnames(ix_counts)) |>
  left_join(ix_comp |> group_by(Stage, DB, SampleID) |> summarise(Kept = sum(Reads), .groups = 'drop'),
    by = c('Stage','DB','SampleID')) |>
  mutate(Reads = as.numeric(colSums(ix_counts)[SampleID]) - coalesce(Kept, 0), Taxon = 'No ITS2 output') |>
  select(Stage, DB, SampleID, Taxon, Reads)
stopifnot(all(ix_missing$Reads >= 0))
ix_comp <- bind_rows(ix_comp, ix_missing) |> left_join(metadata, by = 'SampleID') |>
  mutate(Relative = Reads/as.numeric(colSums(ix_counts)[SampleID]),
    Stage = factor(Stage, levels = names(ix_stage_ids)))
ix_levels <- c(sort(ix_top), 'Others', 'No ITS2 output')
ix_cols <- setNames(c(brewer.pal(12, 'Paired'), brewer.pal(8, 'Dark2'))[seq_along(ix_levels)], ix_levels)
shared <- intersect(names(fungal_genus_cols), names(ix_cols))
ix_cols[shared] <- fungal_genus_cols[shared]
ix_cols[c('Others','No ITS2 output')] <- c('grey75','grey95')
ix_comp_plots <- setNames(lapply(db_levels, function(db) {
  ix_comp |> filter(DB == db) |>
    ggplot(aes(SampleID, Relative, fill = factor(Taxon, levels = ix_levels))) + geom_col(width = .8) +
    facet_grid(Stage ~ Subject_Group, scales = 'free_x', space = 'free_x',
      labeller = labeller(Subject_Group = label_wrap_gen(width = 10))) +
    scale_fill_manual(values = ix_cols, drop = FALSE) + scale_y_continuous(labels = percent) +
    soyeon_theme + theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1)) +
    labs(x = 'SampleID', y = 'Fraction of original Cutadapt ASV reads', fill = 'Genus / status', title = db)
}), db_levels)
ix_exports <- list(ASV_audit = ix_audit, retention_summary = ix_summary, sample_retention = ix_sample,
  boundary_summary = ix_boundary_summary, taxonomy_transitions = ix_transition,
  taxonomy_summary = ix_tax_summary, focal_ASVs = ix_focal, verification = ix_check, composition = ix_comp)
iwalk(ix_exports, ~write_csv(.x, file.path(ix_out, paste0(.y, '.csv'))))
ggsave(file.path(ix_out, 'sample_retention.png'), ix_retention_plot, width = 16, height = 8, dpi = 160)
ggsave(file.path(ix_out, 'length_comparison.png'), ix_length_plot, width = 10, height = 6, dpi = 160)
for (j in seq_along(ix_comp_plots)) ggsave(file.path(ix_out, paste0('composition_DB',j,'.png')),
  ix_comp_plots[[j]], width = 20, height = 10, dpi = 160)
