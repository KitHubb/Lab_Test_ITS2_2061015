suppressPackageStartupMessages({
  library(qiime2R); library(phyloseq); library(dplyr); library(tidyr)
  library(tibble); library(readr); library(stringr); library(purrr)
  library(ggplot2); library(scales); library(RColorBrewer); library(patchwork)
  library(vegan); library(DT); library(plotly)
})
set.seed(42)

out_dir <- file.path("Results", "01g_single_UNITE_database_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

root <- file.path("benchmark_results", "q20_ee22_both_full_cutadapt_forward")
method_levels <- c("Cutadapt Single", "ITSxpress Single")
db_levels <- c("UNITE 99% Fungi", "UNITE Dynamic All-euk")
method_cols <- c("Cutadapt Single" = "#1F78B4", "ITSxpress Single" = "#E31A1C")
db_cols <- c("UNITE 99% Fungi" = "#6A3D9A", "UNITE Dynamic All-euk" = "#33A02C")
pcr_levels <- c("PCR1_33__PCR2_8", "PCR1_33__PCR2_10", "PCR1_35__PCR2_8", "PCR1_35__PCR2_10")

soyeon_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey30"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.key = element_blank()
  )

metadata <- data.frame(
  sample_data(readRDS(file.path("Phyloseq", "phy_ITS2_cutadapt_single_complete_q20_ee2_20261016.rds"))),
  check.names = FALSE
) |>
  rownames_to_column("SampleID") |>
  mutate(
    Sample_Type = case_when(
      str_detect(Subject, regex("mock|zymo", ignore_case = TRUE)) ~ "Mock",
      str_detect(Subject, regex("^(ASwab|CSwab|ETube)", ignore_case = TRUE)) ~ "Control",
      TRUE ~ "True"
    ),
    Control_Type = case_when(
      str_detect(Subject, regex("^ASwab", ignore_case = TRUE)) ~ "Air swab",
      str_detect(Subject, regex("^CSwab", ignore_case = TRUE)) ~ "Clean bench swab",
      str_detect(Subject, regex("^ETube", ignore_case = TRUE)) ~ "Empty tube",
      Sample_Type == "Mock" ~ "Mock", TRUE ~ "Skin"
    ),
    Subject_Group = if_else(Sample_Type == "True", as.character(Subject), as.character(Control_Type)),
    Site_Group = case_when(
      Sample_Type != "True" ~ as.character(Control_Type),
      is.na(Direction) | Direction %in% c("", "-") ~ as.character(Region),
      TRUE ~ paste(Region, Direction)
    ),
    PCR_Cycle_Group = factor(PCR_Cycle_Group, levels = pcr_levels)
  )

paths <- tribble(
  ~Method, ~Table, ~RepSeq, ~FungiTax, ~AllTax,
  "Cutadapt Single",
  file.path(root, "dada2", "cutadapt", "single", "table.qza"),
  file.path(root, "dada2", "cutadapt", "single", "repseq.qza"),
  file.path(root, "taxonomy", "cutadapt", "single", "unite_v10_99_full", "taxonomy.tsv"),
  file.path(root, "taxonomy", "cutadapt", "single", "unite_v10_dynamic_all_full", "taxonomy.tsv"),
  "ITSxpress Single",
  file.path(root, "dada2", "itsxpress", "single", "table.qza"),
  file.path(root, "dada2", "itsxpress", "single", "repseq.qza"),
  file.path(root, "taxonomy", "itsxpress", "single", "unite_v10_99_full", "taxonomy.tsv"),
  file.path(root, "taxonomy", "itsxpress", "single", "unite_v10_dynamic_all_full", "taxonomy.tsv")
)
stopifnot(all(file.exists(unlist(paths[, -1]))))

clean_rank <- function(x) {
  x <- str_replace(as.character(x), "^[dkpcofgs]__", "")
  x[is.na(x) | x %in% c("", "NA", "Unassigned", "unidentified", "uncultured")] <- NA_character_
  x
}

parse_taxonomy <- function(path, db) {
  z <- read_tsv(path, show_col_types = FALSE)
  spl <- str_split_fixed(z$Taxon, ";", 8)[, 1:7, drop = FALSE]
  colnames(spl) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  spl <- as.data.frame(spl, stringsAsFactors = FALSE)
  spl[] <- lapply(spl, clean_rank)
  bind_cols(
    tibble(FeatureID = as.character(z$`Feature ID`), DB = db, Raw_Taxon = as.character(z$Taxon), Confidence = as.numeric(z$Confidence)),
    as_tibble(spl)
  ) |>
    mutate(
      Kingdom_Group = case_when(
        is.na(Kingdom) ~ "Unassigned",
        Kingdom == "Fungi" ~ "Fungi",
        Kingdom == "Viridiplantae" ~ "Viridiplantae",
        TRUE ~ Kingdom
      ),
      Informative_Phylum = !is.na(Phylum) & !str_detect(Phylum, regex("incertae|unidentified", TRUE)),
      Informative_Class = !is.na(Class) & !str_detect(Class, regex("incertae|unidentified", TRUE)),
      Informative_Order = !is.na(Order) & !str_detect(Order, regex("incertae|unidentified", TRUE)),
      Informative_Family = !is.na(Family) & !str_detect(Family, regex("incertae|unidentified", TRUE)),
      Informative_Genus = !is.na(Genus) & !str_detect(Genus, regex("incertae|unidentified", TRUE)),
      Informative_Species = !is.na(Species) & !str_detect(Species, regex("incertae|unidentified|_sp$", TRUE)),
      Deepest_Rank = case_when(
        Informative_Species ~ "Species", Informative_Genus ~ "Genus",
        Informative_Family ~ "Family", Informative_Order ~ "Order",
        Informative_Class ~ "Class", Informative_Phylum ~ "Phylum",
        !is.na(Kingdom) ~ "Kingdom", TRUE ~ "Unassigned"
      ),
      Deepest_N = match(Deepest_Rank, c("Unassigned", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) - 1L,
      Display_Taxon = case_when(
        Informative_Species ~ str_replace_all(Species, "_", " "),
        Informative_Genus ~ paste0(Genus, " sp."),
        Informative_Family ~ paste0(Family, " (family)"),
        Informative_Order ~ paste0(Order, " (order)"),
        Informative_Class ~ paste0(Class, " (class)"),
        Informative_Phylum ~ paste0(Phylum, " (phylum)"),
        !is.na(Kingdom) ~ paste0(Kingdom, " (kingdom)"),
        TRUE ~ "Unassigned"
      )
    )
}

read_method <- function(method, table_path, repseq_path, fungi_path, all_path) {
  counts <- qiime2R::read_qza(table_path)$data
  counts <- as(counts, "matrix")
  if (nrow(counts) < ncol(counts) && all(colnames(counts) %in% rownames(metadata))) counts <- t(counts)
  seqs <- qiime2R::read_qza(repseq_path)$data
  seq_tbl <- tibble(FeatureID = names(seqs), Length = nchar(as.character(seqs)))
  feature <- tibble(
    FeatureID = as.character(rownames(counts)), Total_Reads = as.numeric(rowSums(counts)),
    Prevalence = rowSums(counts > 0)
  ) |>
    left_join(seq_tbl, by = "FeatureID")
  fungi <- parse_taxonomy(fungi_path, "UNITE 99% Fungi")
  all_euk <- parse_taxonomy(all_path, "UNITE Dynamic All-euk")
  comparison <- feature |>
    left_join(fungi |> select(-DB) |> rename_with(~paste0(.x, "_Fungi"), -FeatureID), by = "FeatureID") |>
    left_join(all_euk |> select(-DB) |> rename_with(~paste0(.x, "_All"), -FeatureID), by = "FeatureID") |>
    mutate(
      Method = method,
      Kingdom_Transition = paste0(coalesce(Kingdom_Group_Fungi, "Unassigned"), " → ", coalesce(Kingdom_Group_All, "Unassigned")),
      Rank_Change = Deepest_N_All - Deepest_N_Fungi,
      Label_Changed = Raw_Taxon_Fungi != Raw_Taxon_All
    )
  long_tax <- bind_rows(fungi, all_euk) |> mutate(Method = method)
  list(counts = counts, feature = feature, comparison = comparison, taxonomy = long_tax)
}

res <- pmap(
  paths,
  function(Method, Table, RepSeq, FungiTax, AllTax) read_method(Method, Table, RepSeq, FungiTax, AllTax)
)
names(res) <- paths$Method
comparison <- map_dfr(res, "comparison") |>
  mutate(Method = factor(Method, levels = method_levels))
taxonomy_long <- map_dfr(res, "taxonomy") |>
  mutate(Method = factor(Method, levels = method_levels), DB = factor(DB, levels = db_levels))

summary_table <- comparison |>
  group_by(Method) |>
  summarise(
    Total_ASVs = n(),
    AllDB_Fungal_ASVs = sum(Kingdom_Group_All == "Fungi", na.rm = TRUE),
    AllDB_Fungal_Reads = sum(.data$Total_Reads[Kingdom_Group_All == "Fungi"], na.rm = TRUE),
    AllDB_Nonfungal_ASVs = sum(!Kingdom_Group_All %in% c("Fungi", "Unassigned"), na.rm = TRUE),
    AllDB_Nonfungal_Reads = sum(.data$Total_Reads[!Kingdom_Group_All %in% c("Fungi", "Unassigned")], na.rm = TRUE),
    AllDB_Unassigned_ASVs = sum(is.na(Kingdom_Group_All) | Kingdom_Group_All == "Unassigned"),
    Label_Changed_ASVs = sum(Label_Changed, na.rm = TRUE),
    Median_Confidence_FungiDB = median(Confidence_Fungi, na.rm = TRUE),
    Median_Confidence_AllDB = median(Confidence_All, na.rm = TRUE),
    Total_Reads = sum(.data$Total_Reads),
    .groups = "drop"
  ) |>
  mutate(
    Fungal_Read_Retention = AllDB_Fungal_Reads / Total_Reads,
    Nonfungal_Read_Fraction = AllDB_Nonfungal_Reads / Total_Reads
  )

kingdom_summary <- comparison |>
  mutate(AllDB_Kingdom = coalesce(Kingdom_Group_All, "Unassigned")) |>
  group_by(Method, AllDB_Kingdom) |>
  summarise(ASVs = n(), Reads = sum(Total_Reads), .groups = "drop") |>
  group_by(Method) |>
  mutate(ASV_Fraction = ASVs / sum(ASVs), Read_Fraction = Reads / sum(Reads)) |>
  ungroup()

transition_summary <- comparison |>
  group_by(Method, Kingdom_Transition) |>
  summarise(ASVs = n(), Reads = sum(Total_Reads), .groups = "drop") |>
  group_by(Method) |>
  mutate(Read_Fraction = Reads / sum(Reads), ASV_Fraction = ASVs / sum(ASVs)) |>
  ungroup()

rank_transition <- comparison |>
  count(Method, Deepest_Rank_Fungi, Deepest_Rank_All, wt = Total_Reads, name = "Reads") |>
  group_by(Method) |>
  mutate(Read_Fraction = Reads / sum(Reads)) |>
  ungroup()

rank_change_summary <- map_dfr(c("Kingdom", "Phylum", "Genus"), function(rank) {
  old_col <- paste0(rank, "_Fungi")
  new_col <- paste0(rank, "_All")
  comparison |>
    mutate(
      Rank = rank,
      Old_Label = coalesce(as.character(.data[[old_col]]), "Unassigned"),
      New_Label = coalesce(as.character(.data[[new_col]]), "Unassigned"),
      Change_Status = if_else(Old_Label == New_Label, "Same", "Changed")
    ) |>
    group_by(Method, Rank, Change_Status) |>
    summarise(ASVs = n(), Reads = sum(Total_Reads), .groups = "drop") |>
    group_by(Method, Rank) |>
    mutate(ASV_Fraction = ASVs / sum(ASVs), Read_Fraction = Reads / sum(Reads)) |>
    ungroup()
}) |>
  mutate(Rank = factor(Rank, levels = c("Kingdom", "Phylum", "Genus")),
         Change_Status = factor(Change_Status, levels = c("Same", "Changed")))

rank_transition_detail <- map_dfr(c("Kingdom", "Phylum", "Genus"), function(rank) {
  old_col <- paste0(rank, "_Fungi")
  new_col <- paste0(rank, "_All")
  comparison |>
    mutate(
      Rank = rank,
      Old_Label = coalesce(as.character(.data[[old_col]]), "Unassigned"),
      New_Label = coalesce(as.character(.data[[new_col]]), "Unassigned")
    ) |>
    filter(Old_Label != New_Label) |>
    group_by(Method, Rank, Old_Label, New_Label) |>
    summarise(ASVs = n(), Reads = sum(Total_Reads), .groups = "drop")
}) |>
  group_by(Method, Rank) |>
  mutate(Read_Fraction_Within_Rank = Reads / sum(Reads)) |>
  arrange(Method, Rank, desc(Reads)) |>
  ungroup()

focus_unresolved <- comparison |>
  mutate(
    Focus_Group = case_when(
      Kingdom_Group_Fungi == "Fungi" & Deepest_Rank_Fungi == "Kingdom" ~ "Fungi_k",
      Phylum_Fungi == "Ascomycota" & Deepest_Rank_Fungi == "Phylum" ~ "Ascomycota_p",
      TRUE ~ NA_character_
    ),
    AllDB_Genus = case_when(
      Informative_Genus_All ~ Genus_All,
      TRUE ~ "Unclassified genus"
    )
  ) |>
  filter(!is.na(Focus_Group))

focus_summary <- focus_unresolved |>
  group_by(Method, Focus_Group, AllDB_Genus) |>
  summarise(ASVs = n(), Reads = sum(Total_Reads), .groups = "drop") |>
  group_by(Method, Focus_Group) |>
  mutate(ASV_Fraction = ASVs / sum(ASVs), Read_Fraction = Reads / sum(Reads)) |>
  arrange(Method, Focus_Group, desc(Reads)) |>
  ungroup()

top_reclassified <- comparison |>
  filter(Kingdom_Group_All != "Fungi" | is.na(Kingdom_Group_All)) |>
  arrange(Method, desc(Total_Reads)) |>
  group_by(Method) |>
  slice_head(n = 30) |>
  ungroup() |>
  transmute(
    Method, FeatureID = as.character(FeatureID), ASV = as.character(substr(FeatureID, 1, 10)), Length, Total_Reads, Prevalence,
    Fungi_DB = Display_Taxon_Fungi, All_euk_DB = Display_Taxon_All,
    All_euk_Kingdom = coalesce(Kingdom_Group_All, "Unassigned"),
    Confidence_Fungi, Confidence_All
  )

known_raphanus <- comparison |>
  filter(FeatureID == "5c8f111a2f4a8e9e41a1cb58ebc985f6") |>
  transmute(Method, FeatureID, Length, Total_Reads, Prevalence,
            Fungi_DB = Raw_Taxon_Fungi, All_euk_DB = Raw_Taxon_All,
            Confidence_Fungi, Confidence_All)

sample_long <- imap_dfr(res, function(x, method) {
  common <- intersect(colnames(x$counts), metadata$SampleID)
  tax <- x$taxonomy |> select(FeatureID, DB, Kingdom_Group)
  map_dfr(db_levels, function(db) {
    keep <- tax |> filter(DB == db, Kingdom_Group == "Fungi") |> pull(FeatureID)
    keep <- intersect(keep, rownames(x$counts))
    fungal <- colSums(x$counts[keep, common, drop = FALSE])
    tibble(
      SampleID = common, Method = method, DB = db,
      Total_Reads = colSums(x$counts[, common, drop = FALSE]),
      Fungal_Reads = as.numeric(fungal),
      Fungal_ASVs = colSums(x$counts[keep, common, drop = FALSE] > 0)
    )
  })
}) |>
  left_join(metadata, by = "SampleID") |>
  mutate(
    Method = factor(Method, levels = method_levels), DB = factor(DB, levels = db_levels),
    Fungal_Fraction = if_else(Total_Reads > 0, Fungal_Reads / Total_Reads, NA_real_)
  )

sample_wide <- sample_long |>
  select(SampleID, Method, DB, Fungal_Reads, Fungal_ASVs, Fungal_Fraction) |>
  pivot_wider(names_from = DB, values_from = c(Fungal_Reads, Fungal_ASVs, Fungal_Fraction), names_sep = "__")

sample_loss <- sample_wide |>
  mutate(
    Lost_Fungal_Reads = .data[["Fungal_Reads__UNITE 99% Fungi"]] - .data[["Fungal_Reads__UNITE Dynamic All-euk"]],
    Lost_Fungal_ASVs = .data[["Fungal_ASVs__UNITE 99% Fungi"]] - .data[["Fungal_ASVs__UNITE Dynamic All-euk"]],
    Lost_Read_Fraction = if_else(.data[["Fungal_Reads__UNITE 99% Fungi"]] > 0,
                                Lost_Fungal_Reads / .data[["Fungal_Reads__UNITE 99% Fungi"]], NA_real_)
  ) |>
  left_join(metadata |> select(SampleID, Sample_Type, Subject_Group, Site_Group, PCR_Cycle_Group), by = "SampleID") |>
  arrange(Method, desc(Lost_Read_Fraction), desc(Lost_Fungal_Reads))

sample_kingdom_long <- imap_dfr(res, function(x, method) {
  common <- intersect(colnames(x$counts), metadata$SampleID)
  map_dfr(db_levels, function(db) {
    tx <- x$taxonomy |>
      filter(DB == db) |>
      transmute(FeatureID, Taxonomy = coalesce(Kingdom_Group, "Unassigned"))
    ids <- intersect(rownames(x$counts), tx$FeatureID)
    tax_vec <- tx$Taxonomy[match(ids, tx$FeatureID)]
    long <- as.data.frame(x$counts[ids, common, drop = FALSE]) |>
      rownames_to_column("FeatureID") |>
      mutate(Taxonomy = tax_vec) |>
      pivot_longer(-c(FeatureID, Taxonomy), names_to = "SampleID", values_to = "Reads") |>
      group_by(SampleID, Taxonomy) |>
      summarise(Reads = sum(Reads), .groups = "drop")
    long |> mutate(Method = method, DB = db)
  })
}) |>
  group_by(Method, DB, SampleID) |>
  mutate(Relative_Abundance = Reads / sum(Reads)) |>
  ungroup() |>
  left_join(metadata |> select(SampleID, Sample_Type, Subject_Group, Site_Group, PCR_Cycle_Group), by = "SampleID") |>
  mutate(Method = factor(Method, levels = method_levels), DB = factor(DB, levels = db_levels))

overall_kingdom_db <- sample_kingdom_long |>
  group_by(Method, DB, Taxonomy) |>
  summarise(Raw_Reads = sum(Reads), .groups = "drop") |>
  group_by(Method, DB) |>
  mutate(Relative_Abundance = Raw_Reads / sum(Raw_Reads)) |>
  ungroup()

# Whole-sample phylum composition before and after changing the taxonomy DB.
# Phyla with mean relative abundance >1% in at least one Method x DB group are
# retained in the legend; the remainder are combined as Other.
sample_phylum_raw <- imap_dfr(res, function(x, method) {
  common <- intersect(colnames(x$counts), metadata$SampleID)
  map_dfr(db_levels, function(db) {
    tx <- x$taxonomy |>
      filter(DB == db) |>
      transmute(
        FeatureID,
        Taxonomy = if_else(Informative_Phylum, Phylum, "Unclassified phylum")
      )
    ids <- intersect(rownames(x$counts), tx$FeatureID)
    tax_vec <- tx$Taxonomy[match(ids, tx$FeatureID)]
    as.data.frame(x$counts[ids, common, drop = FALSE]) |>
      rownames_to_column("FeatureID") |>
      mutate(Taxonomy = tax_vec) |>
      pivot_longer(-c(FeatureID, Taxonomy), names_to = "SampleID", values_to = "Reads") |>
      group_by(SampleID, Taxonomy) |>
      summarise(Reads = sum(Reads), .groups = "drop") |>
      group_by(SampleID) |>
      mutate(Relative_Abundance = Reads / sum(Reads)) |>
      ungroup() |>
      mutate(Method = method, DB = db)
  })
})

phylum_keep <- sample_phylum_raw |>
  group_by(Method, DB, Taxonomy) |>
  summarise(Mean_Relative_Abundance = mean(Relative_Abundance), .groups = "drop") |>
  filter(Mean_Relative_Abundance > 0.01 | Taxonomy == "Unclassified phylum") |>
  pull(Taxonomy) |>
  unique()

sample_order <- metadata |>
  arrange(Subject_Group, PCR_Cycle_Group, Site_Group, SampleID) |>
  pull(SampleID)

sample_phylum_long <- sample_phylum_raw |>
  mutate(Taxonomy = if_else(Taxonomy %in% phylum_keep, Taxonomy, "Others")) |>
  group_by(Method, DB, SampleID, Taxonomy) |>
  summarise(Reads = sum(Reads), .groups = "drop") |>
  group_by(Method, DB, SampleID) |>
  mutate(Relative_Abundance = Reads / sum(Reads)) |>
  ungroup() |>
  left_join(metadata |> select(SampleID, Sample_Type, Subject_Group, Site_Group, PCR_Cycle_Group), by = "SampleID") |>
  mutate(
    Method = factor(Method, levels = method_levels),
    DB = factor(DB, levels = db_levels),
    Sample_Type = factor(Sample_Type, levels = c("True", "Control", "Mock")),
    SampleID = factor(SampleID, levels = sample_order)
  )

# Fungal-only genus composition. The denominator is the sum of reads assigned
# to Kingdom Fungi within each sample and database.
sample_fungal_genus_raw <- imap_dfr(res, function(x, method) {
  common <- intersect(colnames(x$counts), metadata$SampleID)
  map_dfr(db_levels, function(db) {
    tx <- x$taxonomy |>
      filter(DB == db, Kingdom_Group == "Fungi") |>
      transmute(
        FeatureID,
        Taxonomy = if_else(Informative_Genus, Genus, "Unclassified genus")
      )
    ids <- intersect(rownames(x$counts), tx$FeatureID)
    if (!length(ids)) return(tibble())
    tax_vec <- tx$Taxonomy[match(ids, tx$FeatureID)]
    as.data.frame(x$counts[ids, common, drop = FALSE]) |>
      rownames_to_column("FeatureID") |>
      mutate(Taxonomy = tax_vec) |>
      pivot_longer(-c(FeatureID, Taxonomy), names_to = "SampleID", values_to = "Reads") |>
      group_by(SampleID, Taxonomy) |>
      summarise(Reads = sum(Reads), .groups = "drop") |>
      group_by(SampleID) |>
      mutate(Fungal_Reads = sum(Reads)) |>
      filter(Fungal_Reads > 0) |>
      mutate(Relative_Abundance = Reads / Fungal_Reads) |>
      ungroup() |>
      mutate(Method = method, DB = db)
  })
})

# Match the 01f composition rule: one shared top-14 list across every plotted
# method/database panel, based on mean sample-level relative abundance.
fungal_genus_keep <- sample_fungal_genus_raw |>
  group_by(Taxonomy) |>
  summarise(Mean_Relative_Abundance = mean(Relative_Abundance), .groups = "drop") |>
  arrange(desc(Mean_Relative_Abundance)) |>
  slice_head(n = 14) |>
  pull(Taxonomy)

sample_fungal_genus_long <- sample_fungal_genus_raw |>
  mutate(Taxonomy = if_else(Taxonomy %in% fungal_genus_keep, Taxonomy, "Others")) |>
  group_by(Method, DB, SampleID, Taxonomy) |>
  summarise(Reads = sum(Reads), .groups = "drop") |>
  group_by(Method, DB, SampleID) |>
  mutate(Relative_Abundance = Reads / sum(Reads)) |>
  ungroup() |>
  left_join(metadata |> select(SampleID, Sample_Type, Subject_Group, Site_Group, PCR_Cycle_Group), by = "SampleID") |>
  mutate(
    Method = factor(Method, levels = method_levels),
    DB = factor(DB, levels = db_levels),
    SampleID = factor(SampleID, levels = sample_order)
  )

paired_tests <- sample_wide |>
  group_by(Method) |>
  summarise(
    N = n(),
    Median_Reads_FungiDB = median(.data[["Fungal_Reads__UNITE 99% Fungi"]]),
    Median_Reads_AllDB = median(.data[["Fungal_Reads__UNITE Dynamic All-euk"]]),
    Median_ASVs_FungiDB = median(.data[["Fungal_ASVs__UNITE 99% Fungi"]]),
    Median_ASVs_AllDB = median(.data[["Fungal_ASVs__UNITE Dynamic All-euk"]]),
    Read_P = wilcox.test(.data[["Fungal_Reads__UNITE Dynamic All-euk"]], .data[["Fungal_Reads__UNITE 99% Fungi"]], paired = TRUE, exact = FALSE)$p.value,
    ASV_P = wilcox.test(.data[["Fungal_ASVs__UNITE Dynamic All-euk"]], .data[["Fungal_ASVs__UNITE 99% Fungi"]], paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  )

ksy_ids <- metadata |> filter(Sample_Type == "True", str_detect(Subject, regex("KSY", TRUE))) |> pull(SampleID)
ksy_fungi_k <- imap_dfr(res, function(x, method) {
  old <- x$taxonomy |> filter(DB == "UNITE 99% Fungi", Kingdom_Group == "Fungi", Deepest_Rank == "Kingdom")
  ids <- intersect(old$FeatureID, rownames(x$counts))
  if (!length(ids)) return(tibble())
  all_tax <- x$taxonomy |> filter(DB == "UNITE Dynamic All-euk") |> select(FeatureID, All_Kingdom = Kingdom_Group, All_Taxon = Display_Taxon, All_Confidence = Confidence)
  tibble(
    FeatureID = ids,
    KSY_True_Reads = rowSums(x$counts[ids, intersect(ksy_ids, colnames(x$counts)), drop = FALSE]),
    Total_Reads = rowSums(x$counts[ids, , drop = FALSE]),
    Length = x$feature$Length[match(ids, x$feature$FeatureID)]
  ) |>
    left_join(all_tax, by = "FeatureID") |>
    mutate(Method = method)
}) |>
  arrange(Method, desc(KSY_True_Reads)) |>
  transmute(Method, FeatureID = as.character(FeatureID), ASV = as.character(substr(FeatureID, 1, 10)), Length,
            KSY_True_Reads, Total_Reads, All_Kingdom, All_Taxon, All_Confidence)

# Interactive True-sample taxonomy flow: Subject -> fungi-only phylum -> all-euk phylum.
true_meta <- metadata |>
  filter(Sample_Type == "True") |>
  select(SampleID, Subject_Group)

true_flow_data <- imap_dfr(res, function(x, method) {
  method_tax <- x$comparison |>
    transmute(
      FeatureID,
      Old_Phylum = case_when(
        Informative_Phylum_Fungi ~ Phylum_Fungi,
        Kingdom_Group_Fungi == "Fungi" ~ "Fungi_k",
        !is.na(Kingdom_Group_Fungi) ~ paste0(Kingdom_Group_Fungi, "_k"),
        TRUE ~ "Unclassified"
      ),
      New_Phylum = case_when(
        Informative_Phylum_All ~ Phylum_All,
        !is.na(Kingdom_Group_All) & Kingdom_Group_All != "Unassigned" ~ paste0(Kingdom_Group_All, "_k"),
        TRUE ~ "Unclassified"
      )
    )
  map_dfr(unique(true_meta$Subject_Group), function(subject_name) {
    ids <- true_meta |> filter(Subject_Group == subject_name) |> pull(SampleID)
    ids <- intersect(ids, colnames(x$counts))
    if (!length(ids)) return(tibble())
    tibble(
      FeatureID = rownames(x$counts),
      Subject = subject_name,
      Reads = as.numeric(rowSums(x$counts[, ids, drop = FALSE]))
    )
  }) |>
    filter(Reads > 0) |>
    left_join(method_tax, by = "FeatureID") |>
    mutate(Method = method) |>
    group_by(Method, Subject, Old_Phylum, New_Phylum) |>
    summarise(Reads = sum(Reads), .groups = "drop") |>
    group_by(Method) |>
    mutate(Fraction_of_True_Reads = Reads / sum(Reads)) |>
    ungroup()
})

make_true_flow <- function(method_name) {
  z <- true_flow_data |> filter(Method == method_name)
  first_links <- z |>
    group_by(Subject, Old_Phylum) |>
    summarise(Reads = sum(Reads), .groups = "drop") |>
    transmute(
      Source = paste0("Subject: ", Subject),
      Target = paste0("Fungi DB: ", Old_Phylum), Reads
    )
  second_links <- z |>
    group_by(Old_Phylum, New_Phylum) |>
    summarise(Reads = sum(Reads), .groups = "drop") |>
    transmute(
      Source = paste0("Fungi DB: ", Old_Phylum),
      Target = paste0("All-euk DB: ", New_Phylum), Reads
    )
  links <- bind_rows(first_links, second_links)
  total_true <- sum(first_links$Reads)
  nodes <- tibble(Name = unique(c(links$Source, links$Target))) |>
    mutate(
      Node_Label = str_replace(Name, ": ", ":<br>"),
      Colour = case_when(
        str_starts(Name, "Subject:") ~ "#FDBF6F",
        str_starts(Name, "Fungi DB:") ~ "#1F78B4",
        TRUE ~ "#33A02C"
      )
    )
  links <- links |>
    mutate(
      source = match(Source, nodes$Name) - 1L,
      target = match(Target, nodes$Name) - 1L,
      Fraction = Reads / total_true,
      Hover = paste0(
        Source, " → ", Target, "<br>",
        scales::comma(Reads), " reads<br>",
        scales::percent(Fraction, accuracy = 0.01), " of True-sample reads"
      )
    )
  plotly::plot_ly(
    type = "sankey", orientation = "h",
    node = list(
      label = nodes$Node_Label, color = nodes$Colour,
      pad = 18, thickness = 22,
      line = list(color = "grey30", width = 0.4),
      hovertemplate = "%{label}<extra></extra>"
    ),
    link = list(
      source = links$source, target = links$target, value = links$Reads,
      customdata = links$Hover,
      hovertemplate = "%{customdata}<extra></extra>",
      color = "rgba(120,120,120,0.32)"
    )
  ) |>
    plotly::layout(
      title = list(text = paste0(method_name, ": True samples taxonomy flow")),
      font = list(size = 12), height = 650,
      margin = list(l = 20, r = 20, t = 55, b = 20)
    )
}

true_flow_widgets <- setNames(lapply(method_levels, make_true_flow), method_levels)

# Family-level version of the True-sample flow. Retain the 20 most abundant
# families on each DB side within a method and combine the remainder.
true_family_flow_raw <- imap_dfr(res, function(x, method) {
  method_tax <- x$comparison |>
    transmute(
      FeatureID,
      Old_Family = if_else(Informative_Family_Fungi, Family_Fungi, "Unclassified family"),
      New_Family = if_else(Informative_Family_All, Family_All, "Unclassified family")
    )
  map_dfr(unique(true_meta$Subject_Group), function(subject_name) {
    ids <- true_meta |> filter(Subject_Group == subject_name) |> pull(SampleID)
    ids <- intersect(ids, colnames(x$counts))
    if (!length(ids)) return(tibble())
    tibble(
      FeatureID = rownames(x$counts),
      Subject = subject_name,
      Reads = as.numeric(rowSums(x$counts[, ids, drop = FALSE]))
    )
  }) |>
    filter(Reads > 0) |>
    left_join(method_tax, by = "FeatureID") |>
    mutate(Method = method) |>
    group_by(Method, Subject, Old_Family, New_Family) |>
    summarise(Reads = sum(Reads), .groups = "drop")
})

true_family_flow_data <- true_family_flow_raw |>
  group_by(Method) |>
  group_modify(~{
    old_keep <- .x |>
      group_by(Old_Family) |>
      summarise(Reads = sum(Reads), .groups = "drop") |>
      slice_max(Reads, n = 20, with_ties = FALSE) |>
      pull(Old_Family)
    new_keep <- .x |>
      group_by(New_Family) |>
      summarise(Reads = sum(Reads), .groups = "drop") |>
      slice_max(Reads, n = 20, with_ties = FALSE) |>
      pull(New_Family)
    .x |>
      mutate(
        Old_Family = if_else(Old_Family %in% old_keep, Old_Family, "Other families"),
        New_Family = if_else(New_Family %in% new_keep, New_Family, "Other families")
      ) |>
      group_by(Subject, Old_Family, New_Family) |>
      summarise(Reads = sum(Reads), .groups = "drop")
  }) |>
  group_by(Method) |>
  mutate(Fraction_of_True_Reads = Reads / sum(Reads)) |>
  ungroup()

make_true_family_flow <- function(method_name) {
  z <- true_family_flow_data |> filter(Method == method_name)
  first_links <- z |>
    group_by(Subject, Old_Family) |>
    summarise(Reads = sum(Reads), .groups = "drop") |>
    transmute(
      Source = paste0("Subject: ", Subject),
      Target = paste0("Fungi DB family: ", Old_Family), Reads
    )
  second_links <- z |>
    group_by(Old_Family, New_Family) |>
    summarise(Reads = sum(Reads), .groups = "drop") |>
    transmute(
      Source = paste0("Fungi DB family: ", Old_Family),
      Target = paste0("All-euk DB family: ", New_Family), Reads
    )
  links <- bind_rows(first_links, second_links)
  total_true <- sum(first_links$Reads)
  nodes <- tibble(Name = unique(c(links$Source, links$Target))) |>
    mutate(
      Node_Label = str_replace(Name, ": ", ":<br>"),
      Colour = case_when(
        str_starts(Name, "Subject:") ~ "#FDBF6F",
        str_starts(Name, "Fungi DB family:") ~ "#1F78B4",
        TRUE ~ "#33A02C"
      )
    )
  links <- links |>
    mutate(
      source = match(Source, nodes$Name) - 1L,
      target = match(Target, nodes$Name) - 1L,
      Fraction = Reads / total_true,
      Hover = paste0(
        Source, " → ", Target, "<br>",
        scales::comma(Reads), " reads<br>",
        scales::percent(Fraction, accuracy = 0.01), " of True-sample reads"
      )
    )
  plotly::plot_ly(
    type = "sankey", orientation = "h",
    node = list(
      label = nodes$Node_Label, color = nodes$Colour,
      pad = 12, thickness = 18,
      line = list(color = "grey30", width = 0.35),
      hovertemplate = "%{label}<extra></extra>"
    ),
    link = list(
      source = links$source, target = links$target, value = links$Reads,
      customdata = links$Hover,
      hovertemplate = "%{customdata}<extra></extra>",
      color = "rgba(120,120,120,0.28)"
    )
  ) |>
    plotly::layout(
      title = list(text = paste0(method_name, ": True samples family-level taxonomy flow")),
      font = list(size = 11), height = 850,
      margin = list(l = 20, r = 20, t = 55, b = 20)
    )
}

true_family_flow_widgets <- setNames(lapply(method_levels, make_true_family_flow), method_levels)

kingdom_palette <- c(
  "Fungi" = "#1B9E77", "Viridiplantae" = "#66A61E", "Metazoa" = "#E7298A",
  "Alveolata" = "#7570B3", "Stramenopiles" = "#E6AB02", "Unassigned" = "grey70"
)
missing_kingdoms <- setdiff(unique(kingdom_summary$AllDB_Kingdom), names(kingdom_palette))
if (length(missing_kingdoms)) {
  extra <- colorRampPalette(brewer.pal(8, "Dark2"))(length(missing_kingdoms))
  kingdom_palette <- c(kingdom_palette, setNames(extra, missing_kingdoms))
}

p_kingdom_relative <- overall_kingdom_db |>
  ggplot(aes(DB, Relative_Abundance, fill = Taxonomy)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  facet_wrap(~Method, nrow = 1) +
  scale_fill_manual(values = kingdom_palette) +
  scale_y_continuous(labels = percent) +
  labs(x = "Database", y = "Relative abundance", fill = "Kingdom",
       title = "Database별 전체 kingdom relative composition") + soyeon_theme

p_rank_change_asv <- rank_change_summary |>
  ggplot(aes(Rank, ASV_Fraction, fill = Change_Status)) +
  geom_col(width = 0.68, colour = "white") +
  facet_wrap(~Method, nrow = 1) +
  scale_fill_manual(values = c(Same = "grey75", Changed = "#E31A1C")) +
  scale_y_continuous(labels = percent) +
  labs(x = "Taxonomic rank", y = "Fraction of ASVs", fill = NULL,
       title = "DB 변경에 따라 label이 달라진 ASV 비율") + soyeon_theme

p_rank_change_reads <- rank_change_summary |>
  ggplot(aes(Rank, Read_Fraction, fill = Change_Status)) +
  geom_col(width = 0.68, colour = "white") +
  facet_wrap(~Method, nrow = 1) +
  scale_fill_manual(values = c(Same = "grey75", Changed = "#E31A1C")) +
  scale_y_continuous(labels = percent) +
  labs(x = "Taxonomic rank", y = "Fraction of reads", fill = NULL,
       title = "DB 변경에 따라 label이 달라진 read-weighted 비율") + soyeon_theme

focus_taxa <- unique(focus_summary$AllDB_Genus)
focus_cols <- setNames(colorRampPalette(c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2")))(length(focus_taxa)), focus_taxa)
p_focus <- focus_summary |>
  ggplot(aes(Method, Read_Fraction, fill = AllDB_Genus)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  facet_wrap(~Focus_Group, nrow = 1) +
  scale_fill_manual(values = focus_cols) +
  scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "Read-weighted fraction", fill = "All-euk genus",
       title = "기존 Fungi_k / Ascomycota_p ASV의 All-euk genus 분류") + soyeon_theme

phylum_taxa <- sort(unique(sample_phylum_long$Taxonomy))
phylum_base <- setdiff(phylum_taxa, c("Others", "Unclassified phylum"))
paired_cols <- brewer.pal(12, "Paired")
phylum_cols <- if (length(phylum_base)) setNames(paired_cols[seq_along(phylum_base)], phylum_base) else character()
phylum_cols <- c(phylum_cols, "Unclassified phylum" = "grey45", "Others" = "grey82")
phylum_levels <- c(phylum_base, "Unclassified phylum", "Others")
sample_phylum_long <- sample_phylum_long |>
  mutate(Taxonomy = factor(Taxonomy, levels = rev(phylum_levels)))

draw_sample_phylum <- function(method_name, db_name) {
  sample_phylum_long |>
    filter(Method == method_name, DB == db_name) |>
    ggplot(aes(SampleID, Relative_Abundance, fill = Taxonomy)) +
    geom_col(width = 0.9) +
    facet_wrap(vars(Subject_Group, PCR_Cycle_Group), scales = "free_x", nrow = 1, drop = TRUE) +
    scale_fill_manual(values = phylum_cols, drop = FALSE) +
    scale_y_continuous(labels = percent, expand = c(0, 0)) +
    labs(
      x = "SampleID", y = "Relative abundance", fill = "Phylum",
      title = paste0(method_name, " — ", db_name, " sample-level Phylum composition"),
      subtitle = "Mean relative abundance >1% phyla are shown separately; remaining phyla are Others"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), legend.position = "right")
}
p_phylum_cutadapt_fungi <- draw_sample_phylum("Cutadapt Single", "UNITE 99% Fungi")
p_phylum_cutadapt_all <- draw_sample_phylum("Cutadapt Single", "UNITE Dynamic All-euk")
p_phylum_itsxpress_fungi <- draw_sample_phylum("ITSxpress Single", "UNITE 99% Fungi")
p_phylum_itsxpress_all <- draw_sample_phylum("ITSxpress Single", "UNITE Dynamic All-euk")

fungal_genus_base <- fungal_genus_keep
fungal_genus_levels <- c(fungal_genus_base, "Others")
fungal_genus_palette <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
fungal_genus_cols <- setNames(fungal_genus_palette[seq_along(fungal_genus_levels)], fungal_genus_levels)
sample_fungal_genus_long <- sample_fungal_genus_long |>
  mutate(Taxonomy = factor(Taxonomy, levels = fungal_genus_levels))

draw_sample_fungal_genus <- function(method_name, db_name) {
  sample_fungal_genus_long |>
    filter(Method == method_name, DB == db_name) |>
    ggplot(aes(SampleID, Relative_Abundance, fill = Taxonomy)) +
    geom_col(width = 0.9) +
    facet_wrap(vars(Subject_Group, PCR_Cycle_Group), scales = "free_x", nrow = 1, drop = TRUE) +
    scale_fill_manual(
      values = fungal_genus_cols,
      limits = fungal_genus_levels,
      breaks = fungal_genus_levels,
      drop = FALSE
    ) +
    scale_y_continuous(labels = percent, expand = c(0, 0)) +
    labs(
      x = "SampleID", y = "Relative abundance within Fungi", fill = "Fungal genus",
      title = paste0(method_name, " — ", db_name, " fungal Genus composition"),
      subtitle = "Kingdom == Fungi only; shared top 14 genera across all panels"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), legend.position = "right")
}
p_fungal_genus_cutadapt_fungi <- draw_sample_fungal_genus("Cutadapt Single", "UNITE 99% Fungi")
p_fungal_genus_cutadapt_all <- draw_sample_fungal_genus("Cutadapt Single", "UNITE Dynamic All-euk")
p_fungal_genus_itsxpress_fungi <- draw_sample_fungal_genus("ITSxpress Single", "UNITE 99% Fungi")
p_fungal_genus_itsxpress_all <- draw_sample_fungal_genus("ITSxpress Single", "UNITE Dynamic All-euk")

# Species-level beta diversity between taxonomy databases. Species labels use
# the same deepest-rank fallback convention as 01f, after Kingdom == Fungi.
species_label <- function(tx) {
  case_when(
    tx$Informative_Species ~ str_replace_all(tx$Species, "_", " "),
    tx$Informative_Genus ~ paste0(tx$Genus, "_g"),
    tx$Informative_Family ~ paste0(tx$Family, "_f"),
    tx$Informative_Order ~ paste0(tx$Order, "_o"),
    tx$Informative_Class ~ paste0(tx$Class, "_c"),
    tx$Informative_Phylum ~ paste0(tx$Phylum, "_p"),
    tx$Kingdom_Group == "Fungi" ~ "Fungi_k",
    TRUE ~ "Unclassified fungi"
  )
}

make_species_matrix <- function(x, db) {
  tx <- x$taxonomy |> filter(DB == db, Kingdom_Group == "Fungi")
  tx$Species_Label <- species_label(tx)
  ids <- intersect(rownames(x$counts), tx$FeatureID)
  z <- x$counts[ids, , drop = FALSE]
  rowsum(z, group = tx$Species_Label[match(ids, tx$FeatureID)], reorder = FALSE)
}

shape_values <- c(
  "Air swab" = 0, "Clean bench swab" = 1, "Empty tube" = 2,
  "KSY" = 5, "PYH" = 6, "Mock" = 4
)

analyse_db_beta <- function(method_name, metric_name) {
  x <- res[[method_name]]
  mats <- setNames(lapply(db_levels, function(db) make_species_matrix(x, db)), db_levels)
  paired_samples <- Reduce(intersect, lapply(mats, colnames))
  paired_samples <- paired_samples[
    colSums(mats[[db_levels[1]]][, paired_samples, drop = FALSE]) > 0 &
      colSums(mats[[db_levels[2]]][, paired_samples, drop = FALSE]) > 0
  ]
  all_species <- Reduce(union, lapply(mats, rownames))
  aligned <- lapply(mats, function(m) {
    z <- matrix(0, nrow = length(paired_samples), ncol = length(all_species),
                dimnames = list(paired_samples, all_species))
    z[, rownames(m)] <- t(m[, paired_samples, drop = FALSE])
    z
  })
  relative <- lapply(aligned, function(m) sweep(m, 1, rowSums(m), "/"))
  combined <- do.call(rbind, lapply(names(relative), function(db) {
    z <- relative[[db]]
    rownames(z) <- paste(db, rownames(z), sep = "||")
    z
  }))
  md <- tibble(RowID = rownames(combined)) |>
    separate(RowID, c("DB", "SampleID"), sep = "\\|\\|", extra = "merge") |>
    left_join(metadata |> select(SampleID, Subject_Group, Sample_Type), by = "SampleID") |>
    mutate(DB = factor(DB, levels = db_levels)) |>
    data.frame()
  rownames(md) <- rownames(combined)
  d <- if (metric_name == "Bray-Curtis") {
    vegan::vegdist(combined, method = "bray")
  } else {
    vegan::vegdist(combined, method = "jaccard", binary = TRUE)
  }
  set.seed(42)
  ad <- vegan::adonis2(d ~ DB, data = md, permutations = 999, strata = md$SampleID)
  pc <- vegan::wcmdscale(d, k = 2, eig = TRUE, add = "lingoes")
  positive_eig <- pc$eig[pc$eig > 0]
  variance <- positive_eig / sum(positive_eig) * 100
  score <- as_tibble(pc$points, rownames = "RowID") |>
    setNames(c("RowID", "Axis1", "Axis2")) |>
    left_join(rownames_to_column(md, "RowID"), by = "RowID")
  paired_distance <- map_dbl(paired_samples, function(s) {
    a <- relative[[db_levels[1]]][s, , drop = FALSE]
    b <- relative[[db_levels[2]]][s, , drop = FALSE]
    if (metric_name == "Bray-Curtis") {
      as.numeric(vegan::vegdist(rbind(a, b), method = "bray"))
    } else {
      as.numeric(vegan::vegdist(rbind(a, b), method = "jaccard", binary = TRUE))
    }
  })
  main <- ggplot(score, aes(Axis1, Axis2, colour = DB, shape = Subject_Group)) +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_line(aes(group = SampleID), colour = "grey72", alpha = 0.55, linewidth = 0.3) +
    geom_point(size = 2.1, alpha = 0.82) +
    scale_colour_manual(values = db_cols) +
    scale_shape_manual(values = shape_values, drop = TRUE) +
    labs(
      title = paste0(metric_name, " — ", method_name),
      subtitle = "Species-level fungal composition; paired samples are connected",
      x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
      y = paste0("PCoA2 (", round(variance[2], 1), "%)"),
      colour = "Database", shape = "Subject / control"
    ) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  xd <- cowplot::axis_canvas(main, axis = "x") +
    geom_density(data = score, aes(Axis1, fill = DB, colour = DB), alpha = 0.2, linewidth = 0.35) +
    scale_fill_manual(values = db_cols) + scale_colour_manual(values = db_cols)
  yd <- cowplot::axis_canvas(main, axis = "y", coord_flip = TRUE) +
    geom_density(data = score, aes(Axis2, fill = DB, colour = DB), alpha = 0.2, linewidth = 0.35) +
    scale_fill_manual(values = db_cols) + scale_colour_manual(values = db_cols) + coord_flip()
  plot <- main |>
    cowplot::insert_xaxis_grob(xd, grid::unit(0.58, "in"), position = "top") |>
    cowplot::insert_yaxis_grob(yd, grid::unit(0.58, "in"), position = "right") |>
    cowplot::ggdraw()
  stat <- tibble(
    Method = method_name, Distance = metric_name, N_Paired_Samples = length(paired_samples),
    PERMANOVA_R2 = ad$R2[1], PERMANOVA_P = ad$`Pr(>F)`[1],
    Median_Paired_Distance = median(paired_distance),
    Q1_Paired_Distance = quantile(paired_distance, 0.25),
    Q3_Paired_Distance = quantile(paired_distance, 0.75)
  )
  list(plot = plot, stat = stat, paired_distance = tibble(
    Method = method_name, Distance = metric_name,
    SampleID = paired_samples, Paired_DB_Distance = paired_distance
  ))
}

db_beta_results <- map(c("Cutadapt Single", "ITSxpress Single"), function(method_name) {
  setNames(map(c("Bray-Curtis", "Jaccard"), ~analyse_db_beta(method_name, .x)), c("Bray-Curtis", "Jaccard"))
}) |>
  setNames(c("Cutadapt Single", "ITSxpress Single"))
db_beta_stats <- map_dfr(db_beta_results, ~map_dfr(.x, "stat")) |>
  mutate(PERMANOVA_P_BH = p.adjust(PERMANOVA_P, "BH"))
db_beta_paired_distance <- map_dfr(db_beta_results, ~map_dfr(.x, "paired_distance"))
p_db_beta_cutadapt <- cowplot::plot_grid(
  db_beta_results[["Cutadapt Single"]][["Bray-Curtis"]]$plot,
  db_beta_results[["Cutadapt Single"]][["Jaccard"]]$plot,
  ncol = 1
)
p_db_beta_itsxpress <- cowplot::plot_grid(
  db_beta_results[["ITSxpress Single"]][["Bray-Curtis"]]$plot,
  db_beta_results[["ITSxpress Single"]][["Jaccard"]]$plot,
  ncol = 1
)

p_sample <- sample_long |>
  filter(DB == "UNITE Dynamic All-euk") |>
  ggplot(aes(PCR_Cycle_Group, Fungal_Fraction, fill = Sample_Type)) +
  geom_boxplot(width = 0.68, outlier.shape = NA, alpha = 0.8) +
  geom_point(aes(shape = Subject_Group), position = position_jitter(width = 0.12), size = 1.8, alpha = 0.8) +
  facet_wrap(~Method, nrow = 1) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_fill_manual(values = c(True = "#1F78B4", Control = "#FF7F00", Mock = "#6A3D9A")) +
  labs(x = "PCR condition", y = "Reads classified as Fungi", fill = "Sample type", shape = "Subject / control",
       title = "Sample별 All-euk DB fungal read fraction") + soyeon_theme

p_paired <- sample_wide |>
  pivot_longer(cols = starts_with("Fungal_Reads__"), names_to = "DB", values_to = "Fungal_Reads") |>
  mutate(DB = str_remove(DB, "^Fungal_Reads__"), DB = factor(DB, levels = db_levels)) |>
  ggplot(aes(DB, Fungal_Reads, group = SampleID)) +
  geom_line(alpha = 0.22, colour = "grey45") +
  geom_point(aes(colour = DB), size = 1.8, alpha = 0.75) +
  facet_wrap(~Method, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = db_cols) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Fungal reads", colour = "Database",
       title = "동일 sample에서 DB 변경 후 유지되는 fungal reads") + soyeon_theme

p_conf <- comparison |>
  ggplot(aes(Confidence_Fungi, Confidence_All, colour = Kingdom_Group_All)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(size = Total_Reads), alpha = 0.58) +
  facet_wrap(~Method, nrow = 1) +
  scale_colour_manual(values = kingdom_palette) +
  scale_size_continuous(range = c(1, 7), trans = "sqrt", labels = comma) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "UNITE 99% Fungi confidence", y = "UNITE Dynamic All-euk confidence",
       colour = "All-euk kingdom", size = "Total reads",
       title = "ASV별 classifier confidence 변화") + soyeon_theme

p_loss <- sample_loss |>
  filter(Lost_Fungal_Reads > 0) |>
  mutate(Sample_Label = factor(
    paste(Method, SampleID, sep = "___"),
    levels = rev(unique(paste(Method, SampleID, sep = "___")))
  )) |>
  ggplot(aes(Sample_Label, Lost_Read_Fraction, fill = Sample_Type)) +
  geom_col(width = 0.75) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  coord_flip() +
  scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(True = "#1F78B4", Control = "#FF7F00", Mock = "#6A3D9A")) +
  labs(x = "SampleID", y = "Fungal read loss fraction", fill = "Sample type",
       title = "All-euk DB 적용 후 fungal read가 감소한 sample") + soyeon_theme +
  theme(axis.text.x = element_text(angle = 0))

top_loss_ids <- sample_loss |>
  filter(Lost_Fungal_Reads > 0) |>
  group_by(Method) |>
  slice_max(Lost_Read_Fraction, n = 12, with_ties = FALSE) |>
  ungroup() |>
  select(Method, SampleID)

loss_composition <- sample_kingdom_long |>
  inner_join(top_loss_ids, by = c("Method", "SampleID")) |>
  mutate(SampleID = factor(SampleID, levels = unique(top_loss_ids$SampleID)))

draw_loss_composition <- function(method) {
  loss_composition |>
    filter(Method == method) |>
    ggplot(aes(DB, Relative_Abundance, fill = Taxonomy)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    facet_wrap(~SampleID, ncol = 4) +
    scale_fill_manual(values = kingdom_palette) +
    scale_y_continuous(labels = percent) +
    labs(x = "Database", y = "Relative abundance", fill = "Kingdom",
         title = paste0(method, ": fungal read 감소 상위 sample의 DB별 kingdom composition")) +
    soyeon_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}
p_loss_cutadapt <- draw_loss_composition("Cutadapt Single")
p_loss_itsxpress <- draw_loss_composition("ITSxpress Single")

ggsave(file.path(out_dir, "01a_kingdom_relative_composition.png"), p_kingdom_relative, width = 11.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "02_sample_fungal_fraction.png"), p_sample, width = 14, height = 7.2, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "03_paired_fungal_reads.png"), p_paired, width = 11.5, height = 6.5, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "04_classifier_confidence.png"), p_conf, width = 12, height = 6.5, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "05a_rank_change_ASV_fraction.png"), p_rank_change_asv, width = 11.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "05b_rank_change_read_fraction.png"), p_rank_change_reads, width = 11.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "06_Fungi_k_Ascomycota_p_destinations.png"), p_focus, width = 13, height = 7, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "06a1_Cutadapt_UNITE99Fungi_sample_phylum.png"), p_phylum_cutadapt_fungi, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06a2_Cutadapt_AllEuk_sample_phylum.png"), p_phylum_cutadapt_all, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06b1_ITSxpress_UNITE99Fungi_sample_phylum.png"), p_phylum_itsxpress_fungi, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06b2_ITSxpress_AllEuk_sample_phylum.png"), p_phylum_itsxpress_all, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06c1_Cutadapt_UNITE99Fungi_fungal_genus.png"), p_fungal_genus_cutadapt_fungi, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06c2_Cutadapt_AllEuk_fungal_genus.png"), p_fungal_genus_cutadapt_all, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06d1_ITSxpress_UNITE99Fungi_fungal_genus.png"), p_fungal_genus_itsxpress_fungi, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06d2_ITSxpress_AllEuk_fungal_genus.png"), p_fungal_genus_itsxpress_all, width = 42, height = 8, dpi = 190, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "06e_species_Bray_Jaccard_PCoA_Cutadapt.png"), p_db_beta_cutadapt, width = 11.5, height = 14, dpi = 210, bg = "white")
ggsave(file.path(out_dir, "06f_species_Bray_Jaccard_PCoA_ITSxpress.png"), p_db_beta_itsxpress, width = 11.5, height = 14, dpi = 210, bg = "white")
ggsave(file.path(out_dir, "07_sample_fungal_read_loss.png"), p_loss, width = 11.5, height = 14, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "08a_Cutadapt_top_loss_sample_composition.png"), p_loss_cutadapt, width = 15, height = 10, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "08b_ITSxpress_top_loss_sample_composition.png"), p_loss_itsxpress, width = 15, height = 10, dpi = 220, bg = "white")

write_csv(summary_table, file.path(out_dir, "01_method_database_summary.csv"))
write_csv(kingdom_summary, file.path(out_dir, "02_all_euk_kingdom_summary.csv"))
write_csv(transition_summary, file.path(out_dir, "03_kingdom_transition_summary.csv"))
write_csv(rank_transition, file.path(out_dir, "04_rank_transition_summary.csv"))
write_csv(rank_change_summary, file.path(out_dir, "04b_rank_change_fraction.csv"))
write_csv(rank_transition_detail, file.path(out_dir, "04c_rank_transition_detail.csv"))
write_csv(focus_summary, file.path(out_dir, "04d_Fungi_k_Ascomycota_p_destination_summary.csv"))
write_csv(sample_phylum_long, file.path(out_dir, "04e_all_samples_phylum_relative_composition.csv"))
write_csv(sample_fungal_genus_long, file.path(out_dir, "04f_fungal_genus_relative_composition.csv"))
write_csv(db_beta_stats, file.path(out_dir, "04g_species_Bray_Jaccard_PCoA_statistics.csv"))
write_csv(db_beta_paired_distance, file.path(out_dir, "04h_species_paired_database_distances.csv"))
write_csv(top_reclassified, file.path(out_dir, "05_top_reclassified_ASVs.csv"))
write_csv(known_raphanus, file.path(out_dir, "06_Raphanus_ASV_database_comparison.csv"))
write_csv(sample_long, file.path(out_dir, "07_sample_fungal_retention.csv"))
write_csv(sample_loss, file.path(out_dir, "07b_sample_fungal_read_loss.csv"))
write_csv(sample_kingdom_long, file.path(out_dir, "07c_sample_database_kingdom_composition.csv"))
write_csv(paired_tests, file.path(out_dir, "08_paired_database_tests.csv"))
write_csv(ksy_fungi_k, file.path(out_dir, "09_KSY_Fungi_k_ASV_reclassification.csv"))
write_csv(true_flow_data, file.path(out_dir, "09b_all_true_sample_taxonomy_flow.csv"))
write_csv(true_family_flow_data, file.path(out_dir, "09c_all_true_sample_family_flow.csv"))
saveRDS(list(comparison = comparison, taxonomy_long = taxonomy_long, summary_table = summary_table,
             kingdom_summary = kingdom_summary, transition_summary = transition_summary,
             rank_transition = rank_transition, rank_change_summary = rank_change_summary,
             rank_transition_detail = rank_transition_detail, focus_summary = focus_summary,
             sample_phylum_long = sample_phylum_long,
             sample_fungal_genus_long = sample_fungal_genus_long,
             db_beta_stats = db_beta_stats, db_beta_paired_distance = db_beta_paired_distance,
             top_reclassified = top_reclassified,
             known_raphanus = known_raphanus, sample_long = sample_long,
             sample_wide = sample_wide, sample_loss = sample_loss,
             sample_kingdom_long = sample_kingdom_long, overall_kingdom_db = overall_kingdom_db,
             paired_tests = paired_tests,
             ksy_fungi_k = ksy_fungi_k, true_flow_data = true_flow_data,
             true_family_flow_data = true_family_flow_data),
        file.path(out_dir, "01g_database_comparison_objects.rds"))
