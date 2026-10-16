qza_path <- file.path(
  "benchmark_results", "q20_ee22_both_full", "qiime_preprocess",
  "itsxpress", "paired", "trimmed_itsxpress_pe.qza"
)

read_fastq_from_qza <- function(qza_path, entry) {
  con <- gzcon(unz(qza_path, entry, open = "rb"))
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  list(
    id = sub(" .*", "", lines[seq.int(1, length(lines), by = 4)]),
    sequence = lines[seq.int(2, length(lines), by = 4)],
    quality = lines[seq.int(4, length(lines), by = 4)]
  )
}

entries <- unzip(qza_path, list = TRUE)$Name
r1_entries <- entries[grepl("/data/.*_R1_001\\.fastq\\.gz$", entries)]
sample_ids <- sub("_[0-9]+_L001_R1_001\\.fastq\\.gz$", "", basename(r1_entries))
motif <- "CACCACTCAAGCTATGCTTGGTATTGGGTGCCGTCCTT"

diagnostic <- lapply(seq_along(r1_entries), function(i) {
  r1_entry <- r1_entries[[i]]
  r2_entry <- entries[
    startsWith(basename(entries), paste0(sample_ids[[i]], "_")) &
      grepl("_R2_001\\.fastq\\.gz$", entries)
  ][[1]]
  r1 <- read_fastq_from_qza(qza_path, r1_entry)
  r2 <- read_fastq_from_qza(qza_path, r2_entry)
  stopifnot(length(r1$sequence) == length(r2$sequence))

  target <- startsWith(r1$sequence, motif)
  r2_length <- nchar(r2$sequence)
  data.frame(
    SampleID = sample_ids[[i]],
    ITSxpress_Pairs = length(target),
    Tremulum_Motif_Pairs = sum(target),
    Tremulum_Motif_Fraction = mean(target),
    Tremulum_R2_lt20 = sum(target & r2_length < 20),
    Tremulum_R2_lt20_Fraction = if (sum(target) > 0) {
      mean(r2_length[target] < 20)
    } else {
      NA_real_
    },
    Tremulum_R1_Median_Length = if (sum(target) > 0) median(nchar(r1$sequence[target])) else NA_real_,
    Tremulum_R2_Median_Length = if (sum(target) > 0) median(r2_length[target]) else NA_real_
  )
})

diagnostic <- do.call(rbind, diagnostic)
sample_context <- read.csv(
  file.path("Results", "01_Preprocessing", "06q_aureobasidium_by_sample.csv"),
  check.names = FALSE
)
diagnostic <- merge(diagnostic, sample_context, by = "SampleID", all.x = TRUE)

write.csv(
  diagnostic,
  file.path("Results", "01_Preprocessing", "Aureobasidium_ITSxpress_pair_diagnostic.csv"),
  row.names = FALSE,
  na = ""
)

print(
  aggregate(
    cbind(Tremulum_Motif_Pairs, Tremulum_R2_lt20) ~ Subject,
    data = diagnostic,
    FUN = sum
  )
)
