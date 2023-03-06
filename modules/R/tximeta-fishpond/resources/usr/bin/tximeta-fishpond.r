#!/usr/bin/env Rscript
# Load libraries ---------------------------
library(tibble)
library(stringr)
library(purrr)
library(readr)
library(tximeta)
library(SummarizedExperiment)
library(Biostrings)
library(fishpond)

# Load samples using command line arguments ---------------------------
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No arguments were passed to the script.", call. = FALSE)
}

# Setup fishpond filtering for determining size factors

# Fetch each flag's position. If that flag is passed, fetch the pos + 1 item and
# use that as the value. Remove the items after to clean up the args
extract_arg <- function(.flag, .default, .args) {
  # Fetch the position of the flag (if there is one)
  flag_pos <- which(.args == .flag)

  # Make sure only one item matches at most
  stopifnot(length(flag_pos) <= 1)

  # IF nothing was passed, return the default
  if (length(flag_pos) == 0) {
    return(list(
      flag = .flag,
      value = .default,
      args = .args
    ))
  }

  # We have a flag passed. Make sure it isn't at the end (i.e., there is a value after)
  stopifnot((flag_pos == length(.args)) == FALSE)

  # Now fetch the item
  param_value <- .args[flag_pos + 1]

  # return the items
  list(
    flag = .flag,
    value = param_value,
    args = .args[-c(flag_pos, flag_pos + 1)]
  )
}

# Extract the min transcript (maybe gene) count
minCountFlag <- extract_arg("--minCount", NULL, args)
args = minCountFlag$args
min_count <- minCountFlag$value

# Extract the min gene count
minNFlag <- extract_arg("--minN", NULL, args)
args = minNFlag$args
min_n <- minNFlag$value

# Extract the min gene count
minCountGeneFlag <- extract_arg("--minCountGene", NULL, args)
args = minCountGeneFlag$args
min_count_gene <- minCountGeneFlag$value

# Extract the min gene number
minNGeneFlag <- extract_arg("--minNGene", NULL, args)
args = minNGeneFlag$args
min_n_gene <- minNGeneFlag$value

# Message user
message("Inferential replicates will be filtered based on the following criteria:")
message(paste("counts >=", min_count))
message(paste("n >=", min_n))
if (is.null(min_count_gene) == FALSE) message(paste("gene_counts >=", min_count_gene))
if (is.null(min_n_gene) == FALSE) message(paste("gene_n >=", min_n_gene))

## List of sample and Directories ---------------------------
# The script accepts four possible variations:
#    1) samples and the directory to find the samples. This option requires a
#       sample flag to be passed to the script.
#    2) Just a single directory. If this option is selected, all folders in the
#       directory will be samples and the sample id will be the folder name.
#    3) A list of directories. Similar to option 2, the directory names will be
#       used as the sample names
#    4) A list of samples and directories

samples_flag_pos <- which(args == "--samples")
directory_flag_pos <- which(args == "--directory")

if (length(samples_flag_pos) == 0 && length(directory_flag_pos) == 0) {
  # Update the user about what method was picked
  message("Using passed list of directories.")

  # Get the names of each folder (sample ids) and the sample's directory as a
  # named vector.
  samples <- basename(args)
  sample_paths <- args
}

## List of directories ---------------------------
if (length(samples_flag_pos) != 0 && length(directory_flag_pos) != 0) {
  # Update the user about what method was picked
  message("Using the passed list of samples and directories")

  # Find out which item is listed first: samples or directories
  # If samples are listed first:
  if (samples_flag_pos < directory_flag_pos) {
    samples <- args[((samples_flag_pos + 1):(directory_flag_pos - 1))]
    sample_paths <- args[((directory_flag_pos + 1):length(args))]
  }

  # If directories are listed first:
  if (directory_flag_pos < samples_flag_pos) {
    samples <- args[((samples_flag_pos + 1):length(args))]
    sample_paths <- args[((directory_flag_pos + 1):(samples_flag_pos - 1))]
  }
}

## Single directory ---------------------------
# Means we have a directory flag (first item) and no samples. That means we
# want to search one directory for the samples
if (args[1] == "--directory" && length(samples_flag_pos) == 0) {
  # Update the user about what method was picked
  message("Searching directory for samples.")

  # Find directory for targets and find get the names of the folders (samples) inside
  directory <- args[2]

  # Scan the directory fetch the names and assemble the path
  samples <- list.files(directory)
  sample_paths <- file.path(directory, samples)
}

## List of samples and a single directory ---------------------------
if (args[1] == "--samples" && length(directory_flag_pos) == 0) {
  # Update the user about what method was picked
  message("Using passed list of samples to search a target directory.")

  # Get the samples. The samples will be everything but the first item (the flag),
  # and the last item (directory)
  samples <- args[-c(1, length(args))]

  # Find directory for targets which is the last item in the arguments
  directory <- args[length(args)]

  # Assemble directories
  sample_paths <- file.path(directory, samples)
}

## Sample path determination ---------------------------
# Make sure we know what the samples are and their locations
# Do this by scanning the ls() items for a regex looking for the variable names.
# There should be only two hits. When converting the logics to integers, we should
# get a sum of two.
stopifnot(sum(str_detect(ls(), "^samples$|^sample_paths$")) == 2)

# Quality control step to make sure every sample has a path
stopifnot(length(samples) == length(sample_paths))

# Construct a vector pointing to their Salmon quant.sf files and name the vector
# according to the sample of origin. Use `fs::path_tidy` to clean path names
sample_paths <-
  fs::path_tidy(file.path(sample_paths, "quant.sf")) |>
  set_names(samples)

# Quality control check to make sure all files exist
stopifnot(all(file.exists(sample_paths)))

# Create a tibble with all the information to pass to `tximeta`
samples_table <- enframe(sample_paths, name = "names", value = "files")

# Update the user
message("Samples located")
message("Samples: ", str_flatten(samples, collapse = " "))
message("Paths: ",   str_flatten(sample_paths, collapse = " "))

# tximeta: import transcript level data ---------------------------
message("Importing transcripts")
se <- tximeta(samples_table)

# We want to normalize the counts, so we are going to label them as an inferential
# replicate so fishpond will normalize it. This is temporary and will be removed
# once the scaling is complete. This value will not affect normalization since
# each replicate is normalized independently of the others.
assay(se, "infRep0") <- assay(se, "counts")
se <- scaleInfReps(se, minCount = min_count, minN = min_n)

assay(se, "scaledTPMCounts") <- assay(se, "infRep0")
assay(se, "infRep0") <- NULL

# Now compute uncertainty
se <- computeInfRV(se)

# Save annotation database (sqlite) ---------------------------
# First assemble the file name using metadata information from tximeta import
db <- retrieveDb(se)
txomeInfo <- metadata(se)$txomeInfo

# Don't forget to clean up the names! Replace spaces with underscores, and make
# it all lowercase. Finally, sanitize the name to remove anything bad.
dbfile <-
  str_c(
    txomeInfo$source, txomeInfo$organism, txomeInfo$genome, txomeInfo$release, "sqlite",
    sep = "."
  ) |>
  str_replace_all("\\s+", "_") |>
  str_to_lower()
  fs::path_sanitize()

# Finally save a copy using the sqlite extension
message("Saving transcript database")
AnnotationDbi::saveDb(db, file = dbfile)

# Build transcript annotation ---------------------------
message("Building transcript annotation")

# Get the genome build for naming purposes
genome_build <-
  rtracklayer::ucscGenomes(TRUE) |>
  dplyr::filter(
    stringr::str_to_lower(organism) == str_to_lower(metadata(se)$txomeInfo$organism),
    stringr::str_detect(name, metadata(se)$txomeInfo$genome),
  ) |>
  dplyr::pull(db)

# Check if we got only one item. If not, fallback
if (length(genome_build) != 1) {
  genome_build <- metadata(se)$txomeInfo$genome |> str_replace_all("\\s+", "_") |> str_to_lower()
}

# Strings to name start and end columns
feature_start_genome = str_c("feature_start_", genome_build)
feature_end_genome = str_c("feature_end_", genome_build)

# Get transcript annotation
transcript_annotation <-
  db |>
  transcripts(
    columns = c("tx_chrom", "tx_start", "tx_end", "tx_strand", "tx_name", "gene_id")
  ) |>
  as_tibble() |>
  dplyr::select(tx_chrom:gene_id) |>
  dplyr::full_join(
    transcriptLengths(db) |> dplyr::select(tx_name, nexon, tx_len)
  ) |>
  dplyr::rename(
    {{ feature_start_genome }} := tx_start,
    {{ feature_end_genome }} := tx_end,
    transcript_id = tx_name,
    chromosome = tx_chrom,
    strand = tx_strand,
    transcript_length = tx_len,
    exon_count = nexon,
  ) |>
  tidyr::unnest_longer(gene_id)

# Get the GC content
message("Finding transcript GC content")
cdna <- retrieveCDNA(se)

transcripts_percent_gc <-
  cdna |> letterFrequency(letters = "GC", as.prob = TRUE) |>
  as.numeric() |>
  set_names(str_split_i(names(cdna), fixed("|"), i = 1)) |>
  enframe(name = "transcript_id", value = "transcript_percent_gc") |>
  dplyr::mutate(transcript_percent_gc = 100 * transcript_percent_gc)

# add to the annotation
transcript_annotation <-
  dplyr::full_join(transcript_annotation, transcripts_percent_gc)

# Build gene annotation ---------------------------
message("Building gene annotation")

# Find out how many transcripts each gene has
transcript_counts_by_gene <-
  transcript_annotation |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(transcript_count = dplyr::n())

gene_annotation <-
  genes(db) |>
  as_tibble() |>
  dplyr::rename(
    {{ feature_start_genome }} := start,
    {{ feature_end_genome }} := end,
    chromosome = seqnames,
  ) |>
  dplyr::full_join(transcript_counts_by_gene) |>
  dplyr::rename(
    gene_length = width
  )

# Type coercing for transcript data ---------------------------
message("Finishing up transcript annotation.")
transcript_annotation <-
  transcript_annotation |>
  dplyr::mutate(
    chromosome = factor(chromosome, levels(gene_annotation$chromosome)),
    strand = factor(strand, levels(gene_annotation$strand)),
  )

# Annotation metadata ---------------------------
# Retrieve metadata about the database and annotation we are using
annotation_metadata <- jsonlite::toJSON(metadata(se)$txomeInfo, pretty = TRUE)
annotation_source   <- jsonlite::toJSON(as.list(metadata(se)$txdbInfo), pretty = TRUE)
tximeta_info <-
  as_tibble(metadata(se)$tximeta) |>
  dplyr::mutate(
    version = as.character(version),
    type = as.character(type),
    import_time = as.character(importTime)
  ) |>
  dplyr::select(-importTime) |>
  as.list() |>
  jsonlite::toJSON(pretty = TRUE)

# Saving annotation data ---------------------------
message("Saving...")

readr::write_rds(cdna, "sequence/cdna.rds", compress = "gz")

readr::write_tsv(transcript_annotation, "annotation/tsv/transcripts.tsv.gz")
readr::write_tsv(gene_annotation, "annotation/tsv/genes.tsv.gz")

readr::write_rds(transcript_annotation, "annotation/rds/transcripts.rds", compress = "gz")
readr::write_rds(transcript_annotation, "annotation/rds/genes.rds", compress = "gz")

readr::write_file(annotation_metadata, "tximeta.json")
readr::write_file(annotation_metadata, "annotation/metadata.json")
readr::write_file(annotation_source, "annotation/source.json")

# Quantification summary ---------------------------
message("Gathering quantification results")

# Get quality data about how the alignment and quantification went
# We remove the `length_class` variable because it didn't work with tibbles
quant_info <-
  metadata(se)[["quantInfo"]] |>
  discard_at(c("quant_errors", "length_classes")) |>
  as_tibble() |>
  dplyr::transmute(
    sample_id = colnames(se),
    algorithm = opt_type,
    library_type = library_types,
    mean_fragment_length = frag_length_mean,
    sd_fragment_length = frag_length_sd,
    valid_targets = num_valid_targets,
    decoy_targets = num_decoy_targets,
    eq_classes = num_eq_classes,
    inf_replicates_type = samp_type,
    inf_replicates_number = num_bootstraps,
    fragments_total = num_processed,
    fragments_mapped = num_mapped,
    fragments_mapped_decoy = num_decoy_fragments,
    fragments_dovetailed = num_dovetail_fragments,
    fragments_low_quality = num_fragments_filtered_vm,
    mappings_low_quality = num_alignments_below_threshold_for_mapped_fragments_vm,
    percent_fragments_mapped = percent_mapped,
    percent_fragments_decoy = round(100 * fragments_mapped_decoy / fragments_total, digits = 5),
    percent_fragments_dovetailed = round(100 * fragments_dovetailed / fragments_total, digits = 5),
    percent_fragments_low_quality = round(100 * fragments_low_quality / fragments_total, digits = 5),
    percent_fragments_handled = round(
      100 * (fragments_mapped + fragments_mapped_decoy + fragments_dovetailed + fragments_low_quality) / fragments_total,
      digits = 5
    ),
  )

# Library format summary ---------------------------
message("Gathering quantification results")

lib_format_counts <-
  file.path(dirname(sample_paths), "lib_format_counts.json") |>
  set_names(samples) |>
  map(jsonlite::read_json) |>
  map(as_tibble) |>
  map2(samples, function(x, name) x |> dplyr::mutate(sample_id = name)) |>
  reduce(dplyr::full_join) |>
  dplyr::relocate(sample_id) |>
  dplyr::select(
    -read_files,
    -compatible_fragment_ratio
  ) |>
  dplyr::rename(
    fragments_compatible = num_compatible_fragments,
    fragments_assigned = num_assigned_fragments,
    fragments_with_consistent_mappings = num_frags_with_concordant_consistent_mappings,
    fragments_with_inconsistent_mappings = num_frags_with_inconsistent_or_orphan_mappings,
  )

# Fetch counting ambiguity ---------------------------
message("Calculating count ambiguity...")

read_mapping_ambiguity <-
  sample_paths |>
  set_names(samples) |>
  map(function(quant_file){
    # Read the results file
    path     <- file.path(dirname(quant_file), "aux_info", "ambig_info.tsv")
    ambig_df <- readr::read_tsv(path, progress = FALSE, show_col_types = FALSE)

    # Read from the first row the names of transcripts
    names <- readr::read_tsv(quant_file, col_select = "Name", progress = FALSE, show_col_types = FALSE)

    # Quality check before merging and returning
    stopifnot(nrow(names) == nrow(ambig_df))
    ambig_df |> dplyr::mutate(transcript_id = names)
  })

# Find the number of mapped reads that were uniquely assigned to a transcript
message("...for uniquely mapping reads")

uniquely_mapping_reads <-
  map2(
    read_mapping_ambiguity,
    names(read_mapping_ambiguity),
    function(df, sample) {
      df |> dplyr::select(transcript_id, UniqueCount) |> dplyr::rename({{ sample }} := UniqueCount)
    }
  ) |>
  reduce(dplyr::full_join)

# Find the number of reads that were assigned to multiple transcripts
message("...for ambiguously mapping reads")

ambiguously_mapping_reads <-
  map2(
    read_mapping_ambiguity,
    names(read_mapping_ambiguity),
    function(df, sample) {
      df |> dplyr::select(transcript_id, AmbigCount) |> dplyr::rename({{ sample }} := AmbigCount)
    }
  ) |>
  reduce(dplyr::full_join)

# Save quant summaries ---------------------------
message("Saving summaries")

readr::write_csv(quant_info, "quantification-summary.csv.gz")
readr::write_csv(lib_format_counts, "library-format.csv.gz")
readr::write_csv(uniquely_mapping_reads, "unique-mappings.csv.gz")
readr::write_csv(ambiguously_mapping_reads, "ambiguous-mappings.csv.gz")

# Summarizing counts to gene ---------------------------
message("Calculating gene counts")

# If filtering values aren't set, take the transcript values
if (is.null(min_count_gene)) min_count_gene <- min_count
if (is.null(min_n_gene)) min_n_gene <- min_n

# Fetch a new, unfiltered copy
gse <- summarizeToGene(tximeta(samples_table))

# Same logic as the transcripts is used for the genes
assay(gse, "infRep0") <- assay(gse, "counts")
gse <- scaleInfReps(gse, minCount = min_count_gene, minN = min_n_gene)

assay(gse, "scaledTPMCounts") <- assay(gse, "infRep0")
assay(gse, "infRep0") <- NULL

# Now compute uncertainty
gse <- computeInfRV(gse)

# Library format summary ---------------------------
message("Calculating transcript counts for differential transcript usage (DTU)")

# Fetch a new SE
se_dtu <- tximeta(coldata)
assay(se_dtu, "infRep0") <- assay(se_dtu, "counts")

# Length mean for each transcript
# Find the mean legnth of each transcript using all samples
length_geomean <- exp(rowMeans(log(assay(se_dtu, "length"))))

# DTU corrected median lengths, grouped by gene
# Group each transcript by gene and find the median length
median_length_geomean <-
  length_geomean |> enframe(name = "transcript_id", value = "length") |>
  dplyr::left_join(transcript_annotation |> dplyr::select(gene_id, transcript_id)) |>
  dplyr::mutate(transcript_id, length = median(length), .by = gene_id) |>
  dplyr::select(-gene_id) |>
  deframe()

# Actually correction phase. 
# 1) Scale the counts using fishpond
# 2) Find the correction factor produce by fishpond
# 3) Determine TPM using the original counts
# 4) Scale the TPM to counts using the correction factor and median lengths
assays(se_dtu) <-
  pmap(
    list(
      scaleInfReps(se_dtu, minCount = min_count, minN = min_n) |> assays(),
      assays(se_dtu),
      assayNames(se_dtu)
    ),
    function(scaledCounts, originalCounts, repName) {

      # Only operate on infreps
      if (str_detect(repName, "infRep") == FALSE) { return(scaledCounts) }

      # Calculate TPM from the scaled counts
      # scaledTPM <- (scaledCounts / length_geomean)
      # scaledAbundance <- 1e6 * t(t(scaledTPM) / colSums(scaledTPM))

      # Fetch size factors
      # libSize <- exp(mean(log( colSums(originalCounts)  )))
      # sf <- libSize / colSums(scaledCounts)

      # Fetch shared library size, corrected for by size factor
      library_size_correction <- colSums(scaledCounts)

      # convert the fragmentation rate based on our median values
      length_corrected_counts <- originalCounts / assay(se_dtu, "length")
      abundance <- t(t(length_corrected_counts) / colSums(length_corrected_counts))

      # TPM (abundance) to counts using the median lengths
      effective_fragmentation <- abundance * median_length_geomean
      t(t(effective_fragmentation) * (library_size_correction / colSums(effective_fragmentation)))
    }
  )

# Save the length
assay(se_dtu, "length") <- 
  matrix(
    rep(median_length_geomean, ncol(assay(se_dtu, "length"))),
    ncol = ncol(assay(se_dtu, "length")), 
    dimnames = dimnames(assay(se_dtu, "length"))
  )

# Save the coutns
assay(se_dtu, "dtuScaledTPMCounts") <- assay(se_dtu, "infRep0")
assay(se_dtu, "infRep0") <- NULL

# Now compute uncertainty
se_dtu <- computeInfRV(se_dtu)

# Save all the counts ---------------------------
message("Saving counts from all experiments")

# Helper function for saving summarized experiments
save_se_assay <- function(se, assay, .rownames, .file) {
  readr::write_csv(
    assay(se, assay) |> as_tibble(rownames = .rownames),
    .file = ""
  )
}

# Save the original point estimates as csv files
save_se_assay(se, "scaledTPMCounts", .rownames = "transcript_id", file = "transcript/scaled-counts.csv")
save_se_assay(gse, "scaledTPMCounts", .rownames = "gene_id", file = "gene/scaled-counts.csv")
save_se_assay(se_dtu, "dtuScaledTPMCounts", .rownames = "transcript_id", file = "dtu/scaled-counts.csv")

# Length files
save_se_assay(se, "length", .rownames = "transcript_id", file = "transcript/effective-lengths.csv.gz")
save_se_assay(gse, "length", .rownames = "gene_id", file = "gene/effective-lengths.csv.gz")
save_se_assay(se_dtu, "length", .rownames = "transcript_id", file = "dtu/effective-lengths.csv.gz")

# mean estimates
save_se_assay(se, "mean", .rownames = "transcript_id", file = "transcript/infreps-mean.csv.gz")
save_se_assay(gse, "mean", .rownames = "gene_id", file = "gene/infreps-mean.csv.gz")
save_se_assay(se_dtu, "mean", .rownames = "transcript_id", file = "dtu/infreps-mean.csv.gz")

# variance estimates
save_se_assay(se, "variance", .rownames = "transcript_id", file = "transcript/infreps-relative-variance.csv.gz")
save_se_assay(gse, "variance", .rownames = "gene_id", file = "gene/infreps-relative-variance.csv.gz")
save_se_assay(se_dtu, "variance", .rownames = "transcript_id", file = "dtu/infreps-relative-variance.csv.gz")

# Save the scaled infReps as an RDS file
save_se_infreps <- function(se, .rownames, .file) {
  infreps <-
    assays(se)[str_which(assayNames(se), "^infRep")] |>
    as.list() |>
    map(function(rep) as_tibble(rep, rownames = .rownames))

  readr::write_rds(infreps, file = .file, compress = "gz")
}

save_se_infreps(se, .rownames = "transcript_id", .file = "transcript/infreps.rds")
save_se_infreps(gse, .rownames = "gene_id", .file = "gene/infreps.rds")
save_se_infreps(se_dtu, .rownames = "transcript_id", .file = "dtu/infreps.rds")
