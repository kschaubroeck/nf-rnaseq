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

# The script accepts three possible variations:
#    1) samples and the directory to find the samples. This option requires a
#       sample flag to be passed to the script.
#    2) Just a single directory. If this option is selected, all folders in the
#       directory will be samples and the sample id will be the folder name.
#    3) A list of directories. Similar to option 2, the directory names will be
#       used as the sample names
#    4) A list of samples and directories
if (length(args) == 0) {
  stop("No arguments were passed to the script.", call. = FALSE)
}

## List of sample and Directories ---------------------------
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
