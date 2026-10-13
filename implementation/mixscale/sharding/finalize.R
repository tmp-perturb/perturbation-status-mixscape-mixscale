#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: finalize.R <config.yaml> <output_root>")

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

config <- read_yaml(args[[1]])
output_root <- normalizePath(args[[2]], mustWork = TRUE)
params <- config$mixscale
shard_dirs <- sort(list.dirs(file.path(output_root, "shards"), recursive = FALSE, full.names = TRUE))
if (!length(shard_dirs)) stop("No completed shard directories")

atomic_write <- function(writer, destination) {
  temporary <- paste0(destination, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  writer(temporary)
  if (!file.rename(temporary, destination)) stop("Could not atomically rename: ", destination)
}

metadata_parts <- list()
tool_result <- list()
ntc_reference <- NULL
target_cells <- character()
for (directory in shard_dirs) {
  required <- file.path(directory, c("completion.json", "mixscale_metadata.tsv.gz", "mixscale_tool.rds"))
  if (!all(file.exists(required))) stop("Incomplete shard output: ", directory)
  metadata <- read.delim(gzfile(required[[2]]), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  score <- metadata[[params$new_class_name]]
  names(score) <- rownames(metadata)
  ntc <- as.logical(metadata$is_ntc)
  if (is.null(ntc_reference)) {
    ntc_reference <- score[ntc]
  } else {
    common <- intersect(names(ntc_reference), names(score)[ntc])
    if (!identical(unname(ntc_reference[common]), unname(score[common]))) stop("NTC scores differ across shards")
  }
  target_score <- score[!ntc]
  if (length(intersect(target_cells, names(target_score)))) stop("Duplicate target cells across shards")
  target_cells <- c(target_cells, names(target_score))
  metadata_parts[[length(metadata_parts) + 1L]] <- target_score
  shard_tool <- readRDS(required[[3]])
  if (length(intersect(names(tool_result), names(shard_tool)))) stop("Duplicate target entries across Tool results")
  tool_result <- c(tool_result, shard_tool)
}

scores <- c(ntc_reference, unlist(metadata_parts, use.names = TRUE))
canonical <- read.delim(gzfile(config$canonical_metadata), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
missing <- setdiff(rownames(canonical), names(scores))
unexpected <- setdiff(names(scores), rownames(canonical))
if (length(missing) || length(unexpected) || anyDuplicated(names(scores))) stop("Merged score coverage validation failed")
canonical[[params$new_class_name]] <- as.numeric(scores[rownames(canonical)])

merged_dir <- file.path(output_root, "merged")
dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
atomic_write(function(path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(canonical, connection, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}, file.path(merged_dir, "mixscale_metadata.tsv.gz"))
atomic_write(function(path) saveRDS(tool_result, path, compress = FALSE), file.path(merged_dir, "mixscale_tool.rds"))
validation <- list(status = "complete", shards = length(shard_dirs), cells = nrow(canonical),
                   ntc = sum(as.logical(canonical$is_ntc)), targets = length(tool_result),
                   missing_cells = length(missing), unexpected_cells = length(unexpected),
                   duplicate_cells = anyDuplicated(names(scores)), ntc_agreement = "verified")
atomic_write(function(path) write(toJSON(validation, auto_unbox = TRUE, pretty = TRUE), path),
             file.path(merged_dir, "merge_validation.json"))
cat(toJSON(validation, auto_unbox = TRUE, pretty = TRUE), "\n")
