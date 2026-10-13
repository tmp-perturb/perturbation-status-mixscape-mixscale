#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: run_shard.R <config.yaml> <frozen_shard.rds> <output_dir>")

suppressPackageStartupMessages({
  library(Seurat)
  library(Mixscale)
  library(yaml)
  library(jsonlite)
})

config <- read_yaml(args[[1]])
checkpoint <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stage_log <- file.path(output_dir, "stage_resources.jsonl")
params <- config$mixscale
set.seed(as.integer(config$seed))

record_stage <- function(stage, event, object = NULL) {
  status <- readLines("/proc/self/status", warn = FALSE)
  extract <- function(key) {
    line <- status[startsWith(status, paste0(key, ":"))]
    if (!length(line)) return(NA_real_)
    as.numeric(sub(".*:[[:space:]]*([0-9]+).*", "\\1", line[[1]])) * 1024
  }
  record <- list(timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), epoch = as.numeric(Sys.time()),
                 stage = stage, event = event,
                 process = list(VmRSS = extract("VmRSS"), VmHWM = extract("VmHWM"), VmSize = extract("VmSize")),
                 object_bytes = if (is.null(object)) NULL else as.numeric(object.size(object)))
  write(toJSON(record, auto_unbox = TRUE, null = "null"), file = stage_log, append = TRUE)
}

atomic_write <- function(writer, destination) {
  temporary <- paste0(destination, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  writer(temporary)
  if (!file.rename(temporary, destination)) stop("Could not atomically rename: ", destination)
}

record_stage("load_frozen_shard", "start")
object <- readRDS(checkpoint)
required_metadata <- c(params$labels, "is_ntc", "assignment_structure")
if (!inherits(object, "Seurat") || !all(c("RNA", params$input_assay) %in% names(object@assays))) stop("Invalid frozen shard")
if (!all(required_metadata %in% colnames(object[[]]))) stop("Frozen shard lacks required metadata")
if (!all(object$assignment_structure %in% config$selection$eligible_assignment_structures)) stop("Frozen shard contains ineligible cells")
if (!any(object$is_ntc)) stop("Frozen shard contains no eligible NTC cells")
if (!all(as.character(object[[params$labels]][object$is_ntc, 1]) == params$nt_class_name)) stop("Invalid NTC labels")
object$con1 <- "con1"
record_stage("load_frozen_shard", "end", object)

record_stage("RunMixscale", "start", object)
object <- RunMixscale(
  object = object, assay = params$input_assay, slot = params$score_slot,
  labels = params$labels, nt.class.name = params$nt_class_name,
  new.class.name = params$new_class_name, min.de.genes = as.integer(params$min_de_genes),
  min.cells = as.integer(params$min_cells), de.assay = params$de_assay,
  logfc.threshold = as.numeric(params$logfc_threshold), verbose = FALSE,
  split.by = "con1", fine.mode = isTRUE(params$fine_mode),
  max.de.genes = as.integer(params$max_de_genes), harmonize = isTRUE(params$harmonize),
  min_prop_ntgd = as.numeric(params$min_prop_ntgd), pval.cutoff = as.numeric(params$pval_cutoff),
  seed = as.integer(config$seed)
)
record_stage("RunMixscale", "end", object)

metadata <- object[[]]
metadata$con1 <- NULL
tool_result <- Tool(object, slot = "RunMixscale")
atomic_write(function(path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(metadata, connection, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}, file.path(output_dir, "mixscale_metadata.tsv.gz"))
atomic_write(function(path) saveRDS(tool_result, path, compress = FALSE), file.path(output_dir, "mixscale_tool.rds"))
atomic_write(function(path) writeLines(capture.output(sessionInfo()), path), file.path(output_dir, "sessionInfo.txt"))
atomic_write(function(path) write(toJSON(list(status = "complete", checkpoint = checkpoint,
                                               cells = ncol(object), targets = length(tool_result)),
                                         auto_unbox = TRUE, pretty = TRUE), path),
             file.path(output_dir, "completion.json"))
record_stage("write_output", "end", object)
