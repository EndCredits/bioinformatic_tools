#!/usr/bin/env Rscript
#
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2025-2026 EndCredits <admin@endcredits.cc>
#
# visualize_tfbs.R
# Transcription Factor Binding Site Visualization Tool
#
# Supports two input formats:
#   1. PlantRegMap TFBS format: Standard PlantRegMap TFBS prediction output (requires FASTA file)
#   2. PlantCARE format: PlantCARE cis-acting element prediction output
#
# Usage (PlantRegMap TFBS format):
#   Rscript visualize_tfbs.R <fasta_file> <prediction_file> <output_prefix> [pvalue_threshold]
#
# Usage (PlantCARE format):
#   Rscript visualize_tfbs.R --plantcare <plantcare_tab_file> <fasta_file|-> <output_prefix> [motif_filter]
#
# Options:
#   --plantcare            Use PlantCARE input format
#   --single               Merge plus/minus strands into single line view
#   --exclude MOTIFS       Comma-separated motifs to exclude (overrides config file)
#   --no-exclude           Disable all motif exclusions (include everything)
#   --config FILE          Path to YAML config file (default: tfbs_config.yaml in script dir)
#   --no-config            Ignore config file, use built-in defaults only
#
# Filter priority: exclude runs first, then motif whitelist. --exclude overrides config file.
#
# Examples:
#   # TFBS format
#   Rscript visualize_tfbs.R promoter.fasta predictions.txt output 1e-6
#
#   # PlantCARE format
#   Rscript visualize_tfbs.R --plantcare plantCARE_output.tab promoter.fa output
#   Rscript visualize_tfbs.R --plantcare plantCARE_output.tab - output   # no FASTA
#   Rscript visualize_tfbs.R --plantcare plantCARE_output.tab promoter.fa output TATA-box,ABRE,G-box
#
#   # Exclude specific motifs
#   Rscript visualize_tfbs.R --plantcare --exclude "G-box,ABRE" plantCARE_output.tab promoter.fa output
#
#   # Single strand view
#   Rscript visualize_tfbs.R --plantcare --single plantCARE_output.tab promoter.fa output
#
# Output:
#   <output_prefix>.pdf - main visualization
#   <output_prefix>.png - PNG version
#   <output_prefix>_summary.json - summary statistics
#   <output_prefix>_coordinates.csv - transformed coordinates
#
# Note: In --single mode, enriched_regions in _summary.json use merged strand labels
#       (e.g. "+/-" for overlapping plus/minus elements), while _coordinates.csv retains
#       the original per-element strand information.
#

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(jsonlite)
# yaml is loaded conditionally only when --config is used

MAX_LANES <- 16L

# ---- 0. Configuration ----

default_config <- function() {
  list(
    exclude_motifs = c(),
    remove_singletons = TRUE,
    significance = list(high = 1e-8, medium = 1e-7, low = 1e-6),
    enriched_min_count = 3,
    # Visualization parameters
    strand_y_offset = 0.4,
    stripe_width = 200,
    label_min_dist_ratio = 0.045,
    alpha_score_range = c(4, 12),
    lane_height = 0.15,
    merge_threshold = 0
  )
}

load_config <- function(config_file = NULL) {
  config <- default_config()

  if (!is.null(config_file) && file.exists(config_file)) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required to read config files. Install it with install.packages('yaml') or use --no-config.")
    }
    message("Loading config from: ", config_file)
    user_config <- yaml::read_yaml(config_file)
    if (!is.null(user_config$exclude_motifs)) {
      config$exclude_motifs <- unlist(user_config$exclude_motifs)
    }
    if (!is.null(user_config$remove_singletons)) {
      config$remove_singletons <- as.logical(user_config$remove_singletons)
    }
    if (!is.null(user_config$significance)) {
      config$significance$high <- as.numeric(user_config$significance$high)
      config$significance$medium <- as.numeric(user_config$significance$medium)
      config$significance$low <- as.numeric(user_config$significance$low)
    }
    if (!is.null(user_config$enriched_min_count)) {
      config$enriched_min_count <- as.integer(user_config$enriched_min_count)
    }
    if (!is.null(user_config$strand_y_offset)) {
      config$strand_y_offset <- as.numeric(user_config$strand_y_offset)
    }
    if (!is.null(user_config$stripe_width)) {
      config$stripe_width <- as.integer(user_config$stripe_width)
    }
    if (!is.null(user_config$label_min_dist_ratio)) {
      config$label_min_dist_ratio <- as.numeric(user_config$label_min_dist_ratio)
    }
    if (!is.null(user_config$alpha_score_range)) {
      config$alpha_score_range <- unlist(user_config$alpha_score_range)
    }
    if (!is.null(user_config$lane_height)) {
      config$lane_height <- as.numeric(user_config$lane_height)
    }
    if (!is.null(user_config$merge_threshold)) {
      config$merge_threshold <- as.integer(user_config$merge_threshold)
    }
    message("   exclude_motifs: ", paste(config$exclude_motifs, collapse = ", "))
    message("   remove_singletons: ", config$remove_singletons)
    message("   significance: ***/", config$significance$high,
            " /**/", config$significance$medium,
            " /*/", config$significance$low)
  } else if (!is.null(config_file)) {
    warning("Config file not found: ", config_file, ". Using defaults.")
  }

  config
}

# ---- Helper: parse first sequence from FASTA (handles multi-sequence files) ----

parse_first_fasta_sequence <- function(fasta_lines) {
  header_lines <- grep("^>", fasta_lines)
  if (length(header_lines) > 1) {
    warning("FASTA file contains ", length(header_lines), " sequences. Only the first sequence will be used.")
    first_seq_end <- header_lines[2] - 1
  } else {
    first_seq_end <- length(fasta_lines)
  }
  header <- fasta_lines[1]
  seq_lines <- trimws(fasta_lines[2:first_seq_end])
  seq_lines <- seq_lines[seq_lines != ""]
  sequence <- paste(seq_lines, collapse = "")
  if (nchar(sequence) == 0) stop("FASTA file contains no sequence data (header only or empty after header).")
  list(header = header, sequence = sequence, seq_length = nchar(sequence))
}

# ---- 1. Data loading & validation ----

load_data <- function(fasta_file, prediction_file) {
  fasta_lines <- readLines(fasta_file, warn = FALSE)
  fasta <- parse_first_fasta_sequence(fasta_lines)

  # Read prediction file (TSV, first line is #header)
  pred_lines <- readLines(prediction_file)
  header_line <- sub("^#", "", pred_lines[1])
  header_names <- unlist(strsplit(header_line, "\t"))

  pred_data <- read.table(prediction_file,
                          sep = "\t", header = FALSE, skip = 1,
                          stringsAsFactors = FALSE,
                          col.names = header_names,
                          check.names = FALSE)

  # Normalize column names
  colnames(pred_data) <- tolower(colnames(pred_data))
  colnames(pred_data) <- gsub("-", "_", colnames(pred_data))
  colnames(pred_data) <- gsub("\\.", "_", colnames(pred_data))

  required_cols <- c("family", "start", "stop", "strand", "score", "p_value", "q_value")
  missing_cols <- setdiff(required_cols, colnames(pred_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in prediction file: ", paste(missing_cols, collapse = ", "))
  }

  pred_data <- pred_data %>%
    rename(end = stop) %>%
    mutate(
      start = as.integer(start),
      end = as.integer(end),
      score = as.numeric(score),
      p_value = as.numeric(p_value),
      q_value = as.numeric(q_value),
      strand = as.character(strand),
      family = as.character(family)
    )

  if (any(pred_data$start < 1 | pred_data$end > fasta$seq_length, na.rm = TRUE)) {
    warning("Some coordinates are outside the sequence length range.")
  }

  list(header = fasta$header, sequence = fasta$sequence, seq_length = fasta$seq_length, pred_data = pred_data)
}

# ---- 1b. PlantCARE data loading ----

load_plantcare_data <- function(plantcare_file, fasta_file = NULL, motif_filter = NULL,
                                exclude_list = NULL, remove_singletons = TRUE) {
  message("Loading PlantCARE data from: ", plantcare_file)

  lines <- readLines(plantcare_file, warn = FALSE)

  # Detect and skip header lines
  # PlantCARE header lines have "organism" in col 7 and "short_function" in col 8
  is_header <- sapply(lines, function(line) {
    fields <- unlist(strsplit(line, "\t"))
    if (length(fields) >= 7 && tolower(fields[7]) == "organism") return(TRUE)
    if (length(fields) >= 8 && grepl("short_function", tolower(fields[8]))) return(TRUE)
    FALSE
  })
  if (any(is_header)) {
    message("   Detected and removing ", sum(is_header), " header line(s)")
    lines <- lines[!is_header]
  }

  parse_line <- function(line) {
    fields <- unlist(strsplit(line, "\t"))
    if (length(fields) < 7) return(NULL)

    start_pos <- as.integer(fields[4])
    length_val <- as.integer(fields[5])
    end_pos <- start_pos + length_val - 1

    data.frame(
      seq_name = fields[1],
      family = fields[2],
      motif_sequence = fields[3],
      start = start_pos,
      end = end_pos,
      length = length_val,
      strand = fields[6],
      organism = fields[7],
      description = if (length(fields) >= 8) paste(fields[8:length(fields)], collapse = "; ") else NA,
      score = as.numeric(length_val),
      p_value = NA_real_,
      q_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  parsed_list <- lapply(lines, parse_line)
  pred_data <- do.call(rbind, Filter(Negate(is.null), parsed_list))
  if (is.null(pred_data) || nrow(pred_data) == 0) {
    stop("No valid elements found in PlantCARE file. Check file format (expected tab-separated: seq_name, family, motif_sequence, start, length, strand, organism, ...).")
  }
  pred_data <- pred_data[!is.na(pred_data$family), ]
  message("   Loaded ", nrow(pred_data), " cis-acting elements")

  # Apply exclusion list
  if (!is.null(exclude_list) && length(exclude_list) > 0) {
    excluded <- pred_data$family %in% exclude_list
    if (sum(excluded) > 0) {
      message("   Excluded ", sum(excluded), " elements: ", paste(intersect(pred_data$family, exclude_list), collapse = ", "))
      pred_data <- pred_data[!excluded, ]
    }
  }

  # Read FASTA if provided, otherwise infer sequence length
  if (!is.null(fasta_file) && file.exists(fasta_file)) {
    fasta_lines <- readLines(fasta_file, warn = FALSE)
    fasta <- parse_first_fasta_sequence(fasta_lines)
    header <- fasta$header
    sequence <- fasta$sequence
    seq_length <- fasta$seq_length
    message("   Sequence length from FASTA: ", seq_length, " bp")
  } else {
    max_end <- max(pred_data$end, na.rm = TRUE)
    if (max_end > 100000) {
      warning("Inferred sequence length is unusually large (", max_end, " bp). Check coordinates in input file.")
    }
    seq_length <- max_end + 50
    message("   Inferred sequence length: ", seq_length, " bp (max element end + 50 bp margin)")
    header <- paste0(">", gsub("^PlantCARE_", "", pred_data$seq_name[1]))
    sequence <- NA
  }

  # Motif whitelist filter
  if (!is.null(motif_filter) && length(motif_filter) > 0 && motif_filter[1] != "") {
    filter_list <- unlist(strsplit(motif_filter, ","))
    matched <- pred_data %>% filter(family %in% filter_list)
    message("   Filtered to ", nrow(matched), " elements matching: ", paste(filter_list, collapse = ", "))
    if (nrow(matched) == 0) {
      warning("No elements matched the filter '", motif_filter, "'. Showing all remaining elements after exclusion/singleton filtering.")
    } else {
      pred_data <- matched
    }
  }

  # Remove singletons (families appearing only once)
  if (remove_singletons) {
    family_counts <- table(pred_data$family)
    singleton_families <- names(family_counts[family_counts == 1])
    if (length(singleton_families) > 0) {
      message("   Filtered out ", length(singleton_families), " singleton element(s)")
      pred_data <- pred_data %>% filter(!family %in% singleton_families)
    }
  }

  gene_name <- gsub("^>\\s*", "", header)
  gene_name <- gsub("\\|.*", "", gene_name)

  list(header = header, sequence = sequence, seq_length = seq_length,
       pred_data = pred_data, is_plantcare = TRUE)
}

# ---- 2. Coordinate transformation ----

process_coordinates <- function(pred_data, seq_length) {
  data <- pred_data %>%
    mutate(
      # For plus strand: plot_start=start (near TSS), plot_end=end (far from TSS)
      # For minus strand: mirror coordinates, plot_start is near TSS end
      plot_start = ifelse(strand == "+", start, seq_length - end + 1),
      plot_end = ifelse(strand == "+", end, seq_length - start + 1)
    )

  invalid_rows <- is.na(data$plot_start) | is.na(data$plot_end) |
    data$plot_start < 1 | data$plot_end > seq_length | data$plot_start > data$plot_end
  if (sum(invalid_rows, na.rm = TRUE) > 0) {
    warning(paste("Filtered out", sum(invalid_rows, na.rm = TRUE), "elements with invalid/NA coordinates"))
    data <- data[!invalid_rows, ]
  }

  data %>% mutate(width = plot_end - plot_start + 1)
}

# ---- 3. Site aggregation ----
# O(n log n) interval merge: sort then linear sweep

aggregate_sites <- function(coord_data, merge_threshold = 0, single_strand = FALSE) {
  aggregate_chain <- function(chain_data) {
    if (nrow(chain_data) == 0) return(data.frame())

    chain_data <- chain_data %>% arrange(plot_start, plot_end)
    n <- nrow(chain_data)

    # Linear sweep: maintain current cluster's max_end; a new interval joins
    # the current cluster iff plot_start <= cluster_max_end + merge_threshold.
    # Since input is sorted by plot_start, any transitive overlap is captured
    # because extending cluster_max_end can only pull in more intervals that
    # appear later in the sorted order.
    cluster_starts <- integer(0)
    cluster_ends <- integer(0)

    cluster_start <- 1L
    cluster_max_end <- chain_data$plot_end[1]

    for (i in if (n >= 2) 2:n else integer(0)) {
      if (chain_data$plot_start[i] <= cluster_max_end + merge_threshold) {
        # Overlaps or touches current cluster — extend
        cluster_max_end <- max(cluster_max_end, chain_data$plot_end[i])
      } else {
        # Gap — emit current cluster, start new one
        cluster_starts <- c(cluster_starts, cluster_start)
        cluster_ends <- c(cluster_ends, i - 1L)
        cluster_start <- i
        cluster_max_end <- chain_data$plot_end[i]
      }
    }
    # Emit final cluster
    cluster_starts <- c(cluster_starts, cluster_start)
    cluster_ends <- c(cluster_ends, n)

    # Build aggregated results from cluster ranges
    result_list <- vector("list", length(cluster_starts))
    for (k in seq_along(cluster_starts)) {
      group <- chain_data[cluster_starts[k]:cluster_ends[k], ]

      # Select primary family: best p-value if available, otherwise highest score
      primary_idx <- if (any(!is.na(group$p_value))) which.min(group$p_value) else which.max(group$score)
      mpv <- min(group$p_value, na.rm = TRUE)
      if (!is.finite(mpv)) mpv <- NA_real_

      # For single_strand mode, preserve strand info from source data
      if (single_strand) {
        unique_strands <- unique(group$strand)
        strand_label <- if (length(unique_strands) > 1) paste(sort(unique_strands), collapse = "/") else unique_strands[1]
      } else {
        strand_label <- chain_data$strand[cluster_starts[k]]
      }

      result_list[[k]] <- data.frame(
        plot_start = min(group$plot_start),
        plot_end = max(group$plot_end),
        strand = strand_label,
        family = paste(unique(group$family), collapse = ","),
        family_primary = group$family[primary_idx[1]],
        count = nrow(group),
        min_p_value = mpv,
        max_score = max(group$score),
        avg_score = mean(group$score),
        stringsAsFactors = FALSE
      )
    }
    bind_rows(result_list)
  }

  if (single_strand) {
    aggregate_chain(coord_data)
  } else {
    plus_data <- coord_data %>% filter(strand == "+")
    minus_data <- coord_data %>% filter(strand == "-")
    bind_rows(aggregate_chain(plus_data), aggregate_chain(minus_data)) %>%
      arrange(strand, plot_start)
  }
}

# ---- 4. Visualization data preparation ----

prepare_visual_data <- function(aggregated_data, seq_length, config, single_strand = FALSE) {
  sig <- config$significance
  is_plantcare <- all(is.na(aggregated_data$min_p_value))

  base <- aggregated_data %>%
    mutate(
      site_center = (plot_start + plot_end) / 2
    )

  # Strand line Y depends on mode
  if (single_strand) {
    base <- base %>% mutate(strand_line_y = 0)
  } else {
    base <- base %>% mutate(strand_line_y = ifelse(strand == "+", config$strand_y_offset, -config$strand_y_offset))
  }

  # Point offset to avoid overlap at identical positions
  base <- base %>%
    mutate(
      position_group = if (single_strand) {
        paste0("single_", round(plot_start / 10) * 10)
      } else {
        paste0(strand, "_", round(plot_start / 10) * 10)
      },
      group_index = ave(seq_along(position_group), position_group, FUN = seq_along),
      max_offset = if (single_strand) 0.03 else 0.08,
      point_offset = (group_index - 1) * max_offset / 2,
      offset_direction = ifelse(group_index %% 2 == 0, 1, -1),
      point_y = strand_line_y + point_offset * offset_direction
    )

  # Significance: use p-value for TFBS, use count for PlantCARE
  if (is_plantcare) {
    enriched_min <- config$enriched_min_count
    # Auto-adapt alpha range from actual motif length data (score == length in PlantCARE)
    score_range <- range(base$max_score, na.rm = TRUE)
    if (score_range[1] == score_range[2]) score_range <- score_range + c(-0.5, 0.5)
    base <- base %>%
      mutate(
        significance = case_when(
          count >= enriched_min * 3 ~ "***",
          count >= enriched_min * 2 ~ "**",
          count >= enriched_min ~ "*",
          TRUE ~ ""
        ),
        # PlantCARE: alpha based on motif length (score == length in PlantCARE)
        alpha_val = scales::rescale(max_score,
                                    to = c(0.3, 1.0),
                                    from = score_range)
      )
  } else {
    base <- base %>%
      mutate(
        significance = case_when(
          !is.na(min_p_value) & min_p_value < sig$high ~ "***",
          !is.na(min_p_value) & min_p_value < sig$medium ~ "**",
          !is.na(min_p_value) & min_p_value < sig$low ~ "*",
          TRUE ~ ""
        ),
        alpha_val = ifelse(!is.na(min_p_value), {
          from_range <- c(-log10(sig$low), max(-log10(sig$low) + 0.5, -log10(sig$high * 0.01)))
          scales::rescale(-log10(min_p_value), to = c(0.3, 1.0), from = from_range)
        }, scales::rescale(max_score, to = c(0.3, 1.0), from = config$alpha_score_range))
      )
  }
  base <- base %>% mutate(alpha_val = pmin(pmax(alpha_val, 0.3), 1.0))

  # Anti-overlap: assign Y-axis lanes
  base$lane <- 0

  if (single_strand) {
    # Sweep-line anti-overlap for single strand
    base <- base %>%
      mutate(
        label_len = nchar(family_primary),
        influence_radius = pmax(40, label_len * 15)
      )

    n <- nrow(base)
    base$lane_dir <- 0
    lane_safety <- list()

    for (i in seq_len(n)) {
      curr_center <- base$site_center[i]
      curr_radius <- base$influence_radius[i]

      best_lane <- NULL
      best_dir <- NULL
      for (lane in 0:(MAX_LANES - 1L)) {
        for (dir in c(1, -1)) {
          key <- paste0(lane, "_", dir)
          last_end <- if (is.null(lane_safety[[key]])) -Inf else lane_safety[[key]]
          if (curr_center - last_end > curr_radius) {
            best_lane <- lane
            best_dir <- dir
            break
          }
        }
        if (!is.null(best_lane)) break
      }

      if (is.null(best_lane)) {
        warning("Label lane overflow at position ", curr_center, ": all lanes occupied, labels may overlap")
        best_lane <- 0
        best_dir <- 1
      }
      base$lane[i] <- best_lane
      base$lane_dir[i] <- best_dir

      key <- paste0(best_lane, "_", best_dir)
      curr_end <- curr_center + curr_radius
      if (is.null(lane_safety[[key]]) || curr_end > lane_safety[[key]]) {
        lane_safety[[key]] <- curr_end
      }
    }
  } else {
    # Per-strand anti-overlap for dual strand mode (sweep-line, O(n))
    # Use label width as influence radius, like single-strand mode
    base <- base %>%
      mutate(
        label_len = nchar(family_primary),
        influence_radius = pmax(seq_length * config$label_min_dist_ratio, label_len * 8)
      )

    for (s in c("+", "-")) {
      idx <- which(base$strand == s)
      if (length(idx) < 2) next

      # Sort by site_center within this strand
      idx <- idx[order(base$site_center[idx])]
      lane_safety <- rep(-Inf, MAX_LANES)  # last right-boundary per lane

      for (k in seq_along(idx)) {
        curr_idx <- idx[k]
        curr_center <- base$site_center[curr_idx]
        curr_radius <- base$influence_radius[curr_idx]

        # Find first lane where current label's left edge is past previous label's right edge
        assigned_lane <- NA
        for (lane in 0:(MAX_LANES - 1L)) {
          if (curr_center - curr_radius > lane_safety[lane + 1]) {
            assigned_lane <- lane
            break
          }
        }
        if (is.na(assigned_lane)) {
          warning("Lane overflow: more than ", MAX_LANES, " overlapping labels on strand '", s, "', labels may overlap")
          assigned_lane <- 0
        }
        base$lane[curr_idx] <- assigned_lane
        lane_safety[assigned_lane + 1] <- curr_center + curr_radius
      }
    }
  }

  # Label positioning
  if (single_strand) {
    base <- base %>%
      mutate(
        base_label_offset = 0.08,
        dynamic_offset = base_label_offset + lane * config$lane_height,
        label_y = case_when(
          lane_dir == 1  ~ strand_line_y + dynamic_offset,
          lane_dir == -1 ~ strand_line_y - dynamic_offset,
          TRUE           ~ strand_line_y + dynamic_offset
        ),
        arrow_start_y = case_when(
          lane_dir == 1  ~ label_y - 0.015,
          lane_dir == -1 ~ label_y + 0.015,
          TRUE           ~ label_y - 0.015
        ),
        arrow_end_y = point_y,
        significance_y = ifelse(lane_dir == 1, label_y + 0.03, label_y - 0.03),
        count_label_y = ifelse(lane_dir == 1, point_y - 0.03, point_y + 0.03)
      )
  } else {
    base <- base %>%
      mutate(
        base_label_offset = 0.15,
        dynamic_offset = base_label_offset + lane * config$lane_height,
        label_y = ifelse(strand == "+",
                         strand_line_y + dynamic_offset + point_offset,
                         strand_line_y - dynamic_offset - point_offset),
        arrow_start_y = ifelse(strand == "+", label_y - 0.03, label_y + 0.03),
        arrow_end_y = point_y,
        significance_y = ifelse(strand == "+", label_y + 0.06, label_y - 0.06),
        count_label_y = ifelse(strand == "+", point_y - 0.06, point_y + 0.06)
      )
  }

  # Remove NA family_primary
  base <- base %>% filter(!is.na(family_primary))

  # Color mapping
  families <- unique(base$family_primary)
  n_families <- length(families)
  if (n_families == 0) {
    family_colors <- c("gray"); names(family_colors) <- c("Unknown")
    base$family_primary <- factor("Unknown")
  } else if (n_families <= 12) {
    family_colors <- setNames(viridis(n_families), families)
  } else {
    family_colors <- setNames(rep(viridis(12), ceiling(n_families / 12))[1:n_families], families)
  }

  list(plot_data = base, family_colors = family_colors, seq_length = seq_length, single_strand = single_strand)
}

# ---- 5. Plot generation ----

create_plot <- function(viz_data, header, total_count, mode = "tfbs", config) {
  plot_data <- viz_data$plot_data
  family_colors <- viz_data$family_colors
  seq_length <- viz_data$seq_length
  single_strand <- viz_data$single_strand
  enriched_min <- config$enriched_min_count

  gene_name <- gsub("^>\\s*", "", header)
  gene_name <- gsub("\\|.*", "", gene_name)

  plot_data$family_primary <- factor(plot_data$family_primary, levels = names(family_colors))

  # Title and labels depend on mode
  n_motifs <- length(unique(plot_data$family_primary))
  if (single_strand) {
    if (mode == "plantcare") {
      title_text <- paste("PlantCARE Cis-Acting Element Prediction -", gene_name)
      caption_text <- "Cis-acting elements on single line (strand merged) | Opacity based on motif length | n=x indicates multiple elements at same position"
    } else {
      sig <- config$significance
      title_text <- paste("Transcription Factor Binding Site Predictions -", gene_name)
      caption_text <- paste0("TFBS on single line (strand merged) | ",
                             "*** p<", sig$high, ", ** p<", sig$medium, ", * p<", sig$low,
                             " | n=x indicates multiple predictions at same site")
    }
    subtitle_text <- paste0("Length: ", seq_length, " bp | Total elements: ", total_count, " | Element types: ", n_motifs)
  } else if (mode == "plantcare") {
    title_text <- paste("PlantCARE Cis-Acting Element Prediction -", gene_name)
    subtitle_text <- paste0("Length: ", seq_length, " bp | Total elements: ", total_count, " | Element types: ", n_motifs)
    caption_text <- "Arrows indicate cis-acting elements | Opacity based on motif length | n=x indicates multiple elements at same position"
  } else {
    title_text <- paste("Transcription Factor Binding Site Predictions -", gene_name)
    sig <- config$significance
    subtitle_text <- paste0("Length: ", seq_length, " bp | Filter: p < ", format(sig$low, scientific = TRUE), " | Total predictions: ", total_count)
    caption_text <- paste0("Arrows indicate transcription factor binding sites | ",
                           "*** p<", sig$high, ", ** p<", sig$medium, ", * p<", sig$low,
                           " | n=x indicates multiple predictions at same site")
  }

  p <- ggplot()

  # Background stripes
  stripe_breaks <- seq(0, seq_length, by = config$stripe_width)
  stripe_data <- data.frame(
    xmin = stripe_breaks,
    xmax = pmin(stripe_breaks + config$stripe_width, seq_length)
  ) %>% filter(xmin < seq_length)

  if (nrow(stripe_data) > 0) {
    p <- p + geom_rect(
      data = stripe_data,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = rep(c("gray95", "white"), length.out = nrow(stripe_data)),
      alpha = 0.5
    )
  }

  # TSS marker
  p <- p + geom_vline(xintercept = 1, color = "red", linetype = "dashed", linewidth = 0.5)

  # DNA strand lines
  if (single_strand) {
    p <- p + geom_segment(
      aes(x = 0, xend = seq_length, y = 0, yend = 0),
      color = "black", linewidth = 0.8
    )
  } else {
    dna_strands <- data.frame(strand = c("+", "-"), y = c(config$strand_y_offset, -config$strand_y_offset))
    p <- p + geom_segment(
      data = dna_strands,
      aes(x = 0, xend = seq_length, y = y, yend = y),
      color = "black", linewidth = 0.8
    )
  }

  # Site markers
  p <- p + geom_segment(
    data = plot_data,
    aes(x = site_center, xend = site_center, y = strand_line_y, yend = point_y, color = family_primary),
    linewidth = 0.6, alpha = 0.7
  )

  # Site dots
  p <- p + geom_point(
    data = plot_data,
    aes(x = site_center, y = point_y, color = family_primary, alpha = alpha_val),
    size = 2.5
  )

  # Family labels
  p <- p + geom_text(
    data = plot_data,
    aes(x = site_center, y = label_y, label = family_primary, color = family_primary),
    size = 3, fontface = "bold", vjust = 0.5, hjust = 0.5
  )

  # Arrows from labels to sites
  p <- p + geom_segment(
    data = plot_data,
    aes(x = site_center, xend = site_center, y = arrow_start_y, yend = arrow_end_y, color = family_primary),
    linewidth = 0.3,
    arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
    alpha = 0.6
  )

  # Significance markers
  p <- p + geom_text(
    data = plot_data %>% filter(significance != ""),
    aes(x = site_center, y = significance_y, label = significance),
    size = 2.5, color = "red", fontface = "bold"
  )

  # Enriched region count labels
  if (enriched_min > 0) {
    p <- p + geom_text(
      data = plot_data %>% filter(count >= enriched_min),
      aes(x = site_center, y = count_label_y, label = paste0("n=", count)),
      size = 2, color = "darkgray", fontface = "italic"
    )
  }

  # Dynamic Y-axis range based on actual label positions
  max_label_y <- max(plot_data$label_y, na.rm = TRUE)
  min_label_y <- min(plot_data$label_y, na.rm = TRUE)
  y_floor <- if (single_strand) 0.8 else 1.2
  y_limit <- max(y_floor, max(abs(max_label_y), abs(min_label_y)) * 1.15)

  # Y-axis configuration
  if (single_strand) {
    p <- p +
      scale_y_continuous(name = "", breaks = NULL, limits = c(-y_limit, y_limit))
  } else {
    p <- p +
      scale_y_continuous(
        name = "",
        breaks = c(config$strand_y_offset, -config$strand_y_offset),
        labels = c("Plus strand", "Minus strand"),
        limits = c(-y_limit, y_limit)
      )
  }

  p <- p +
    scale_color_manual(values = family_colors) +
    scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    scale_x_continuous(
      name = "Position (bp)",
      breaks = seq(0, seq_length, by = ifelse(seq_length > 2000, 500, 100)),
      limits = c(0, seq_length),
      expand = expansion(mult = 0.02)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50")
    ) +
    labs(title = title_text, subtitle = subtitle_text, caption = caption_text)

  # Single strand extras
  if (single_strand) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        plot.caption = element_text(hjust = 0.5, color = "gray50"),
        plot.margin = margin(t = 10, r = 10, b = 25, l = 10, unit = "pt")
      )
  }

  p
}

# ---- 6. Result export ----

export_results <- function(processed_data, aggregated_data, viz_data,
                           output_prefix, seq_length, total_pred_count,
                           mode = "tfbs", raw_data = NULL, config, plot_obj = NULL) {
  enriched_min <- config$enriched_min_count

  # Ensure output directory exists
  out_dir <- dirname(output_prefix)
  if (out_dir != "." && !dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir)
  }

  # Coordinates
  coords_file <- paste0(output_prefix, "_coordinates.csv")
  write.csv(processed_data, coords_file, row.names = FALSE)
  message("Saved transformed coordinates to: ", coords_file)

  # Raw mapping (PlantCARE only)
  mapping_file <- NULL
  if (mode == "plantcare" && !is.null(raw_data)) {
    mapping_file <- paste0(output_prefix, "_raw_mapping.csv")
    write.csv(raw_data, mapping_file, row.names = FALSE)
    message("Saved raw data mapping to: ", mapping_file)
  }

  # Helper: convert data.frame to list-of-lists for JSON serialization
  df_to_records <- function(df) {
    if (nrow(df) == 0) return(list())
    # Convert factors to character to avoid serializing as integer indices
    df <- as.data.frame(lapply(df, function(col) if (is.factor(col)) as.character(col) else col),
                         stringsAsFactors = FALSE)
    unname(split(df, seq_len(nrow(df)))) %>% lapply(function(row) as.list(row))
  }

  # Summary
  if (mode == "plantcare") {
    motif_counts <- raw_data %>%
      group_by(family) %>%
      summarise(count = n(), avg_length = mean(length),
                min_pos = min(start), max_pos = max(end),
                organisms = paste(unique(organism), collapse = "; "),
                .groups = "drop") %>%
      arrange(desc(count))

    strand_summary <- raw_data %>% group_by(strand) %>% summarise(count = n(), .groups = "drop")
    plus_count <- sum(strand_summary$count[strand_summary$strand == "+"])
    minus_count <- sum(strand_summary$count[strand_summary$strand == "-"])

    enriched_df <- aggregated_data %>%
      filter(count >= enriched_min) %>%
      select(plot_start, plot_end, strand, family_primary, count) %>%
      as.data.frame()

    motif_summary_list <- list()
    if (nrow(motif_counts) > 0) {
      motif_summary_list <- motif_counts %>%
        mutate(avg_length = round(avg_length, 2)) %>%
        as.data.frame() %>%
        split(.$family) %>%
        lapply(function(x) list(count = x$count, avg_length = x$avg_length, organisms = x$organisms))
    }

    summary_stats <- list(
      sequence_length = seq_length,
      total_elements = total_pred_count,
      plus_strand_count = plus_count,
      minus_strand_count = minus_count,
      unique_motif_types = length(unique(raw_data$family)),
      motif_distribution = { ft <- table(raw_data$family); setNames(as.list(as.integer(ft)), names(ft)) },
      motif_summary = motif_summary_list,
      enriched_regions = df_to_records(enriched_df),
      timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  } else {
    enriched_df <- aggregated_data %>%
      filter(count >= enriched_min) %>%
      select(plot_start, plot_end, strand, family_primary, count, min_p_value) %>%
      mutate(p_value = format(min_p_value, scientific = TRUE)) %>%
      select(-min_p_value) %>%
      as.data.frame()

    summary_stats <- list(
      sequence_length = seq_length,
      total_predictions = total_pred_count,
      filtered_predictions = nrow(processed_data),
      plus_strand_count = sum(processed_data$strand == "+"),
      minus_strand_count = sum(processed_data$strand == "-"),
      family_distribution = { ft <- table(processed_data$family); setNames(as.list(as.integer(ft)), names(ft)) },
      enriched_regions = df_to_records(enriched_df),
      timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  }

  summary_file <- paste0(output_prefix, "_summary.json")
  write_json(summary_stats, summary_file, pretty = TRUE, auto_unbox = TRUE)
  message("Saved summary statistics to: ", summary_file)

  # Plot dimensions
  if (mode == "plantcare") {
    plot_width <- max(10, ceiling(seq_length / 300) * 1.5)
  } else {
    plot_width <- 10
  }
  plot_height <- 6

  # Use explicit plot_obj parameter; fall back to viz_data$plot_obj for compatibility
  p <- if (!is.null(plot_obj)) plot_obj else viz_data$plot_obj
  if (is.null(p)) stop("No plot object provided to export_results. Pass plot_obj or set viz_data$plot_obj.")

  pdf_file <- paste0(output_prefix, ".pdf")
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height, device = "pdf")
  message("Saved PDF visualization to: ", pdf_file)

  png_file <- paste0(output_prefix, ".png")
  ggsave(png_file, plot = p, width = plot_width, height = plot_height, dpi = 300, device = "png")
  message("Saved PNG visualization to: ", png_file)

  out <- list(coordinates = coords_file, summary = summary_file, pdf = pdf_file, png = png_file)
  if (!is.null(mapping_file)) out$raw_mapping <- mapping_file
  out
}

# ---- 7. CLI argument parsing ----

parse_args <- function(args) {
  opts <- list(
    plantcare = FALSE,
    single = FALSE,
    exclude = NULL,
    config = NULL,
    no_config = FALSE,
    positional = character()
  )

  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--plantcare") {
      opts$plantcare <- TRUE
    } else if (a == "--single") {
      opts$single <- TRUE
    } else if (a == "--no-exclude") {
      opts$exclude <- character(0)  # exclude nothing
    } else if (a == "--exclude") {
      if (i < length(args)) {
        raw <- args[i + 1]
        if (nchar(trimws(raw)) == 0) {
          opts$exclude <- character(0)  # empty string means exclude nothing
        } else {
          opts$exclude <- unlist(strsplit(raw, ","))
        }
        i <- i + 1
      } else {
        stop("--exclude requires a comma-separated argument")
      }
    } else if (a == "--config") {
      if (i < length(args)) {
        opts$config <- args[i + 1]
        i <- i + 1
      } else {
        stop("--config requires a file path argument")
      }
    } else if (a == "--no-config") {
      opts$no_config <- TRUE
    } else if (startsWith(a, "--")) {
      stop("Unknown option: ", a, ". Use --help for usage information.")
    } else {
      opts$positional <- c(opts$positional, a)
    }
    i <- i + 1
  }

  opts
}

# ---- Main ----

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("Usage: Rscript visualize_tfbs.R [options] <input_files...>\n")
    cat("\nModes:\n")
    cat("  TFBS mode:      <fasta_file> <prediction_file> <output_prefix>\n")
    cat("  PlantCARE mode: --plantcare <plantcare_tab> <fasta_file|-> <output_prefix>\n")
    cat("\nOptions:\n")
    cat("  --plantcare          Use PlantCARE input format\n")
    cat("  --single             Merge strands into single line view\n")
    cat("  --exclude MOTIFS     Comma-separated motifs to exclude (overrides config)\n")
    cat("  --no-exclude         Disable all motif exclusions (include everything)\n")
    cat("  --config FILE        Path to YAML config file (default: tfbs_config.yaml in script dir)\n")
    cat("  --no-config          Ignore config file, use built-in defaults only\n")
    cat("\nFilter priority: exclude runs first, then motif whitelist. --exclude overrides config file.\n")
    cat("\nExamples:\n")
    cat("  # TFBS format\n")
    cat("  Rscript visualize_tfbs.R promoter.fasta predictions.txt output 1e-6\n")
    cat("\n  # PlantCARE format\n")
    cat("  Rscript visualize_tfbs.R --plantcare plantCARE_output.tab promoter.fa output\n")
    cat("  Rscript visualize_tfbs.R --plantcare plantCARE_output.tab - output\n")
    cat("  Rscript visualize_tfbs.R --plantcare plantCARE_output.tab promoter.fa output TATA-box,ABRE,G-box\n")
    cat("\n  # Exclude motifs via CLI\n")
    cat("  Rscript visualize_tfbs.R --plantcare --exclude \"G-box,ABRE\" plantCARE_output.tab promoter.fa output\n")
    cat("\n  # Single strand view\n")
    cat("  Rscript visualize_tfbs.R --plantcare --single plantCARE_output.tab promoter.fa output\n")
    quit(status = 1)
  }

  opts <- parse_args(args)

  # Load config
  if (opts$no_config) {
    config <- default_config()
  } else {
    config_file <- if (!is.null(opts$config)) opts$config else {
      # Determine script directory for default config location
      script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
      if (is.null(script_dir) || is.na(script_dir)) script_dir <- "."
      default_path <- file.path(script_dir, "tfbs_config.yaml")
      if (file.exists(default_path)) default_path else NULL
    }
    config <- load_config(config_file)
  }

  # Apply --exclude override
  if (!is.null(opts$exclude)) {
    config$exclude_motifs <- opts$exclude
  }

  pos <- opts$positional

  if (opts$plantcare) {
    # PlantCARE mode
    if (length(pos) < 3) {
      stop("PlantCARE mode requires: --plantcare <plantcare_tab> <fasta_file|-> <output_prefix> [motif_filter]")
    }

    plantcare_file <- pos[1]
    fasta_file <- if (pos[2] == "-") NULL else pos[2]
    output_prefix <- pos[3]
    motif_filter <- if (length(pos) >= 4) pos[4] else NULL

    # Validate input files
    if (!file.exists(plantcare_file)) stop("PlantCARE file not found: ", plantcare_file)
    if (!is.null(fasta_file) && !file.exists(fasta_file)) stop("FASTA file not found: ", fasta_file)

    message("Starting PlantCARE visualization...")
    message("  PlantCARE file: ", plantcare_file)
    if (!is.null(fasta_file)) message("  FASTA file: ", fasta_file)
    else message("  FASTA file: (none, will infer from data)")
    message("  Output prefix: ", output_prefix)
    if (!is.null(motif_filter)) message("  Motif filter: ", motif_filter)
    if (opts$single) message("  Mode: Single strand (strand merged)")

    # 1. Load data
    message("\n1. Loading PlantCARE data...")
    data <- load_plantcare_data(plantcare_file, fasta_file, motif_filter,
                                exclude_list = config$exclude_motifs,
                                remove_singletons = config$remove_singletons)
    filtered_count <- nrow(data$pred_data)
    message(paste0("   Using all ", filtered_count, " cis-acting elements"))

    # 2. Coordinate transformation
    message("2. Processing coordinates...")
    processed_data <- process_coordinates(data$pred_data, data$seq_length)

    # 3. Aggregation
    message("3. Aggregating overlapping sites...")
    aggregated_data <- aggregate_sites(processed_data,
                                       merge_threshold = config$merge_threshold,
                                       single_strand = opts$single)

    # 4. Visualization data
    message("4. Preparing visualization data...")
    viz_data <- prepare_visual_data(aggregated_data, data$seq_length, config,
                                    single_strand = opts$single)

    # 5. Create plot
    message("5. Creating visualization plot...")
    p <- create_plot(viz_data, data$header, filtered_count, mode = "plantcare", config)
    viz_data$plot_obj <- p

    # 6. Export
    message("6. Exporting results...")
    output_files <- export_results(processed_data, aggregated_data, viz_data,
                                   output_prefix, data$seq_length, filtered_count,
                                   mode = "plantcare", raw_data = data$pred_data, config,
                                   plot_obj = p)

  } else {
    # TFBS mode
    if (length(pos) < 3) {
      stop("TFBS mode requires: <fasta_file> <prediction_file> <output_prefix> [pvalue_threshold]")
    }

    fasta_file <- pos[1]
    prediction_file <- pos[2]
    output_prefix <- pos[3]
    pvalue_threshold <- ifelse(length(pos) >= 4, as.numeric(pos[4]), config$significance$low)

    # If user specified p-value threshold on CLI, align significance thresholds
    if (length(pos) >= 4) {
      config$significance$low <- pvalue_threshold
      config$significance$medium <- pvalue_threshold / 10
      config$significance$high <- pvalue_threshold / 100
    }

    # Validate input files
    if (!file.exists(fasta_file)) stop("FASTA file not found: ", fasta_file)
    if (!file.exists(prediction_file)) stop("Prediction file not found: ", prediction_file)

    message("Starting TFBS visualization...")
    message("  FASTA file: ", fasta_file)
    message("  Prediction file: ", prediction_file)
    message("  Output prefix: ", output_prefix)
    message("  P-value threshold: ", pvalue_threshold)

    # 1. Load data
    message("\n1. Loading data...")
    data <- load_data(fasta_file, prediction_file)

    original_count <- nrow(data$pred_data)
    filtered_data <- data$pred_data %>% filter(p_value < pvalue_threshold)
    filtered_count <- nrow(filtered_data)
    message(paste0("   Loaded ", original_count, " predictions, ",
                   filtered_count, " passed p-value filter (p < ", pvalue_threshold, ")"))

    if (filtered_count == 0) {
      stop("No predictions passed the p-value threshold (p < ", pvalue_threshold, "). Nothing to visualize.")
    }

    # 2. Coordinate transformation
    message("2. Processing coordinates...")
    processed_data <- process_coordinates(filtered_data, data$seq_length)

    # 3. Aggregation
    message("3. Aggregating overlapping sites...")
    aggregated_data <- aggregate_sites(processed_data, merge_threshold = config$merge_threshold,
                                       single_strand = opts$single)

    # 4. Visualization data
    message("4. Preparing visualization data...")
    viz_data <- prepare_visual_data(aggregated_data, data$seq_length, config,
                                    single_strand = opts$single)

    # 5. Create plot
    message("5. Creating visualization plot...")
    p <- create_plot(viz_data, data$header, filtered_count, mode = "tfbs", config)
    viz_data$plot_obj <- p

    # 6. Export
    message("6. Exporting results...")
    output_files <- export_results(processed_data, aggregated_data, viz_data,
                                   output_prefix, data$seq_length, filtered_count,
                                   mode = "tfbs", config = config, plot_obj = p)
  }

  message("\nDone! Output files:")
  for (file_type in names(output_files)) {
    message(paste0("  ", file_type, ": ", output_files[[file_type]]))
  }

  if (interactive()) {
    print(p)
  }
}

if (sys.nframe() == 0) {
  main()
}
