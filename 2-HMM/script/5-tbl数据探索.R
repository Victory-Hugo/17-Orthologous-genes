#!/usr/bin/env Rscript

# HMM matching visualization helper (ggplot2 edition)

suppressWarnings(options(stringsAsFactors = FALSE))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list(input = NULL, output = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-i", "--input")) {
      if (i + 1 > length(args)) stop("Value missing for --input")
      parsed$input <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (key %in% c("-o", "--output")) {
      if (i + 1 > length(args)) stop("Value missing for --output")
      parsed$output <- args[[i + 1]]
      i <- i + 2
      next
    }
    stop(sprintf("Unknown argument: %s", key))
  }
  if (is.null(parsed$input) || is.null(parsed$output)) {
    stop("Usage: Rscript 5-tbl数据探索.R --input <tbl.tsv> --output <dir>")
  }
  parsed
}

read_tbl <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file not found: %s", path))
  }
  cat(sprintf("[INFO] Reading data: %s\n", path))
  df <- tryCatch(
    read.delim(path, check.names = FALSE),
    error = function(e) stop(sprintf("Failed to read TSV: %s", e$message))
  )
  if (anyDuplicated(names(df))) {
    dup_names <- unique(names(df)[duplicated(names(df))])
    cat(sprintf("[WARN] Duplicate column names detected: %s\n",
                paste(dup_names, collapse = ", ")))
    names(df) <- make.unique(names(df))
  }
  cat(sprintf("[OK] Loaded %d rows / %d columns\n", nrow(df), ncol(df)))
  df
}

safe_log10 <- function(values) {
  finite_mask <- is.finite(values) & values > 0
  if (!any(finite_mask)) {
    return(rep(-Inf, length(values)))
  }
  min_positive <- min(values[finite_mask])
  safe_vals <- values
  safe_vals[!finite_mask] <- min_positive * 0.1
  log10(safe_vals)
}

prepare_evalue_data <- function(df) {
  breaks <- c(1e-120, 1e-110, 1e-100, 1e-90, 1e-80, 1e-70, 1e-60, 1e-50,
              1e-40, 1e-30, 1e-20, 1e-10, 1e-5, 1e-3, Inf)
  labels <- c("<1e-110", "1e-110~1e-100", "1e-100~1e-90", "1e-90~1e-80",
              "1e-80~1e-70", "1e-70~1e-60", "1e-60~1e-50", "1e-50~1e-40",
              "1e-40~1e-30", "1e-30~1e-20", "1e-20~1e-10", "1e-10~1e-5",
              "1e-5~1e-3", ">=1e-3")
  df %>%
    mutate(`E-value bin` = cut(`E-value`, breaks = breaks, labels = labels,
                               include.lowest = TRUE, right = FALSE)) %>%
    count(`E-value bin`, name = "count", .drop = FALSE) %>%
    mutate(`E-value bin` = factor(`E-value bin`, levels = labels),
           count = ifelse(count > 0, count, NA_real_))
}

plot_evalue <- function(df) {
  palette <- c("#085FE3", "#085FE3", "#085FE3", "#00AFFF", "#00AFFF",
               "#23A5AC", "#23A5AC", "#DB9C15", "#DB9C15", "#7A1616",
               "#7A1616", "#7A1616", "#7A1616", "#7A1616")
  ggplot(df, aes(x = `E-value bin`, y = count, fill = `E-value bin`)) +
    geom_col(width = 0.8, alpha = 0.5) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), trans = "log10") +
    scale_fill_manual(values = palette, guide = "none") +
    labs(title = "Distribution of E-value",
         x = NULL,
         y = "Matching count (log scale)") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1))
}

plot_score <- function(df) {
  ggplot(df, aes(x = score, fill = `query name`, color = `query name`)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 60) +
    scale_fill_brewer(palette = "Set2", name = "Query gene") +
    scale_color_brewer(palette = "Set2", guide = "none") +
    labs(title = "Distribution of score",
         x = "Score",
         y = "Matching count") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = "top")
}

plot_bias <- function(df) {
  ggplot(df, aes(x = bias)) +
    geom_histogram(fill = "#23A5AC", color = "white", bins = 80, alpha = 0.5) +
    labs(title = "Distribution of bias",
         x = "Bias",
         y = "Matching count") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())
}

plot_scatter <- function(df) {
  if (nrow(df) == 0) {
    return(ggplot() + geom_blank() + labs(title = "Score vs bias (no data)") + theme_void())
  }
  set.seed(42)
  sample_size <- min(5000, nrow(df))
  df_sample <- df %>%
    slice_sample(n = sample_size)
  df_sample <- df_sample %>%
    mutate(log10_evalue = safe_log10(`E-value`))
  ggplot(df_sample, aes(x = bias, y = score, color = log10_evalue)) +
    geom_point(size = 1.2, alpha = 0.5) +
    scale_color_gradientn(colors = c("#085FE3", "#00AFFF", "#23A5AC", "#DB9C15", "#7A1616"),
                          name = "log10(E-value)") +
    labs(title = "Score vs bias",
         x = "Bias",
         y = "Score") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())
}

build_dashboard <- function(df) {
  p1 <- df %>% prepare_evalue_data() %>% plot_evalue()
  p2 <- plot_score(df)
  p3 <- plot_bias(df)
  p4 <- plot_scatter(df)
  (p1 | p2) / (p3 | p4) + plot_annotation(title = "HMM matching quality assessment")
}

create_figures <- function(df, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  composed <- build_dashboard(df)
  pdf_path <- file.path(output_dir, "tbl_data_exploration.pdf")
  png_path <- file.path(output_dir, "tbl_data_exploration.png")

  cat("[INFO] Rendering PDF figure...\n")
  ggsave(pdf_path, composed, width = 14, height = 10, units = "in")
  cat(sprintf("[OK] Saved PDF: %s\n", pdf_path))

  cat("[INFO] Rendering PNG figure...\n")
  ggsave(png_path, composed, width = 14, height = 10, units = "in", dpi = 300)
  cat(sprintf("[OK] Saved PNG: %s\n", png_path))
}

main <- function() {
  args <- parse_args()
  df <- read_tbl(args$input)
  create_figures(df, args$output)
  cat("[DONE] Visualization completed.\n")
}

if (!interactive()) {
  tryCatch(
    main(),
    error = function(e) {
      message("[ERROR] ", e$message)
      quit(status = 1)
    }
  )
}
