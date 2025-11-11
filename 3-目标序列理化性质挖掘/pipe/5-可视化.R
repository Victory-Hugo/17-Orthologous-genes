#!/usr/bin/env Rscript
# Author: BigLin
# Requirements: R (base), capable of running hist/density; no extra packages needed.

input_path <- "/mnt/c/Users/Administrator/Desktop/output/merged_physicochemical.tsv"
output_pdf <- "/mnt/c/Users/Administrator/Desktop/output/merged_physicochemical.pdf"

if (!file.exists(input_path)) {
  stop(sprintf("找不到输入文件: %s", input_path))
}

data <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
colors <- c("#83AEF0", "#7FD6FF", "#90D1D5", "#BC8A8A")
numeric_cols <- setdiff(names(data), c("Source_Fasta", "Sequence_ID"))
numeric_cols <- numeric_cols[sapply(numeric_cols, function(col) {
  suppressWarnings(any(is.finite(as.numeric(data[[col]]))))
})]

if (length(numeric_cols) == 0) {
  stop("没有可用的数值列用于绘图。")
}

rows <- ceiling(sqrt(length(numeric_cols)))
cols <- ceiling(length(numeric_cols) / rows)

pdf(output_pdf, width = 3.5 * cols, height = 3.5 * rows)
par(mfrow = c(rows, cols), mar = c(4, 4, 2, 1))

for (idx in seq_along(numeric_cols)) {
  col <- numeric_cols[idx]
  values <- suppressWarnings(as.numeric(data[[col]]))
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    plot.new()
    title(main = sprintf("%s\n(无有效数值)", col), cex.main = 0.9)
    next
  }

  color_idx <- ((idx - 1) %% length(colors)) + 1
  hist(values,
       main = col,
       xlab = col,
       col = adjustcolor(colors[color_idx], alpha.f = 0.6),
       border = colors[color_idx],
       freq = TRUE)
}

dev.off()

message(sprintf("直方图 PDF 已生成: %s", output_pdf))
