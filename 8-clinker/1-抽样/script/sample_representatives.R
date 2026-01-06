#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phytools)
})

option_list <- list(
  make_option(c("--tree"), type = "character", help = "Path to Newick tree"),
  make_option(c("--meta"), type = "character", help = "Path to meta TSV"),
  make_option(c("--outdir"), type = "character", help = "Output directory"),
  make_option(c("--n"), type = "integer", default = 50, help = "Total representatives"),
  make_option(c("--n_clade"), type = "integer", default = 18, help = "Min reps in Rv1819c clade"),
  make_option(c("--seed"), type = "integer", default = 123, help = "Random seed"),
  make_option(c("--rv_tip"), type = "character", help = "Target tip label (e.g., Rv1819c)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$tree) || is.null(opt$meta) || is.null(opt$outdir)) {
  stop("--tree, --meta, --outdir are required")
}
if (is.null(opt$rv_tip)) {
  stop("--rv_tip is required")
}

set.seed(opt$seed)

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

rv_tip <- opt$rv_tip

tree <- read.tree(opt$tree)
if (is.null(tree$tip.label) || length(tree$tip.label) == 0) {
  stop("Tree has no tip labels")
}
if (!rv_tip %in% tree$tip.label) {
  stop("Rv1819c not found in tree tip labels")
}

tips <- tree$tip.label

meta <- read.delim(opt$meta, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(meta) == 0) stop("meta is empty")

norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

find_id_col <- function(cols) {
  candidates <- c("sample", "sampleid", "id", "label", "tip", "tiplabel", "tiplab")
  ncols <- norm_name(cols)
  for (cand in candidates) {
    hit <- which(ncols == cand)
    if (length(hit) > 0) return(cols[hit[1]])
  }
  for (cand in candidates) {
    hit <- which(grepl(cand, ncols, fixed = TRUE))
    if (length(hit) > 0) return(cols[hit[1]])
  }
  return(NA_character_)
}

id_col <- find_id_col(colnames(meta))
if (is.na(id_col)) {
  stop(paste0("Cannot detect sample ID column. Available columns: ", paste(colnames(meta), collapse = ", ")))
}

rank_syn <- list(
  domain = c("domain", "superkingdom", "kingdom", "d", "d__"),
  phylum = c("phylum", "p", "p__"),
  class = c("class", "c", "c__"),
  order = c("order", "o", "o__"),
  family = c("family", "f", "f__"),
  genus = c("genus", "g", "g__"),
  species = c("species", "sp", "s", "s__")
)

find_rank_col <- function(cols, rank) {
  ncols <- norm_name(cols)
  for (syn in rank_syn[[rank]]) {
    hit <- which(ncols == syn)
    if (length(hit) > 0) return(cols[hit[1]])
  }
  for (syn in rank_syn[[rank]]) {
    hit <- which(grepl(syn, ncols, fixed = TRUE))
    if (length(hit) > 0) return(cols[hit[1]])
  }
  return(NA_character_)
}

rank_cols <- sapply(names(rank_syn), function(r) find_rank_col(colnames(meta), r), USE.NAMES = TRUE)

idx <- match(tips, meta[[id_col]])
meta_missing_count <- sum(is.na(idx))

aligned <- data.frame(tip_label = tips, stringsAsFactors = FALSE)
for (r in names(rank_cols)) {
  col <- rank_cols[[r]]
  if (!is.na(col)) {
    aligned[[r]] <- meta[[col]][idx]
  } else {
    aligned[[r]] <- NA_character_
  }
}

make_species_key <- function(df) {
  if (all(is.na(df$species))) {
    if (!all(is.na(df$genus))) {
      return(df$genus)
    }
    last_rank <- c("genus", "family", "order", "class", "phylum", "domain")
    for (r in last_rank) {
      if (!all(is.na(df[[r]]))) return(df[[r]])
    }
    return(rep(NA_character_, nrow(df)))
  }
  return(df$species)
}

species_key <- make_species_key(aligned)

maxmin_pick <- function(candidates, n, dist_mat, selected, key_primary = NULL, key_secondary = NULL, relax_counts) {
  used_primary <- character(0)
  used_secondary <- character(0)
  if (!is.null(key_primary)) used_primary <- unique(na.omit(key_primary[selected]))
  if (!is.null(key_secondary)) used_secondary <- unique(na.omit(key_secondary[selected]))

  while (length(selected) < n && length(candidates) > 0) {
    pool <- setdiff(candidates, selected)
    if (length(pool) == 0) break

    if (!is.null(key_primary)) {
      pool1 <- pool[is.na(key_primary[pool]) | !(key_primary[pool] %in% used_primary)]
      if (length(pool1) > 0) {
        pool <- pool1
      } else {
        relax_counts$primary <- relax_counts$primary + 1
      }
    }

    if (!is.null(key_secondary)) {
      pool2 <- pool[is.na(key_secondary[pool]) | !(key_secondary[pool] %in% used_secondary)]
      if (length(pool2) > 0) {
        pool <- pool2
      } else {
        relax_counts$secondary <- relax_counts$secondary + 1
      }
    }

    dmin <- sapply(pool, function(t) {
      if (length(selected) == 0) return(Inf)
      min(dist_mat[t, selected, drop = TRUE])
    })
    max_d <- max(dmin)
    picks <- names(dmin)[dmin == max_d]
    pick <- sample(picks, 1)
    selected <- c(selected, pick)
    if (!is.null(key_primary) && !is.na(key_primary[pick])) used_primary <- unique(c(used_primary, key_primary[pick]))
    if (!is.null(key_secondary) && !is.na(key_secondary[pick])) used_secondary <- unique(c(used_secondary, key_secondary[pick]))
  }

  list(selected = selected, relax_counts = relax_counts)
}

# Define Rv1819c clade
rv_idx <- which(tips == rv_tip)
parent_node <- tree$edge[tree$edge[, 2] == rv_idx, 1]

get_parent <- function(node) {
  hit <- tree$edge[tree$edge[, 2] == node, 1]
  if (length(hit) == 0) return(NA_integer_)
  hit[1]
}

clade_tips <- extract.clade(tree, parent_node)$tip.label
if (length(clade_tips) < opt$n_clade) {
  current <- parent_node
  while (length(clade_tips) < opt$n_clade) {
    current <- get_parent(current)
    if (is.na(current)) break
    clade_tips <- extract.clade(tree, current)$tip.label
  }
  warning("Rv1819c clade size < n_clade, moved up to larger clade")
}

n_in_clade <- min(opt$n_clade, length(clade_tips))

dist_mat <- as.matrix(cophenetic(tree))

relax_species_clade <- list(primary = 0, secondary = 0)
clade_sel <- maxmin_pick(
  candidates = clade_tips,
  n = n_in_clade,
  dist_mat = dist_mat,
  selected = rv_tip,
  key_primary = species_key,
  key_secondary = NULL,
  relax_counts = relax_species_clade
)

clade_selected <- unique(clade_sel$selected)
relax_species_clade <- clade_sel$relax_counts

# Outgroup sampling
out_candidates <- setdiff(tips, clade_tips)

n_out <- opt$n - length(clade_selected)
if (n_out < 0) n_out <- 0
n_out <- min(n_out, length(out_candidates))

relax_genus_out <- list(primary = 0, secondary = 0)

out_sel <- maxmin_pick(
  candidates = out_candidates,
  n = length(clade_selected) + n_out,
  dist_mat = dist_mat,
  selected = clade_selected,
  key_primary = aligned$genus,
  key_secondary = species_key,
  relax_counts = relax_genus_out
)

out_selected <- out_sel$selected
relax_genus_out <- out_sel$relax_counts

out_only <- setdiff(out_selected, clade_selected)
selected_all <- unique(c(clade_selected, out_only))

selection_stage <- rep(NA_character_, length(tips))
selection_stage[tips %in% clade_selected] <- "clade"
selection_stage[tips %in% out_only] <- "outgroup"
selection_stage[tips == rv_tip] <- "forced_Rv1819c"

in_clade <- tips %in% clade_tips

distance_to_rv <- dist_mat[, rv_tip]

out_df <- data.frame(
  tip_label = tips,
  in_rv1819c_clade = in_clade,
  domain = aligned$domain,
  phylum = aligned$phylum,
  class = aligned$class,
  order = aligned$order,
  family = aligned$family,
  genus = aligned$genus,
  species = aligned$species,
  selection_stage = selection_stage,
  distance_to_Rv1819c = distance_to_rv,
  stringsAsFactors = FALSE
)

out_df <- out_df[out_df$tip_label %in% selected_all, , drop = FALSE]

write.table(out_df, file = file.path(opt$outdir, "representatives.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

report <- c(
  paste0("seed=", opt$seed),
  paste0("n_total=", opt$n),
  paste0("n_clade_target=", opt$n_clade),
  paste0("clade_tips=", length(clade_tips)),
  paste0("n_in_clade=", length(clade_selected)),
  paste0("n_out=", length(out_only)),
  paste0("relax_species_clade=", relax_species_clade$primary),
  paste0("relax_genus_out=", relax_genus_out$primary),
  paste0("relax_species_out=", relax_genus_out$secondary),
  paste0("meta_missing_tip=", meta_missing_count)
)
writeLines(report, con = file.path(opt$outdir, "sampling_report.txt"))

# Optional tree plot
pdf(file.path(opt$outdir, "representatives_tree.pdf"), width = 8, height = 10)
plot(tree, cex = 0.4, no.margin = TRUE)
sel_tips <- which(tips %in% selected_all)
rv_tips <- which(tips == rv_tip)
tiplabels(pch = 19, col = "red", tip = sel_tips, cex = 0.6)
tiplabels(pch = 19, col = "blue", tip = rv_tips, cex = 0.8)
dev.off()
