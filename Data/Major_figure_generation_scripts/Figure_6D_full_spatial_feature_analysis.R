suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
})


args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args) >= 1) args[[1]] else "data/visium_mouse/processed"
output_dir <- if (length(args) >= 2) args[[2]] else "Figure_6D_spatial_feature_maps"
h5ad_dir <- file.path(data_root, "mouse_main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("20220401-1_3", "20220401-2_8", "20221010-2_7", "20221010-3_11")
feature_genes <- c("Ltf", "Lrp1", "Hfe", "Tfrc")

major_map <- c(
  `TB-EMT` = "Tumor_cells", TBasal1 = "Tumor_cells", TBasal2 = "Tumor_cells",
  TFlike = "Tumor_cells", `TL-Alv` = "Tumor_cells", `TLA-EMT` = "Tumor_cells",
  TMlike = "Tumor_cells", TProliferating = "Tumor_cells",
  CHMacrophage = "Macrophages", InflMacrophage = "Macrophages",
  Spp1Macrophage = "Macrophages", NKcell = "NK_cells",
  Tcell = "T_cells", CD4T = "T_cells", CD8T = "T_cells", Treg = "T_cells",
  Bcell = "B_cells", Plasmacell = "B_cells", APC = "Dendritic_cells",
  cDC = "Dendritic_cells", pDC = "Dendritic_cells",
  Endothelialcell = "Endothelial_cells", Fibroblast = "Fibroblasts"
)

major_colors <- c(
  Tumor_cells = "#D73027", Macrophages = "#3288BD", Fibroblasts = "#66BD63",
  Endothelial_cells = "#91BFDB", T_cells = "#FDAE61", B_cells = "#FEE08B",
  NK_cells = "#F46D43", Dendritic_cells = "#8073AC"
)


read_csr_genes <- function(h5ad_file, matrix_path, genes, gene_symbols) {
  shape <- h5readAttributes(h5ad_file, matrix_path)$shape
  n_obs <- as.integer(shape[[1]])
  data <- h5read(h5ad_file, paste0(matrix_path, "/data"))
  indices <- h5read(h5ad_file, paste0(matrix_path, "/indices"))
  indptr <- h5read(h5ad_file, paste0(matrix_path, "/indptr"))
  wanted <- match(genes, gene_symbols)
  out <- matrix(0, nrow = n_obs, ncol = length(genes),
                dimnames = list(NULL, genes))
  for (i in seq_len(n_obs)) {
    a <- indptr[[i]] + 1L
    b <- indptr[[i + 1L]]
    if (a > b) next
    positions <- a:b
    hit <- match(indices[positions] + 1L, wanted, nomatch = 0L)
    keep <- which(hit > 0L)
    if (length(keep)) out[i, hit[keep]] <- data[positions[keep]]
  }
  out
}

read_cell2loc <- function(h5ad_file, n_obs) {
  listing <- h5ls(h5ad_file, recursive = TRUE)
  subtypes <- listing$name[listing$group == "/obsm/cell2loc" &
                            listing$otype == "H5I_DATASET" &
                            listing$name != "_index"]
  stopifnot(length(subtypes) > 0L)
  out <- sapply(subtypes, function(x) {
    as.numeric(h5read(h5ad_file, paste0("/obsm/cell2loc/", x)))
  })
  out <- as.matrix(out)
  colnames(out) <- subtypes
  stopifnot(nrow(out) == n_obs)
  out
}

read_one_sample <- function(sample_id) {
  h5ad_file <- file.path(h5ad_dir, paste0(sample_id, ".h5ad"))
  stopifnot(file.exists(h5ad_file))

  obs_id <- as.character(h5read(h5ad_file, "/obs/_index"))
  spatial <- h5read(h5ad_file, "/obsm/spatial")
  if (nrow(spatial) == 2L) spatial <- t(spatial)
  spatial <- as.data.table(spatial)
  setnames(spatial, c("x", "y"))
  spatial[, spot_id_raw := obs_id]

  gene_symbols <- as.character(h5read(h5ad_file, "/var/gene_symbols"))
  expr_norm <- read_csr_genes(h5ad_file, "/X", feature_genes, gene_symbols)
  expr_counts <- read_csr_genes(h5ad_file, "/raw/X", feature_genes, gene_symbols)

  cell2loc <- read_cell2loc(h5ad_file, length(obs_id))
  dominant_subtype <- colnames(cell2loc)[max.col(cell2loc, ties.method = "first")]
  major <- unname(major_map[dominant_subtype])
  major[is.na(major)] <- "Unmapped"

  out <- data.table(
    sample = sample_id,
    spot_id = sub("-[0-9]+$", "", obs_id),
    spot_id_raw = obs_id,
    x = spatial$x,
    y = spatial$y,
    dominant_subtype = dominant_subtype,
    major_cell_type = major
  )
  out <- cbind(out, as.data.table(expr_norm))
  setnames(out, feature_genes, paste0(feature_genes, "_norm"))
  out <- cbind(out, as.data.table(expr_counts))
  setnames(out, feature_genes, paste0(feature_genes, "_count"))
  out
}

source_list <- lapply(samples, function(s) {
  message("Reading ", s)
  read_one_sample(s)
})
source_dt <- rbindlist(source_list, use.names = TRUE)

fwrite(source_dt, file.path(output_dir, "Figure_6D_source_data_all_MRD.csv"))
for (sample_id in samples) {
  fwrite(source_dt[sample == sample_id],
         file.path(output_dir, paste0("Figure_6D_source_data_", sample_id, ".csv")))
}

gene_manifest <- data.table(
  gene = feature_genes,
  normalized_column = paste0(feature_genes, "_norm"),
  count_column = paste0(feature_genes, "_count"),
  matrix_used_for_plot = "/X",
  raw_count_matrix = "/raw/X"
)
fwrite(gene_manifest, file.path(output_dir, "Figure_6D_gene_manifest.csv"))

sample_summary <- source_dt[, .(
  n_spots = .N,
  n_major_mapped = sum(major_cell_type != "Unmapped"),
  n_tumor_cells = sum(major_cell_type == "Tumor_cells"),
  n_macrophages = sum(major_cell_type == "Macrophages"),
  n_fibroblasts = sum(major_cell_type == "Fibroblasts")
), by = sample]
fwrite(sample_summary, file.path(output_dir, "Figure_6D_run_summary.csv"))

fwrite(data.table(
  source = c("/X", "/raw/X", "/obsm/spatial", "/obsm/cell2loc"),
  use = c("normalized expression used for feature maps", "raw counts retained for audit",
          "Visium spatial coordinates", "cell2location subtype abundances used for dominant labels")
), file.path(output_dir, "Figure_6D_input_manifest.csv"))

fwrite(data.table(
  panel_row = c("upper", "lower"),
  sample = c("20221010-2_7", "20220401-1_3"),
  features = c("Ltf; Lrp1; Major cell types", "Hfe; Tfrc; Major cell types")
), file.path(output_dir, "Figure_6D_paper_panel_sample_mapping.csv"))

