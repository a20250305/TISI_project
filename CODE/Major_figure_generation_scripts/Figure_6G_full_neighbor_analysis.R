if (!nzchar(Sys.getenv("BASILISK_EXTERNAL_DIR"))) {
  Sys.setenv(BASILISK_EXTERNAL_DIR = file.path(tempdir(), "basilisk"))
}

suppressPackageStartupMessages({
  library(data.table)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
  library(FNN)
})

args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args) >= 1) args[[1]] else "data/visium_mouse/processed"
output_dir <- if (length(args) >= 2) args[[2]] else "Figure_6G_full_neighbor_analysis"

h5ad_dir <- file.path(data_root, "mouse_main")
tisi_dir <- file.path(data_root, "TISI_scores_visum_011626")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("20220401-1_3", "20220401-2_8", "20221010-2_7", "20221010-3_11")
k_neighbors <- 6L

macro_genes <- c(
  "Slc7a11", "Slc11a2", "Tfrc", "Hfe", "Steap3", "Ncoa4", "Fth1", "Alas1",
  "Hmox1", "Blvrb", "Blvra", "Hmox2", "Slc48a1", "Lrp1", "Flvcr1", "Bach1",
  "Col1a2", "Mmp3", "Fn1", "Col1a1"
)

fibroblast_genes <- c(
  "Hfe", "Steap3", "Tfrc", "Fech", "Hmox1", "Blvrb", "Hmox2", "Lrp1",
  "Col1a2", "Mmp3", "Fn1", "Col1a1", "Selenbp1", "C3", "Lif"
)

tumor_subtypes <- c(
  "TB-EMT", "TBasal1", "TBasal2", "TFlike", "TL-Alv", "TLA-EMT",
  "TMlike", "TProliferating"
)

clean_barcode <- function(x) {
  sub("-[0-9]+-[0-9]+$", "", x)
}

read_tisi <- function(sample_id) {
  file <- file.path(tisi_dir, paste0(sample_id, "_TISI_scores.csv"))
  stopifnot(file.exists(file))
  dt <- fread(file)
  stopifnot(all(c("spot_id", "TISI_scaled") %in% names(dt)))
  dt[, spot := clean_barcode(spot_id)]
  dt[, .(TISI_scaled = mean(TISI_scaled, na.rm = TRUE)), by = spot]
}

make_neighbor_table <- function(sample_id, sce) {
  expr <- assay(sce, assayNames(sce)[1])
  colnames(expr) <- clean_barcode(colnames(expr))

  spatial <- as.matrix(reducedDim(sce, "spatial"))
  colnames(spatial) <- c("x", "y")
  if (is.null(rownames(spatial))) rownames(spatial) <- colnames(sce)
  rownames(spatial) <- clean_barcode(rownames(spatial))

  cell2loc <- as.matrix(reducedDim(sce, "cell2loc"))
  if (is.null(rownames(cell2loc))) rownames(cell2loc) <- colnames(sce)
  rownames(cell2loc) <- clean_barcode(rownames(cell2loc))

  common_spots <- Reduce(intersect, list(rownames(spatial), rownames(cell2loc), colnames(expr)))
  stopifnot(length(common_spots) > k_neighbors)
  spatial <- spatial[common_spots, , drop = FALSE]
  cell2loc <- cell2loc[common_spots, , drop = FALSE]
  expr <- expr[, common_spots, drop = FALSE]

  dominant_subtype <- colnames(cell2loc)[max.col(cell2loc, ties.method = "first")]
  # This is the original Figure 6G convention: tumor subtypes are collapsed,
  # whereas non-tumor subtypes retain their original labels.  In particular,
  # the fibroblast-neighbor filter below therefore uses "Fibroblast".
  celltype <- dominant_subtype
  celltype[dominant_subtype %in% tumor_subtypes] <- "Tumor_cells"

  nn <- get.knn(spatial, k = k_neighbors)$nn.index

  niche <- rep(NA_character_, length(common_spots))
  if ("chr_aa" %in% reducedDimNames(sce)) {
    chr_aa <- as.matrix(reducedDim(sce, "chr_aa"))
    if (is.null(rownames(chr_aa))) rownames(chr_aa) <- colnames(sce)
    rownames(chr_aa) <- clean_barcode(rownames(chr_aa))
    chr_aa <- chr_aa[common_spots, , drop = FALSE]
    niche <- paste0("Niche_", max.col(chr_aa, ties.method = "first"))
  }

  neighbor_table <- data.table(
    spot = common_spots,
    celltype = unname(celltype),
    niche = unname(niche)
  )
  for (i in seq_len(k_neighbors)) {
    neighbor_table[[paste0("nbr", i, "_spot")]] <- common_spots[nn[, i]]
    neighbor_table[[paste0("nbr", i, "_celltype")]] <- celltype[nn[, i]]
    neighbor_table[[paste0("nbr", i, "_niche")]] <- niche[nn[, i]]
  }

  tumor_idx <- which(celltype == "Tumor_cells")
  stopifnot(length(tumor_idx) > 0)

  list(
    expr = expr,
    celltype = celltype,
    nn = nn,
    common_spots = common_spots,
    neighbor_table = neighbor_table,
    tumor_idx = tumor_idx
  )
}

neighbor_matrix <- function(expr, celltype, nn, tumor_idx, genes, mode) {
  genes <- intersect(genes, rownames(expr))
  out <- lapply(tumor_idx, function(i) {
    neighbor_idx <- nn[i, ]
    if (mode == "fibroblasts") {
      neighbor_idx <- neighbor_idx[celltype[neighbor_idx] == "Fibroblast"]
    }
    neighbor_spots <- colnames(expr)[neighbor_idx]
    if (length(neighbor_spots) == 0) {
      return(setNames(rep(NA_real_, length(genes)), genes))
    }
    # The published Figure 6G macrophage block uses all six neighbors.
    # The fibroblast block uses only neighbors labelled Fibroblast.
    as.numeric(rowMeans(expr[genes, neighbor_spots, drop = FALSE])) |>
      setNames(genes)
  })
  m <- do.call(cbind, out)
  if (is.null(dim(m))) m <- matrix(m, nrow = length(genes))
  rownames(m) <- genes
  colnames(m) <- colnames(expr)[tumor_idx]
  m
}

make_figure_matrix <- function(sample_id, info, mode, genes) {
  tisi <- read_tisi(sample_id)
  tumor_spots <- info$common_spots[info$tumor_idx]
  tisi <- tisi[match(tumor_spots, tisi$spot)]
  stopifnot(all(tisi$spot == tumor_spots))

  m <- neighbor_matrix(
    expr = info$expr,
    celltype = info$celltype,
    nn = info$nn,
    tumor_idx = info$tumor_idx,
    genes = genes,
    mode = mode
  )

  dt <- data.table(spot = colnames(m), TISI_scaled = tisi$TISI_scaled)
  dt <- cbind(dt, as.data.table(t(m)))
  dt[, sample := sample_id]
  setorder(dt, TISI_scaled)
  dt
}

summaries <- list()
macro_list <- list()
fibro_list <- list()

for (sample_id in samples) {
  h5ad_file <- file.path(h5ad_dir, paste0(sample_id, ".h5ad"))
  stopifnot(file.exists(h5ad_file))
  message("Processing ", sample_id)
  sce <- readH5AD(h5ad_file)
  info <- make_neighbor_table(sample_id, sce)

  fwrite(
    info$neighbor_table,
    file.path(output_dir, paste0("Figure_6G_raw_", sample_id, "_spot_neighbor_table.csv"))
  )

  tumor_neighbors <- rbindlist(lapply(info$tumor_idx, function(i) {
    data.table(
      sample = sample_id,
      tumor_spot = info$common_spots[i],
      neighbor_index = seq_len(k_neighbors),
      neighbor_spot = info$common_spots[info$nn[i, ]],
      neighbor_celltype = info$celltype[info$nn[i, ]],
      neighbor_rank = seq_len(k_neighbors)
    )
  }))
  fwrite(
    tumor_neighbors,
    file.path(output_dir, paste0("Figure_6G_", sample_id, "_tumor_neighbor_assignments.csv"))
  )

  macro_dt <- make_figure_matrix(sample_id, info, "macrophages", macro_genes)
  fibro_dt <- make_figure_matrix(sample_id, info, "fibroblasts", fibroblast_genes)
  fwrite(macro_dt, file.path(output_dir, paste0("Figure_6G_", sample_id, "_macrophages_matrix.csv")))
  fwrite(fibro_dt, file.path(output_dir, paste0("Figure_6G_", sample_id, "_fibroblasts_matrix.csv")))
  macro_list[[sample_id]] <- macro_dt
  fibro_list[[sample_id]] <- fibro_dt

  fibro_genes_present <- intersect(fibroblast_genes, names(fibro_dt))

  summaries[[sample_id]] <- data.table(
    sample = sample_id,
    total_spots = length(info$common_spots),
    tumor_dominant_spots = length(info$tumor_idx),
    macrophage_matrix_spots = nrow(macro_dt),
    fibroblast_matrix_spots = nrow(fibro_dt),
    fibroblast_nonmissing_spots = sum(rowSums(!is.na(fibro_dt[, ..fibro_genes_present])) > 0),
    macro_genes_available = sum(macro_genes %in% rownames(info$expr)),
    fibroblast_genes_available = sum(fibroblast_genes %in% rownames(info$expr))
  )
}

macro_all <- rbindlist(macro_list, fill = TRUE)
fibro_all <- rbindlist(fibro_list, fill = TRUE)
setorder(macro_all, TISI_scaled)
setorder(fibro_all, TISI_scaled)

macro_all <- macro_all[, c("spot", "sample", "TISI_scaled", macro_genes), with = FALSE]
fibro_all <- fibro_all[, c("spot", "sample", "TISI_scaled", fibroblast_genes), with = FALSE]

fwrite(macro_all, file.path(output_dir, "Figure_6G_source_data_macrophages.csv"))
fwrite(fibro_all, file.path(output_dir, "Figure_6G_source_data_fibroblasts.csv"))

manifest <- rbindlist(list(
  data.table(Module = "Macrophages", Gene = macro_genes, GeneOrder = seq_along(macro_genes)),
  data.table(Module = "Fibroblasts", Gene = fibroblast_genes, GeneOrder = seq_along(fibroblast_genes))
))
fwrite(manifest, file.path(output_dir, "Figure_6G_gene_manifest.csv"))
fwrite(rbindlist(summaries), file.path(output_dir, "Figure_6G_run_summary.csv"))

cat("Completed Figure 6G neighbor analysis.\n")
cat("Output directory: ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")
