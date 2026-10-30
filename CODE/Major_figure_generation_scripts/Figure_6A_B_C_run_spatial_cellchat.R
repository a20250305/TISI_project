suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(Seurat)
  library(SpatialCellChat)
  library(Matrix)
  library(data.table)
  library(future)
})


h5ad_dir <- "path/to/visium_h5ad"
db_file <- "path/to/CellChatDB.mouse_TISI_IronNiche_extended.rds"
output_dir <- "analysis_output"
sample_ids <- c(
  "20221010-1_1",
  "20221010-2_6",
  "20220401-2_8",
  "20221010-2_8"
)

major_levels <- c(
  "Tumor_cells", "Macrophages", "Endothelial_cells", "Fibroblasts",
  "T_cells", "Dendritic_cells", "B_cells"
)

subtype_to_major <- c(
  `TB-EMT` = "Tumor_cells", TBasal1 = "Tumor_cells",
  TBasal2 = "Tumor_cells", TFlike = "Tumor_cells", `TL-Alv` = "Tumor_cells",
  `TLA-EMT` = "Tumor_cells", TMlike = "Tumor_cells",
  TProliferating = "Tumor_cells",
  CHMacrophage = "Macrophages", InflMacrophage = "Macrophages",
  Spp1Macrophage = "Macrophages",
  NKcell = "T_cells", Tcell = "T_cells", CD4T = "T_cells",
  CD8T = "T_cells", Treg = "T_cells",
  Bcell = "B_cells", Plasmacell = "B_cells",
  APC = "Dendritic_cells", cDC = "Dendritic_cells", pDC = "Dendritic_cells",
  Endothelialcell = "Endothelial_cells", Fibroblast = "Fibroblasts"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(db_file))

sparse3d_to_dt <- function(x, value_name) {
  if (is.null(x)) return(NULL)
  a <- as.array(x)
  idx <- which(a != 0, arr.ind = TRUE)
  if (!nrow(idx)) return(NULL)
  dn <- dimnames(a)
  out <- data.table(
    source = dn[[1]][idx[, 1]],
    target = dn[[2]][idx[, 2]],
    pair = dn[[3]][idx[, 3]]
  )
  out[[value_name]] <- a[idx]
  out
}

run_one <- function(sample_id) {
  h5ad_file <- file.path(h5ad_dir, paste0(sample_id, ".h5ad"))
  stopifnot(file.exists(h5ad_file))

  sample_dir <- file.path(output_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  sce <- read_h5ad(h5ad_file, as = "SingleCellExperiment")
  expr <- as(assay(sce, "X"), "dgCMatrix")
  coords <- as.matrix(reducedDim(sce, "spatial"))[, 1:2, drop = FALSE]
  colnames(coords) <- c("x", "y")
  cell2loc <- as.matrix(reducedDim(sce, "cell2loc"))
  stopifnot(identical(colnames(expr), rownames(coords)))
  stopifnot(identical(colnames(expr), rownames(cell2loc)))

  seu <- CreateSeuratObject(counts = expr)
  seu <- NormalizeData(seu, verbose = FALSE)
  subtype <- colnames(cell2loc)[max.col(cell2loc, ties.method = "first")]
  seu$celltype <- unname(subtype_to_major[subtype])
  if (anyNA(seu$celltype)) stop("Unmapped cell subtype in ", sample_id)
  seu$celltype <- factor(seu$celltype, levels = major_levels)

  chat <- createSpatialCellChat(
    object = seu,
    group.by = "celltype",
    assay = "RNA",
    datatype = "spatial",
    coordinates = coords,
    spatial.factors = list(ratio = 1, tol = 65 / 2)
  )
  chat@idents <- droplevels(chat@idents)
  chat@DB <- readRDS(db_file)

  chat <- subsetData(chat)
  chat <- preProcessing(chat)
  chat <- identifyOverExpressedGenes(chat, selection.method = "meringue")
  chat <- identifyOverExpressedInteractions(chat, variable.both = FALSE)

  plan(sequential)
  options(future.globals.maxSize = 8 * 1024^3)
  chat <- computeCommunProb(
    chat,
    distance.use = TRUE,
    scale.distance = 0.2,
    contact.dependent = TRUE,
    interaction.range = 250,
    contact.range = 10
  )
  chat <- filterProbability(chat)
  chat <- filterCommunication(chat, min.links = 10, min.cells.sr = 10)

  spot_prob <- sparse3d_to_dt(chat@net$prob, "prob")
  spot_pval <- sparse3d_to_dt(chat@net$pval, "pval")
  if (!is.null(spot_prob)) {
    if (!is.null(spot_pval)) {
      spot <- merge(spot_prob, spot_pval,
                    by = c("source", "target", "pair"), all = TRUE)
    } else {
      spot <- copy(spot_prob)
      spot[, pval := NA_real_]
    }
    setnames(spot, "pair", "lr")
    fwrite(spot, file.path(sample_dir, paste0(sample_id, "_spot_LR_prob_pval.csv")))
  }
  saveRDS(chat, file.path(sample_dir, paste0(sample_id, "_chat_spot_level.rds")))

  groups <- levels(chat@idents)
  cell2loc_major <- sapply(unique(subtype_to_major[colnames(cell2loc)]), function(ct) {
    rowSums(cell2loc[, subtype_to_major[colnames(cell2loc)] == ct, drop = FALSE])
  })
  cell2loc_major <- as.matrix(cell2loc_major)
  cell2loc_use <- cell2loc_major[, groups, drop = FALSE]

  chat <- computeAvgCommunProb_Visium(
    chat,
    cell.type.decomposition = cell2loc_use,
    avg.type = "avg",
    nboot = 100,
    do.permutation = TRUE
  )
  chat <- filterCommunication(chat, min.cells = 10)
  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)

  group_lr <- as.data.table(subsetCommunication(chat))
  group_path <- as.data.table(subsetCommunication(chat, slot.name = "netP"))
  fwrite(group_lr, file.path(sample_dir, paste0(sample_id, "_group_LR_raw.csv")))
  fwrite(group_path, file.path(sample_dir, paste0(sample_id, "_group_pathway_raw.csv")))
  fwrite(group_lr[, .(
    total_weight = sum(prob, na.rm = TRUE),
    n_links = .N,
    mean_pval = mean(pval, na.rm = TRUE)
  ), by = .(source, target)],
  file.path(sample_dir, paste0(sample_id, "_group_pair_LR_summary.csv")))
  fwrite(group_path[, .(
    total_weight = sum(prob, na.rm = TRUE),
    n_links = .N,
    mean_pval = mean(pval, na.rm = TRUE)
  ), by = .(source, target)],
  file.path(sample_dir, paste0(sample_id, "_group_pair_pathway_summary.csv")))

  write.csv(chat@net$weight,
            file.path(sample_dir, paste0(sample_id, "_group_LR_weight_matrix.csv")))
  saveRDS(chat, file.path(sample_dir, paste0(sample_id, "_SpatialCellChat_final.rds")))
  invisible(chat)
}

invisible(lapply(sample_ids, run_one))
message("SpatialCellChat analysis finished: ", output_dir)
