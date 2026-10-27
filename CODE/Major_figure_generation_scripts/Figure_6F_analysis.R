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
output_dir <- if (length(args) >= 2) args[[2]] else "Figure_6F_full_neighbor_analysis"

h5ad_dir <- file.path(data_root, "mouse_main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("20220401-1_3", "20220401-2_8", "20221010-2_7", "20221010-3_11")
k_neighbors <- 6L
n_permutations <- 1000L
set.seed(1L)

tumor_subtypes <- c(
  "TB-EMT", "TBasal1", "TBasal2", "TFlike", "TL-Alv", "TLA-EMT",
  "TMlike", "TProliferating"
)

major_map <- c(
  `TB-EMT` = "Tumor_cells", TBasal1 = "Tumor_cells", TBasal2 = "Tumor_cells",
  TFlike = "Tumor_cells", `TL-Alv` = "Tumor_cells", `TLA-EMT` = "Tumor_cells",
  TMlike = "Tumor_cells", TProliferating = "Tumor_cells",
  CHMacrophage = "Macrophages", InflMacrophage = "Macrophages", Spp1Macrophage = "Macrophages",
  NKcell = "NK_cells", Tcell = "T_cells", CD4T = "T_cells", CD8T = "T_cells", Treg = "T_cells",
  Bcell = "B_cells", Plasmacell = "B_cells",
  APC = "Dendritic_cells", cDC = "Dendritic_cells", pDC = "Dendritic_cells",
  Endothelialcell = "Endothelial_cells", Fibroblast = "Fibroblasts"
)

major_order <- c(
  "Tumor_cells", "Macrophages", "NK_cells", "T_cells", "B_cells",
  "Dendritic_cells", "Endothelial_cells", "Fibroblasts"
)

gene_sets <- list(
  iron_uptake = c("Trf", "Tfrc", "Tfr2", "Slc11a2", "Slc39a14", "Slc39a8", "Slc7a11", "Slc25a37", "Slc25a28"),
  iron_storage = c("Fth1", "Ftl", "Ftl1", "Slc40a1", "Hamp", "Hfe", "Hjv", "Steap2", "Steap3", "Cybrd1", "Ncoa4", "Aco1", "Ireb2", "Cp", "Heph", "Ftmt"),
  heme_globin = c("Hmox1", "Hmox2", "Blvrb", "Blvra", "Hp", "Hpr", "Lrp1", "Lrp2", "Alas1", "Alas2", "Fech", "Uros", "Urod", "Hmbs", "Cpox", "Ppox", "Abcb6", "Abcb10", "Abcg2", "Flvcr1", "Flvcr2", "Slc48a1", "Slc46a1", "Hba-a1", "Hba-a2", "Hbb", "Hbb-bs", "Hbb-bt", "Hbq1", "Hbd", "Hbz", "Ahsp", "Bpgm", "Alad"),
  redox = c("Gpx4", "Aifm2", "Nfe2l2", "Keap1", "Gclc", "Gclm", "Slc3a2", "Sat1", "Acsl4", "Lpcat3", "Sod2", "Prdx2", "Txnrd1", "Cat", "Gss", "Gch1", "Coq2"),
  iMac = c("Lcn2", "Slc22a17", "Vcam1", "Cd163", "Siglec1", "Ccr3", "Pparg", "Maf", "Spic", "Mpo", "Csf1r", "Cybb", "Nox1", "Nox4"),
  iron_SASP = c("Il6", "Il11", "Lif", "Cxcl1", "Ccl2", "Ccl12", "Serpine1", "Timp1", "Il1b", "Stat3", "Csf1", "Il34", "Ccl5"),
  fibrosis = c("Col1a1", "Col1a2", "Fn1", "Acta2", "Mmp3", "Mmp8", "Mmp9", "Vegfa", "Pecam1", "Ednrb", "Tns1", "Bmp2k"),
  iron_TME = c("Hmox1", "Slc40a1", "Fth1", "Ftl1", "Tfrc", "Lcn2", "Ccl2", "Kdm6b", "Ednrb", "Bach1", "Stat3", "Vegfa", "Pvr", "Pdcd1"),
  metal_transport_extra = c("Slc30a1", "Slc30a2", "Slc30a3", "Slc30a4", "Slc30a5", "Slc30a6", "Slc30a7", "Slc30a8", "Slc30a10", "Slc39a1", "Slc39a2", "Slc39a4", "Slc39a5", "Slc39a6", "Slc39a10", "Slc39a11", "Slc39a12", "Slc39a13", "Slc26a3", "Slc26a6", "Slc26a9", "Slc31a1", "Atp7a", "Atp7b", "Atox1", "Heph", "Hephl1", "Atp1a1", "Atp1a2", "Atp1a3", "Atp1b1", "Atp2b1", "Atp2b2", "Atp2b3", "Atp2b4", "Trpm6", "Trpm7", "Trpv6", "Slc4a1", "Slc5a1", "Slc6a8", "Slc6a9", "Slc8a1", "Slc8a2", "Slc8a3", "Slc9a3"),
  erythroid_RBC = c("Gata1", "Klf1", "Klf3", "Tal1", "Lmo2", "Epor", "Ermap", "Gypa", "Gypb", "Gypc", "Gype", "Rhag", "Rhce", "Rhd", "Ank1", "Spta1", "Sptb", "Mpp1", "Icam4", "Kel", "Bcam", "Xk"),
  autophagy_ferroptosis = c("Atg4a", "Atg5", "Atg7", "Map1lc3a", "Map1lc3b", "Map1lc3c", "Bnip3l", "Optn", "Cisd1", "Prnp", "Chmp5", "Chmp6"),
  redox_extended = c("Hspb1", "Vdac2", "Vdac3", "Pcbp1", "Pcbp2", "Akr1c1", "Akr1c2", "Akr1c3", "Smox", "Ucp2"),
  metabolism_misc = c("Acsl1", "Acsl3", "Acsl5", "Acsl6", "Fdft1", "Hmgcr", "Por", "Pc", "Ftcd", "Lpin2", "Mboat2", "Agpat4", "Gde1", "Gmps", "Minpp1"),
  immune_misc = c("C3", "Acp5", "Dpp4", "Selenbp1", "Sidt2", "Igsf3", "Mospd1", "Mfhas1", "Hdgf"),
  others = c("Ltf", "Cd44", "Hrg1")
)

gene_order <- unique(unlist(gene_sets, use.names = FALSE))

clean_barcode <- function(x) {
  sub("-[0-9]+-[0-9]+$", "", x)
}

neighbor_enrich <- function(target_idx, labels, nn, nperm) {
  types <- sort(unique(labels))
  observed <- setNames(numeric(length(types)), types)
  for (i in target_idx) {
    observed[labels[nn[i, ]]] <- observed[labels[nn[i, ]]] + 1
  }

  permuted <- matrix(0, nrow = nperm, ncol = length(types), dimnames = list(NULL, types))
  for (p in seq_len(nperm)) {
    perm_labels <- sample(labels)
    for (i in target_idx) {
      permuted[p, perm_labels[nn[i, ]]] <- permuted[p, perm_labels[nn[i, ]]] + 1
    }
  }

  mu <- colMeans(permuted)
  sdv <- apply(permuted, 2, sd)
  z <- (observed - mu) / sdv
  p <- vapply(seq_along(types), function(i) {
    mean(abs(permuted[, i] - mu[i]) >= abs(observed[i] - mu[i]))
  }, numeric(1))

  data.table(
    CellType = types,
    Observed = as.numeric(observed),
    Z = as.numeric(z),
    P = p,
    FDR = p.adjust(p, "fdr")
  )
}

weighted_neighbor_expression <- function(spot_nb, cell_type, genes) {
  weights <- spot_nb[[cell_type]]
  keep <- is.finite(weights) & weights > 0
  if (!any(keep)) {
    out <- as.list(setNames(rep(NA_real_, length(genes)), genes))
    out[["Weight"]] <- cell_type
    return(as.data.table(out))
  }

  out <- spot_nb[keep, lapply(.SD, function(x) {
    sum(x * weights[keep], na.rm = TRUE) / sum(weights[keep], na.rm = TRUE)
  }), .SDcols = genes]
  out[, Weight := cell_type]
  out
}

run_one_sample <- function(sample_id) {
  h5ad_file <- file.path(h5ad_dir, paste0(sample_id, ".h5ad"))
  stopifnot(file.exists(h5ad_file))

  message("Processing ", sample_id)
  sample_out <- output_dir
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)

  sce <- readH5AD(h5ad_file)
  expr <- assay(sce, assayNames(sce)[1])
  colnames(expr) <- clean_barcode(colnames(expr))

  spatial <- reducedDim(sce, "spatial")
  spatial <- as.matrix(spatial)
  colnames(spatial) <- c("x", "y")
  if (is.null(rownames(spatial))) rownames(spatial) <- colnames(sce)
  rownames(spatial) <- clean_barcode(rownames(spatial))

  cell2loc <- as.matrix(reducedDim(sce, "cell2loc"))
  if (is.null(rownames(cell2loc))) rownames(cell2loc) <- colnames(sce)
  rownames(cell2loc) <- clean_barcode(rownames(cell2loc))

  chr_aa <- as.matrix(reducedDim(sce, "chr_aa"))
  rownames(chr_aa) <- rownames(cell2loc)
  colnames(chr_aa) <- paste0("Niche_", seq_len(ncol(chr_aa)))

  common_spots <- Reduce(intersect, list(rownames(spatial), rownames(cell2loc), colnames(expr)))
  stopifnot(length(common_spots) > 0)
  spatial <- spatial[common_spots, , drop = FALSE]
  cell2loc <- cell2loc[common_spots, , drop = FALSE]
  chr_aa <- chr_aa[common_spots, , drop = FALSE]
  expr <- expr[, common_spots, drop = FALSE]

  dominant_subtype <- colnames(cell2loc)[max.col(cell2loc, ties.method = "first")]
  celltype_sub <- dominant_subtype
  celltype_sub[dominant_subtype %in% tumor_subtypes] <- "Tumor_cells"
  tumor_idx <- which(celltype_sub == "Tumor_cells")
  stopifnot(length(tumor_idx) >= 10)

  nn <- get.knn(spatial, k = k_neighbors)$nn.index
  dominant_niche <- paste0("Niche_", max.col(chr_aa, ties.method = "first"))

  neighbor_table <- data.table(
    spot = common_spots,
    celltype = unname(celltype_sub),
    niche = dominant_niche
  )
  for (i in seq_len(k_neighbors)) {
    neighbor_table[[paste0("nbr", i, "_spot")]] <- common_spots[nn[, i]]
    neighbor_table[[paste0("nbr", i, "_celltype")]] <- celltype_sub[nn[, i]]
    neighbor_table[[paste0("nbr", i, "_niche")]] <- dominant_niche[nn[, i]]
  }
  fwrite(neighbor_table, file.path(sample_out, paste0(sample_id, "_spot_neighbor_table.csv")))

  fwrite(neighbor_enrich(tumor_idx, celltype_sub, nn, n_permutations),
         file.path(sample_out, paste0(sample_id, "_Tumor_neighbor_subtype_celltype.csv")))
  fwrite(neighbor_enrich(tumor_idx, dominant_niche, nn, n_permutations),
         file.path(sample_out, paste0(sample_id, "_Tumor_neighbor_CNV_niche.csv")))

  gene_available <- intersect(gene_order, rownames(expr))
  neighbor_spots <- unique(unlist(neighbor_table[tumor_idx, paste0("nbr", seq_len(k_neighbors), "_spot"), with = FALSE], use.names = FALSE))
  neighbor_spots <- intersect(neighbor_spots, colnames(expr))

  expr_neighbor_matrix <- expr[gene_available, neighbor_spots, drop = FALSE]
  expr_neighbors <- as.data.table(t(expr_neighbor_matrix))
  expr_neighbors[, spot := neighbor_spots]
  cell2loc_neighbors <- as.data.table(cell2loc[neighbor_spots, , drop = FALSE])
  cell2loc_neighbors[, spot := neighbor_spots]
  spot_nb <- merge(cell2loc_neighbors, expr_neighbors, by = "spot", sort = FALSE)

  subtype_cols <- intersect(names(major_map), colnames(spot_nb))
  for (major_name in major_order) {
    subtype_for_major <- subtype_cols[unname(major_map[subtype_cols]) == major_name]
    if (length(subtype_for_major) > 0) {
      spot_nb[, (major_name) := rowSums(.SD), .SDcols = subtype_for_major]
    }
  }

  major_present <- major_order[major_order %in% names(spot_nb)]
  major_matrix <- spot_nb[, c("spot", major_present, gene_available), with = FALSE]
  fwrite(major_matrix, file.path(sample_out, paste0(sample_id, "_NeighborSpot_majorCell_Iron_matrix.csv")))

  pb_major <- rbindlist(lapply(major_present, function(x) {
    weighted_neighbor_expression(spot_nb, x, gene_available)
  }), fill = TRUE)
  fwrite(pb_major, file.path(sample_out, paste0(sample_id, "_PB_weighted_neighbor_byMajorCell.csv")))

  data.table(
    Sample = sample_id,
    Total_spots = length(common_spots),
    Tumor_dominant_spots = length(tumor_idx),
    Neighbor_spots_used = length(neighbor_spots),
    Genes_available = length(gene_available),
    Major_cell_types = length(major_present)
  )
}

sample_summary <- rbindlist(lapply(samples, run_one_sample), fill = TRUE)
fwrite(sample_summary, file.path(output_dir, "Figure_6F_full_neighbor_run_summary.csv"))

long_dt <- rbindlist(lapply(samples, function(sample_id) {
  file <- file.path(output_dir, paste0(sample_id, "_PB_weighted_neighbor_byMajorCell.csv"))
  x <- fread(file)
  x[, Sample := sample_id]
  melt(x, id.vars = c("Weight", "Sample"), variable.name = "Gene", value.name = "Expr")
}), fill = TRUE)
fwrite(long_dt, file.path(output_dir, "Figure_6F_full_Combined_long.csv"))

matrix_dt <- dcast(long_dt, Gene ~ Sample + Weight, value.var = "Expr")
fwrite(matrix_dt, file.path(output_dir, "Figure_6F_full_Combined_matrix.csv"))

figure_gene_sets <- list(
  iron_uptake = c("Trf", "Tfrc", "Slc11a2", "Slc39a14", "Slc39a8", "Slc7a11", "Slc25a37", "Slc25a28"),
  iron_storage = c("Fth1", "Hfe", "Steap3", "Ncoa4", "Ireb2", "Cp"),
  heme_metabolism = c("Alas1", "Uros", "Ppox", "Fech", "Hmox1", "Hmox2", "Blvra", "Blvrb"),
  iron_trafficking = c("Lcn2", "Slc22a17", "Ltf", "Slc46a1", "Slc48a1", "Abcb10", "Flvcr1", "Hp", "Lrp1"),
  iron_TME = c("Col1a1", "Col1a2", "Fn1", "Acta2", "Mmp3", "Vegfa", "Stat3", "Bach1", "Kdm6b", "Il1b", "Lif", "C3", "Selenbp1", "Cd44")
)

figure_manifest <- rbindlist(lapply(names(figure_gene_sets), function(gs) {
  data.table(GeneSet = gs, Gene = figure_gene_sets[[gs]], GeneOrder = seq_along(figure_gene_sets[[gs]]))
}))
fwrite(figure_manifest, file.path(output_dir, "Figure_6F_full_gene_set_manifest.csv"))

figure_genes <- figure_manifest$Gene
figure_source <- matrix_dt[match(figure_genes, Gene)]
figure_source[, GeneSet := figure_manifest$GeneSet[match(Gene, figure_manifest$Gene)]]
setcolorder(figure_source, c("Gene", "GeneSet", setdiff(names(figure_source), c("Gene", "GeneSet"))))
fwrite(figure_source, file.path(output_dir, "Figure_6F_full_figure_source_data.csv"))

cat("Completed full neighbor-to-source-data analysis.\n")
cat("Input:", h5ad_dir, "\n")
cat("Output:", normalizePath(output_dir, mustWork = FALSE), "\n")
cat("Samples:", paste(samples, collapse = ", "), "\n")
cat("Neighbor definition: six nearest spots from H5AD spatial coordinates.\n")
cat("Figure 6F focus genes:", nrow(figure_source), "\n")
