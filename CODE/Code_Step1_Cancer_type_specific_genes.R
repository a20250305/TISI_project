library(dplyr)
library(DESeq2)
library(Seurat)
library(mclust)     
library(ggplot2)

############################################################
## 0. Define paths
############################################################

counts_file <- "./Counts_Global.csv"
raw_rdata   <- "./raw_full.rdata"

outdir <- "./Machine_Learning/Step1_Cancer_type_specific_genes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "marker_sets"), showWarnings = FALSE)


############################################################
## 1. Load TPM expression & metadata
############################################################
load(raw_rdata)   # exprData, meta, gene_symbol

############################################################
## 2. Load raw counts
############################################################
counts <- read.csv(counts_file, check.names = FALSE)
cat("Counts:", nrow(counts), "genes ×", ncol(counts)-2, "samples\n\n")

## fix gene id
gene_id <- sub("\\..*", "", counts$Gene_ID)

## match counts rows to exprData rows
idx <- match(rownames(exprData), gene_id)

counts_mat <- counts[idx, -(1:2)]
rownames(counts_mat) <- gene_id[idx]

## match sample order
counts_mat <- counts_mat[, colnames(exprData)]

cat("Counts matched:", nrow(counts_mat), "×", ncol(counts_mat), "\n\n")


############################################################
## 3. Select non-senescent samples 
############################################################
meta_NS   <- meta %>% filter(label == "non-senescent")
counts_NS <- counts_mat[, rownames(meta_NS)]


############################################################
## 4: Cancer_type-specific DESeq2 + save marker sets
############################################################

# create marker output dir
marker_dir <- file.path(outdir, "marker_sets_Cancer_type")
dir.create(marker_dir, showWarnings = FALSE)

meta_NS <- meta %>% filter(label == "non-senescent")
counts_NS <- counts_mat[, rownames(meta_NS)]
stopifnot(all(colnames(counts_NS) == rownames(meta_NS)))

# Cancer_type with >=2 samples
tab_ct <- table(meta_NS$Cancer_type)
valid_ct <- names(tab_ct)[tab_ct >= 2]

DEGs_CT <- list()

for (ct in valid_ct) {
  cat(">>> Cancer_type:", ct, "\n")

  m <- meta_NS
  m$condition <- ifelse(m$Cancer_type == ct, ct, "other")

  if (sum(m$condition == ct) < 2 || sum(m$condition == "other") < 2) {
    cat("  Skip (not enough samples)\n\n")
    next
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts_NS,
    colData   = m,
    design    = ~ condition
  )
  
  dds <- DESeq(dds, fitType = "mean")
  res <- as.data.frame(results(dds, contrast = c("condition", ct, "other")))
  res <- res[!is.na(res$padj), ]
  
  # save full DESeq2 table
  out_full <- file.path(marker_dir, paste0("DESeq2_full_", ct, ".csv"))
  write.csv(res, out_full, row.names = TRUE)
  
  # filter markers (log2FC>1 & padj<0.05)
  res2 <- res %>% filter(log2FoldChange > 1, padj < 0.05)
  res2 <- res2[order(res2$log2FoldChange, decreasing = TRUE), ]
  
  DEGs_CT[[ct]] <- res2
  
  # save filtered marker list
  out_filtered <- file.path(marker_dir, paste0("Markers_", ct, "_filtered.csv"))
  write.csv(res2, out_filtered, row.names = TRUE)
  
  cat("  Total DEGs:", nrow(res), " | Markers kept:", nrow(res2), "\n")
  cat("  Saved:", out_full, "\n")
  cat("  Saved:", out_filtered, "\n\n")
}

############################################################
## 5: ARI sweep for Cancer_type (Seurat v5 FIXED)
############################################################

deg_file <- file.path(outdir, "DEGs_CancerType_list.rds")
DEGs_CT  <- readRDS(deg_file)

idx_sen   <- which(meta$label == "senescent")
X_sen     <- exprData[, idx_sen]
sen_meta  <- meta[idx_sen, ]

gene_nums <- c(1:10, seq(15,100,5), 120,150,200,250,300)
ari_res <- data.frame()

for (gene_num in gene_nums) {

    # Select marker genes per cancer type
    MarkerSets <- lapply(DEGs_CT, function(x){
        if (nrow(x) == 0) return(character(0))
        rownames(x)[1:min(nrow(x), gene_num)]
    })

    rm_genes <- unique(unlist(MarkerSets))
    X_new <- X_sen[!(rownames(X_sen) %in% rm_genes), ]

    # -------------------------------------------------------
    # Create SeuratObject (log1p TPM)
    # -------------------------------------------------------
    mat_log <- log1p(X_new)

    sce <- CreateSeuratObject(
        counts = mat_log,
        meta.data = sen_meta,
        min.cells = 3,
        min.features = 200
    )

    # ************* Seurat v5 FIX ***************
    # Assay5 stores counts/data as $counts / $data
    sce[["RNA"]]$data <- sce[["RNA"]]$counts
    # ********************************************

    # -------------------------------------------------------
    # Standard Seurat pipeline
    # -------------------------------------------------------
    sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
    sce <- ScaleData(sce, features = rownames(sce))
    sce <- RunPCA(sce, features = VariableFeatures(sce))
    sce <- FindNeighbors(sce, dims = 1:20)
    sce <- FindClusters(sce, resolution = 0.6)
    
    # -------------------------------------------------------
    # Compute ARI
    # -------------------------------------------------------
    ari_val <- adjustedRandIndex(
        sen_meta$Cancer_type,
        sce$RNA_snn_res.0.6
    )

    ari_res <- rbind(ari_res,
                     data.frame(gene_num = gene_num,
                                ARI = ari_val))
}

# Save results
write.csv(ari_res, file.path(outdir, "ARI_results.csv"),
          row.names = FALSE)

# Find elbow point (max second derivative)
d2 <- diff(diff(ari_res$ARI))
elbow_idx <- which.min(d2) + 1
elbow_gene_num <- ari_res$gene_num[elbow_idx]

write(elbow_gene_num,file.path(outdir, "ARI_elbow_point.txt"))


############################################################
## 6 — Remove Cancer_type-specific markers 
############################################################

raw_rdata <- "./TISI_preprocess/raw_full.rdata"
load(raw_rdata)   # exprData, meta, gene_symbol

DEGs_CT <- readRDS(deg_file)

TopN <-50   

MarkerSets <- lapply(DEGs_CT, function(x){
    if (nrow(x) == 0) return(character(0))
    x <- x[order(x$log2FoldChange, decreasing = TRUE), ]
    return(rownames(x)[1:min(nrow(x), TopN)])
})

rm_genes <- unique(unlist(MarkerSets))

writeLines(rm_genes, file.path(outdir, paste0("removed_marker_genes_Top", TopN, ".txt")))

exprData_filtered <- exprData[!(rownames(exprData) %in% rm_genes), ]

save(exprData_filtered, meta, file = file.path(outdir, paste0("exprData_dropCancerType_Top", TopN, ".rdata")))

