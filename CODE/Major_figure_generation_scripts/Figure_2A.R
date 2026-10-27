############################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(reshape2)
  library(tibble)
  library(tidyr)
  library(pROC)
  library(GSVA)
  library(pbapply)
  library(readxl)
  library(cogena)
  library(gelnet)
  library(matrixStats)
  library(pcaPP)      
  library(edgeR)      
  library(reticulate) 
  library(SummarizedExperiment) 
})

############################################################
## 0. Path & log
############################################################

##. dataset path
dataset_path <- "./TISI/Step3.1_in_vitro_validation/Benchmark_datasets/"

##. model path
tisi10_path     <- "./TISI/Step3.1_in_vitro_validation/Benchmark_models/SenOCLR_L2_10_final_dot.rdata"

##. signature path
data_dir <- file.path("./TISI/Step3.1_in_vitro_validation/Data")
extra_sig_file <- "./TISI/Step3.1_in_vitro_validation/Data/Sene_signatures.csv"

##. save path
outdir <- "./TISI/Step3.1_in_vitro_validation/Output/Figure_2A"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

##. log file
log_file <- file.path(outdir, "Validation_log.txt")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = " "), "\n")
  cat(msg)
  write(msg, file = log_file, append = TRUE)
}
log_msg("===== Validation of TISI START =====")


############################################################
## 1. Define functions
############################################################

GetAssayData2 <- function(object) {
  if ("Seurat" %in% class(object)) {
    return(Seurat::GetAssayData(object, layer = "counts"))
  } else {
    return(as.matrix(object))
  }
}

normalize_exp_no_log <- function(x){
  x_scaled <- scale(x)
  row_means <- rowMeans(x_scaled, na.rm = TRUE)
  x_centered <- x_scaled - matrix(
    row_means,
    nrow = nrow(x_scaled),
    ncol = ncol(x_scaled),
    byrow = FALSE
  )
  return(x_centered)
}

run_ssgsea <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  if (is.null(rownames(expr_mat))) {
    stop("ssGSEA input must have gene names as row names.")
  }

  gene_sets <- lapply(gene_sets, function(gs) {
    unique(as.character(gs[gs %in% rownames(expr_mat)]))
  })
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  if (length(gene_sets) == 0) {
    stop("No genes from the requested signature are present in the expression matrix.")
  }

  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(exprs = expr_mat)
    )
    param <- GSVA::ssgseaParam(
      exprData = se,
      geneSets = gene_sets,
      normalize = TRUE
    )
    res <- GSVA::gsva(param = param, verbose = FALSE)
    return(SummarizedExperiment::assay(res, "es"))
  }

  ## Legacy GSVA (<2.0) fallback.
  GSVA::gsva(
    expr_mat,
    gene_sets,
    method = "ssgsea",
    ssgsea.norm = TRUE,
    verbose = FALSE
  )
}

scoreOCLR <- function(profile, model, m = 'spearman') {
  w <- model$w
  if (is.null(names(w)) || anyNA(names(w))) {
    stop("OCLR model weights must be a named vector.")
  }
  genes <- intersect(names(w), rownames(profile))
  if (length(genes) == 0) {
    stop("No OCLR model genes are present in the expression matrix.")
  }
  if (length(genes) < length(w)) {
    warning(sprintf(
      "Only %d/%d OCLR genes are present; the score uses a partial model.",
      length(genes), length(w)
    ))
  }
  w <- w[genes]
  profile_mat <- as.matrix(profile[genes, , drop = FALSE])
  if (m == 'spearman') {
    scores <- list(apply(profile_mat, 2, function(s) cor(w, s, method = "spearman")))
  } else if (m == 'pearson') {
    scores <- list(apply(profile_mat, 2, function(s) cor(w, s, method = "pearson")))
  } else if (m == 'dot') {
    scores <- list(apply(profile_mat, 2, function(x) (w %*% x)))
  } else if (m == 'logistic') {
    scores <- list(apply(profile_mat, 2, function(x) plogis(sum(w * x))))
  } else if (m == "gsva") {

    n_values <- unique(pmin(c(100, 200, 500, 1000), length(w)))
    n_values <- n_values[n_values > 0]
    TISI <- pbapply::pblapply(n_values, function(n) {
      top_w <- sort(w, decreasing = TRUE)[seq_len(n)]
      run_ssgsea(profile_mat, list(signature = names(top_w)))[1, ]
    }, cl = 1)
    names(TISI) <- as.character(n_values)
    scores <- TISI
  } else if (m == 'mean') {
    n_max <- min(500, length(w))
    pb <- txtProgressBar(min = 0, max = n_max, style = 3)
    scores <- list()
    for (n in seq_len(n_max)) {
      setTxtProgressBar(pb, n)
      top_w <- sort(w, decreasing = TRUE)[seq_len(n)]
      if (n == 1) {
        score <- profile_mat[names(top_w), ]
      } else {
        score <- colMeans(profile_mat[names(top_w), , drop = FALSE])
      }
      scores[[as.character(n)]] <- score
    }
    close(pb)
  } else {
    stop("Unsupported OCLR scoring method: ", m)
  }
  return(scores)
}

minmax <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(!is.finite(x))) {
    return(rep(NA_real_, length(x)))
  }
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  x_scale <- (x - rng[1]) / (rng[2] - rng[1])
  return(x_scale)
}

calc_laf <- function(mat, genes) {
  genes <- genes %>% dplyr::filter(Gene %in% rownames(mat))
  if (nrow(genes) == 0) return(rep(NA_real_, ncol(mat)))
  x <- mat[genes$Gene, , drop = FALSE]
  x <- t(t(x) / matrixStats::colMedians(mat)) > 1
  ## V3 FIX: invoke_map_lgl() is deprecated in purrr >= 1.0.
  res <- vapply(seq_len(ncol(x)), function(j) {
    ifelse(
      genes$Direction == "Increase",
      x[, j],
      !x[, j]
    )
  }, logical(nrow(x)))
  dimnames(res) <- list(genes$Gene, colnames(mat))
  colSums(res) / nrow(res)
}

loadSetdata <- function(data_dir) {
  EnrichSet <- cogena::gmt2list(file.path(data_dir, "gene_50signatures_merge.gmt"))
  set_map <- c(
    SenMayo   = "SenMayo",
    CellAge   = "CellAge",
    GenAge    = "GenAge",
    ASIG      = "ASIG_all",
    SASP      = "SASP_Pathway",
    AgingAtlas = "AgingAtlas",
    SenUp     = "over_e"
  )
  missing_sets <- setdiff(unname(set_map), names(EnrichSet))
  if (length(missing_sets) > 0) {
    stop("Missing expected GMT signatures: ", paste(missing_sets, collapse = ", "))
  }
  EnrichSet <- EnrichSet[unname(set_map)]
  names(EnrichSet) <- names(set_map)
  rep_sene_genes <- readxl::read_excel(
    file.path(data_dir, "SigRS.xls"),
    sheet = 5, skip = 2
  )$Gene_names
  if (is.null(rep_sene_genes)) {
    stop("SigRS.xls sheet 5 must contain a Gene_names column.")
  }
  EnrichSet$SigRS <- unique(as.character(rep_sene_genes[!is.na(rep_sene_genes)]))
  return(EnrichSet)
}

loadlafdata <- function(data_dir) {
  laf_DAS <- fread(file.path(data_dir, "LaffertyWhyte2010_DAS.csv"))
  laf_mSS <- fread(file.path(data_dir, "LaffertyWhyte2010_mSS.csv"))
  laf_DAS_mSS <- rbind(laf_DAS, laf_mSS)
  laf <- list(DAS = laf_DAS, mSS = laf_mSS, DAS_mSS = laf_DAS_mSS)
  return(laf)
}

loadpredictorFix <- function() {
  CSS <- c(0.2921, 0.1810, 0.0524, 0.1285)
  names(CSS) <- c('C1QTNF6', 'SQOR', 'LYPD3', 'FAM83B')
  lassoCS <- c(-0.158, -0.153, 0.17, 0.347, -0.215, 0.288, -0.196,
               -0.245, -0.103, 0.127)
  names(lassoCS) <- c("ITGA8", "SEMA3G", "DPYSL3", "IFITM1", "ZNF521",
                      "SOCS3", "PCSK6", "DUSP1", "SLC44A4", "IL20RB")
  ECcores36 <- c(0.756, 0.739, 0.629, 0.472, 0.439, 0.429, 0.371,
                 0.332, 0.322, 0.3, 0.281, 0.263, 0.251, 0.235, 0.189,
                 0.186, 0.176, 0.14, 0.134, 0.128, 0.128, 0.102, 0.095,
                 0.067, 0.046, 0.042, 0.04, 0.013, -0.022, -0.078,
                 -0.162, -0.186, -0.192, -0.215, -0.464, -0.593)
  names(ECcores36) <- c('IGFBP5', 'IFI27', 'PLAT', 'MX1', 'IFIT1',
                        'SEMA3A', 'TAGLN', 'LAMC2', 'C19orf33', 'OAS2',
                        'OAS1', 'LHX6', 'CD44', 'ABI3BP', 'KCNMA1',
                        'P3H2', 'IFITM1', 'DIRAS3', 'PLA2G4C', 'IL32',
                        'ADIRF', 'ISG15', 'ADGRF5', 'SAMD9L', 'LYVE1',
                        'IFI6', 'SERPINE2', 'DDX58', 'NUF2', 'CD200',
                        'BMX', 'TOP2A', 'FST', 'PEG10', 'COLEC12',
                        'PHGDH')
  HCCcores19 <- c(-0.24208, -0.32203, -0.15081, 0.66330, 0.35072,
                  -0.46527, 0.64921, -0.07426, 0.37272, 0.16174,
                  0.63189, -0.61190, -0.33728, 0.37964, 0.45727,
                  0.06315, 0.20073, -0.16337, -0.48263)
  names(HCCcores19) <- c("CDCA5", "CENPF", "CENPW", "CDCA8", "SPC25",
                         "CDKN3", "CENPA", "BUB1", "DLGAP5", "IGSF3",
                         "HMMR", "TOP2A", "RAD54L", "TTK", "GINS1",
                         "PTTG1", "ETV4", "GINS2", "PKMYT1")
  coefs <- list(
    CSS = CSS,
    lassoCS = lassoCS,
    ECcores36 = ECcores36,
    HCCcores19 = HCCcores19
  )
  return(coefs)
}

cal_SENCAN <- function(object) {
  ## NOTE (V5-SENCAN): SENCAN must receive a non-negative, count-like
  ## expression matrix. Do not pass the z-scored/row-centered matrix used by
  ## TISI and the other score calculations.
  dm <- GetAssayData2(object) %>% data.matrix()
  if (is.null(rownames(dm)) || is.null(colnames(dm))) {
    stop("SENCAN input must have gene and sample names.")
  }
  if (any(!is.finite(dm)) || any(dm < 0)) {
    stop("SENCAN requires finite, non-negative, count-like expression values.")
  }
  required_objects <- c(
    "dge", "dmCPM", "gene_mean", "gene_sd", "full_classifier",
    "use_gene_ix", "edger_ref_column", "calcNormFactors_WSPref"
  )
  missing_objects <- required_objects[
    !vapply(required_objects, exists, logical(1), inherits = TRUE)
  ]
  if (length(missing_objects) > 0) {
    stop("SENCAN classifier is incomplete; missing: ",
         paste(missing_objects, collapse = ", "))
  }
  genes <- rownames(dm)
  gene_ref  <- dge$genes$gene_name
  count_ref <- dge$counts
  map <- match(gene_ref, genes)
  dmMapped  <- dm[map, ]
  dmImputed <- dmMapped
  for (i in 1:ncol(dmImputed)) {
    sample <- dmImputed[, i]
    if (any(is.na(sample))) {
      cors <- rep(NA, ncol(count_ref))
      names(cors) <- rownames(dge$samples)
      for (j in 1:ncol(count_ref)) {
        have_val <- which(!is.na(sample))
        cors[j] <- pcaPP::cor.fk(sample[have_val], count_ref[have_val, j])
      }
      nearest  <- head(order(cors, decreasing = TRUE, na.last = NA), n = 2)
      if (length(nearest) == 0) {
        stop("SENCAN imputation failed: no valid reference correlations.")
      }
      need_val <- which(is.na(sample))
      dmImputed[need_val, i] <- apply(count_ref[need_val, nearest], 1, mean)
    }
  }

  dge_pred <- edgeR::DGEList(counts = dmImputed,
                             genes = data.frame(ensembl_id = dge$genes$ensembl_id))
  dge_pred <- calcNormFactors_WSPref(dge_pred, edger_ref_column)
  dmPredCPM <- edgeR::cpm(dge_pred, log = TRUE)
  dmPredCPM[dmPredCPM < min(dmCPM)] <- min(dmCPM)
  pred_standardized_expression <- (dmPredCPM - gene_mean) / gene_sd

  ## NOTE (V5-SENCAN): the classifier file stores the fitted glmnet model in
  ## full_classifier$glmnet.fit and stores the selected lambda separately.
  ## Predict on that fitted object rather than on the wrapper list itself.
  if (is.list(full_classifier) && !is.null(full_classifier$glmnet.fit)) {
    classifier_fit <- full_classifier$glmnet.fit
    classifier_s <- full_classifier$lambda.1se
    if (is.null(classifier_s)) {
      stop("SENCAN classifier does not contain lambda.1se.")
    }
  } else if (inherits(full_classifier, "cv.glmnet")) {
    classifier_fit <- full_classifier
    classifier_s <- "lambda.1se"
  } else {
    stop("Unsupported SENCAN classifier format; expected glmnet.fit/lambda.1se.")
  }

  res <- predict(
    classifier_fit,
    newx = t(pred_standardized_expression[use_gene_ix, , drop = FALSE]),
    s = classifier_s
  )[, 1]
  sencan_score <- 1 / (1 + exp(-res))
  return(sencan_score)
}

cal_SENCID <- function(object, denoise = FALSE) {
  tryCatch({
    SenCIDapi <- reticulate::import("SenCID.api")
    sc       <- reticulate::import("scanpy")
    np       <- reticulate::import("numpy")
    
    counts <- GetAssayData2(object) %>% as.matrix()
    adata  <- sc$AnnData(t(counts))
    adata$var_names <- rownames(counts)
    adata$obs_names <- colnames(counts)
    res <- SenCIDapi$SenCID(
      adata    = adata,
      sidnums  = as.integer(c(1, 2, 3, 4, 5, 6)),
      denoising = denoise,
      binarize  = TRUE,
      threads   = as.integer(1),
      savetmp   = TRUE
    )
    recSID <- res[[2]]
    SID <- gsub('rec_', '', which.max(colMeans(recSID[, -1])) %>% names())
    print(sprintf('%s was used in this dataset', SID))
    pred_dict <- res[[1]][[SID]]
    s <- pred_dict$SID_Score
    return(s)
  },
  error = function(e) {
    message("SENCID failed: ", conditionMessage(e))
    return(NULL)
  })
}

calc_scores <- function(methodslist, object = NULL, raw_counts = NULL, type = NULL) {
  mat <- GetAssayData2(object)
  data <- data.frame(cellname = colnames(mat), stringsAsFactors = FALSE)
  if (!is.null(type)) {
    methodslist <- methodslist[type]
  }
  for (m in names(methodslist)) {
    message(sprintf("Start calculating senescence scores by %s ...", m))
    if (m == 'SenMarker') {
      scorelist <- intersect(methodslist[[m]], rownames(mat))
      message(sprintf("%d senescence marker was found in this dataset",
                      length(scorelist)))
      if (length(scorelist) > 0) {
        score <- t(mat[scorelist, , drop = FALSE]) %>% data.frame()
        rownames(score) <- colnames(mat)
        data <- cbind(data, score)
      }
    }
    
    if (m == 'SenSet') {
      scorelist <- methodslist[[m]]
      for (sl in scorelist) {
        message(sl)
        if (!sl %in% names(SenSetlist)) {
          warning("Signature ", sl, " not found in SenSetlist; skip.")
          next
        }
        score <- tryCatch(
          run_ssgsea(as.matrix(mat), list(signature = SenSetlist[[sl]]))[1, ],
          error = function(e) {
            warning("ssGSEA failed for ", sl, ": ", conditionMessage(e))
            rep(NA_real_, ncol(mat))
          }
        )
        data[[sl]] <- score
      }
    }
    
    if (m == "MachineScore") {
      for (sl in methodslist[[m]]) {
        message(sl)
        
        ## ---------------------- traditional ML --------------------------
        if (sl %in% c("CSS", "lassoCS", "ECcores36", "HCCcores19")) {
          all_genes <- names(coefslist[[sl]])
          genes <- intersect(all_genes, rownames(mat))
          
          if (length(genes) == 0) {
            message(sprintf("⚠ %s: no genes found!", sl))
            next
          }
          
          if (length(genes) < length(all_genes)) {
            warning(sprintf(
              "%s uses only %d/%d coefficient genes in this dataset.",
              sl, length(genes), length(all_genes)
            ))
          }
          tmp <- sweep(
            as.matrix(mat[genes, , drop = FALSE]),
            MARGIN = 1,
            STATS = coefslist[[sl]][genes],
            FUN = "*"
          )
          score <- colSums(tmp)
          
          data[[sl]] <- score
          next
        }
        
        ## ---------------------- SENCAN -------------------------------
        if (sl == "SENCAN") {
          ## NOTE (V5-SENCAN): object contains normalized values for the
          ## other methods. SENCAN must use the separately supplied raw matrix.
          if (is.null(raw_counts)) {
            stop("SENCAN requires raw_counts; no valid SENCAN input was supplied.")
          }
          data[[sl]] <- cal_SENCAN(raw_counts)
          next
        }
        
        ## ---------------------- SENCID -------------------------------
        if (sl == "SENCID") {
          score <- cal_SENCID(object)
          if (!is.null(score)) {
            data[[sl]] <- score
          }
          next
        }
      }
      next
    }
    
    if (m == 'TraditionScore') {
      scorelist <- methodslist[[m]]
      for (sl in scorelist) {
        message(sl)
        if (sl %in% c('DAS', 'mSS', 'DAS_mSS')) {
          score <- calc_laf(as.matrix(mat), laflist[[sl]])
          data[[sl]] <- score
        }
        if (sl == 'CS_score') {
          cs_sets <- CSsccorelist[c('Underexpressed', 'Overexpressed')]
          cs_sets <- lapply(cs_sets, function(gs) {
            unique(as.character(gs[gs %in% rownames(mat)]))
          })
          if (any(lengths(cs_sets) == 0)) {
            warning("CS_score skipped: one or both CS gene sets have no genes in this dataset.")
            data[[sl]] <- rep(NA_real_, ncol(mat))
            next
          }
          scorecs <- tryCatch(
            run_ssgsea(
              as.matrix(mat),
              cs_sets
            ),
            error = function(e) {
              warning("ssGSEA failed for CS_score: ", conditionMessage(e))
              matrix(NA_real_, nrow = 2, ncol = ncol(mat),
                     dimnames = list(c("Underexpressed", "Overexpressed"),
                                     colnames(mat)))
            }
          )
          data[[sl]] <- scorecs[2, ] - scorecs[1, ]
        }
      }
    }
  }
  
  return(data)
}

############################################################
## 2. load objectif：SenSetlist / laflist / CSsccorelist / coefslist / SENCAN classifier
############################################################

log_msg("Load SenSet / LAF / CS score / coefs / SENCAN classifier ...")

laflist    <- loadlafdata(data_dir)
SenSetlist <- loadSetdata(data_dir)

## SENCAN classifier
sencan_file <- file.path(data_dir, "SENCAN_classifier.rda")
if (file.exists(sencan_file)) {
  log_msg("Load SENCAN classifier from: ", sencan_file)
  load(sencan_file)
} else {
  log_msg("WARNING: SENCAN_classifier.rda not found, SENCAN will fail.")
}

CSsccore <- read.csv(file.path(data_dir, "CS score_PMID31560156.csv"),
                     header = TRUE)
CSsccorelist <- list(
  Underexpressed = dplyr::filter(CSsccore, expression == "Underexpressed")$Gene.symbol,
  Overexpressed = dplyr::filter(CSsccore, expression == "Overexpressed")$Gene.symbol
)

coefslist <- loadpredictorFix()

## 8 extra SenSet

log_msg("Load extra SenSet from: ", extra_sig_file)
extra_sig_df <- fread(extra_sig_file)
if (ncol(extra_sig_df) < 2) {
  stop("extra signature file has <2 columns, please check format.")
}
extra_sig_df <- extra_sig_df[, 1:2, with = FALSE]
setnames(extra_sig_df, c("Signature", "Gene"))
extra_sets <- split(extra_sig_df$Gene, extra_sig_df$Signature)
SenSetlist <- c(SenSetlist, extra_sets)
extra_sig_names <- names(extra_sets)
log_msg("Extra signatures added: ", paste(extra_sig_names, collapse = ", "))


############################################################
## 3. load TISI-OCLR model
############################################################

log_msg("Load TISI (L2=10, dot) from: ", tisi10_path)
load(tisi10_path)     
tisi10_model <- SenOCLR_final
rm(SenOCLR_final)

############################################################
## 4. Methods list
############################################################

Methods <- list(
  SenMarker = c("GLB1", "TP53", "CDKN1A", "CDKN2A", "LMNB1", "IL1A",
                "RB1", "CDK1", "CDK4", "CDK6", "MKI67", "CDKN2B", "SERPINE1"),
  SenSet = c("SenMayo", "CellAge", "GenAge", "ASIG", "SASP",
             "AgingAtlas", "SenUp", "SigRS", extra_sig_names),
  MachineScore = c("CSS", "lassoCS", "ECcores36", "HCCcores19",
                   "SENCAN", "SENCID"),
  TraditionScore = c("DAS", "mSS", "DAS_mSS", "CS_score")
)

############################################################
## 5. Dataset
############################################################

datasets <- list(
  GSE108895 = list(expr = "GSE108895_norm_log2.csv", meta = "GSE108895_meta.csv"),
  GSE210334 = list(expr="GSE210334_CPM.csv",      meta="GSE210334_meta.csv"),
  GSE221180 = list(expr="GSE221180_nodup.csv",    meta="GSE221180_meta.csv"),  
  GSE221419 = list(expr="GSE221419_CPM.csv",      meta="GSE221419_meta.csv"),
  GSE50474  = list(expr = "GSE50474_norm_log2.csv",  meta = "GSE50474_meta.csv"),
  GSE94805  = list(expr = "GSE94805_nodup.csv",      meta = "GSE94805_meta.csv"),
  PC9_Gem_Eto = list(expr="PC9_Gem_Eto_nodup.csv", meta="PC9_Gem_Eto_meta.csv"),
  U251_TMZ = list(expr="U251_TMZ_nodup.csv",       meta="U251_TMZ_meta.csv")
)

############################################################
## 6. AUC calculation
############################################################

compute_auc_vec <- function(dat, label_factor, methods_vec) {
  sapply(methods_vec, function(meth) {
    if (!meth %in% colnames(dat)) return(NA_real_)
    vals <- dat[[meth]]
    keep <- is.finite(vals) & !is.na(label_factor)
    if (sum(keep) < 2 || length(unique(label_factor[keep])) < 2) {
      return(NA_real_)
    }
    roc_obj <- try(
      pROC::roc(label_factor[keep], vals[keep],
                levels = c("non-senescent", "senescent"),
                direction = "<"),
      silent = TRUE
    )
    if (inherits(roc_obj, "try-error")) return(NA_real_)
    as.numeric(pROC::auc(roc_obj))
  })
}

############################################################
## 7. Main cycle for every datset
############################################################

AUClist <- list()

for (ds in names(datasets)) {
  log_msg("---- Start dataset: ", ds, " ----")
  
  expr_file <- file.path(dataset_path, datasets[[ds]]$expr)
  meta_file <- file.path(dataset_path, datasets[[ds]]$meta)
  
  expr_df <- fread(expr_file)
  genes <- expr_df[[1]]
  if (anyNA(genes) || anyDuplicated(genes)) {
    stop("Dataset ", ds, " contains missing or duplicated gene identifiers.")
  }
  expr_mat <- as.matrix(expr_df[, -1, with = FALSE])
  rownames(expr_mat) <- genes
  storage.mode(expr_mat) <- "numeric"
  if (any(!is.finite(expr_mat))) {
    stop("Dataset ", ds, " contains non-finite expression values.")
  }
  
  meta_df <- fread(meta_file)

  if ("Sample" %in% colnames(meta_df)) {
    rownames(meta_df) <- meta_df$Sample
  } else if ("Gene_name" %in% colnames(meta_df)) {
    rownames(meta_df) <- meta_df$Gene_name
  } else {
    stop("Metadata has no recognizable sample ID column.")
  }
  
  meta_df <- as.data.frame(meta_df)
  
  common <- intersect(colnames(expr_mat), rownames(meta_df))
  if (length(common) < 2) {
    log_msg("Dataset ", ds, ": less than 2 overlapping samples, skip.")
    next
  }
  expr_mat <- expr_mat[, common, drop = FALSE]
  meta_df  <- meta_df[match(common, rownames(meta_df)), , drop = FALSE]
  
  ## Only keep Sene / non-Sene samples for "Group"
  if (!"Group" %in% colnames(meta_df)) {
    stop("Meta for ", ds, " has no 'Group' column!")
  }
  
  valid_samples <- rownames(meta_df)[meta_df$Group %in% c("Sene", "non-Sene")]
  
  samples_use <- intersect(valid_samples, colnames(expr_mat))
  
  expr_use <- expr_mat[, samples_use, drop = FALSE]
  meta_use <- meta_df[samples_use, , drop = FALSE]

  ## V3 FIX: avoid creating an empty/one-class object and an invalid ROC.
  if (length(samples_use) < 2 || length(unique(meta_use$Group)) < 2) {
    log_msg(
      "Dataset ", ds,
      ": fewer than 2 usable samples or fewer than 2 groups; skip."
    )
    next
  }
  if (!identical(colnames(expr_use), rownames(meta_use))) {
    stop("Expression and metadata sample order is not aligned for ", ds, ".")
  }
  
  ## NOTE (V5-SENCAN): keep the pre-normalization matrix separately for
  ## SENCAN. The normalized matrix is still used by the existing TISI and
  ## other score calculations below.
  expr_raw <- expr_use
  expr_use <- normalize_exp_no_log(expr_raw)
  
  label <- ifelse(meta_use$Group == "Sene", "senescent", "non-senescent")
  label_factor <- factor(label, levels = c("non-senescent", "senescent"))
  
  log_msg(ds, ": #samples used = ", length(samples_use),
          " (Sene=", sum(label == "senescent"),
          ", non-Sene=", sum(label == "non-senescent"), ")")
  
  obj <- CreateSeuratObject(
    counts    = expr_use,
    meta.data = data.frame(meta_use, check.names = FALSE)
  )
  
  dat <- calc_scores(Methods, object = obj, raw_counts = expr_raw)
  dat$label <- label
  rownames(dat) <- dat$cellname
  
  ## --- TISI（L2=10，14 dataset）dot ---
  profile_mat <- GetAssayData2(obj) %>% as.matrix()
  tisi10 <- scoreOCLR(profile = profile_mat, model = tisi10_model, m = "dot")[[1]]
  dat$TISI <- minmax(tisi10)
  
  ## save the scores
  score_outfile <- file.path(outdir, paste0(ds, "_AllScores.csv"))
  write.csv(dat, score_outfile, row.names = FALSE)
  log_msg(ds, ": scores saved to ", score_outfile)
  
  core_methods <- c( "TISI", unlist(Methods[c("SenSet", "MachineScore", "TraditionScore")]))
  auc1 <- compute_auc_vec(dat, label_factor, core_methods)
  
  markers <- Methods$SenMarker[Methods$SenMarker %in% colnames(dat)]
  auc2 <- compute_auc_vec(dat, label_factor, markers)
  
  auc_all <- c(auc1, auc2)
  AUClist[[ds]] <- auc_all
  
  auc_df <- data.frame(
    Method = names(auc_all),
    AUC    = as.numeric(auc_all),
    stringsAsFactors = FALSE
  )
  
  auc_outfile <- file.path(outdir, paste0(ds, "_AUC.csv"))
  write.csv(auc_df, auc_outfile, row.names = FALSE)
  log_msg(ds, ": AUC table saved to ", auc_outfile)
  log_msg("---- Finish dataset: ", ds, " ----")
}


############################################################
## 8. AUC summary and bar graph
############################################################

log_msg("Summarizing AUC across datasets ...")

if (length(AUClist) == 0) {
  stop("No dataset completed successfully; AUC summary cannot be generated.")
}

auc_all_df <- do.call(rbind, lapply(names(AUClist), function(ds) {
  data.frame(
    Dataset = ds,
    Method  = names(AUClist[[ds]]),
    AUC     = as.numeric(AUClist[[ds]]),
    stringsAsFactors = FALSE
  )
}))

auc_sum_file <- file.path(outdir, "AllDatasets_AUC_summary.csv")
write.csv(auc_all_df, auc_sum_file, row.names = FALSE)
log_msg("All datasets AUC summary saved to ", auc_sum_file)

auc_all_df2 <- auc_all_df[!is.na(auc_all_df$AUC), , drop = FALSE]

pdf(file.path(outdir, "AllDatasets_AUC_barplot.pdf"), width = 8, height = 5)
print(
  ggplot(auc_all_df2,
         aes(x = reorder(Method, -AUC, FUN = mean, na.rm = TRUE),
             y = AUC, fill = Dataset)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_bw() +
    labs(title = "AUC per method",
         x = "Method", y = "AUC")
)
dev.off()
log_msg("AllDatasets_AUC_barplot.pdf saved.")
