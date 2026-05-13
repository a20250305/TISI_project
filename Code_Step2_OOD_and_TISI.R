suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gelnet)
  library(GSVA)
  library(pbapply)
  library(pROC)
})


#### path and functions

in_rdata_drop <- "./data/exprData_dropCancerType_Top50.rdata"

outdir <- "./OOD_and_TISI"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

normalize_exp <- function(x){
  x_log  <- log2(x + .00001) %>% scale()
  x_centre <- x_log - (apply(x_log, 1, mean))
  return(x_centre)
}


scoreOCLR <- function(profile, model, m = "spearman") {
  w <- model$w
  genes <- intersect(names(w), rownames(profile))
  w <- w[genes]

  profile <- profile[genes, , drop = FALSE]
  profile_mat <- as.matrix(profile)
  attr(profile_mat, "scaled:center") <- NULL
  attr(profile_mat, "scaled:scale")  <- NULL

  if (m == "spearman") {
    TISI <- list(apply(profile_mat, 2, function(s) cor(w, s, method = "spearman")))

  } else if (m == "pearson") {
    TISI <- list(apply(profile_mat, 2, function(s) cor(w, s, method = "pearson")))

  } else if (m == "dot") {
    TISI <- list(apply(profile_mat, 2, function(x) (w %*% x)))

  } else if (m == "logistic") {
    TISI <- list(apply(profile_mat, 2, function(x) plogis(sum(w * x))))

  } else if (m == "gsva") {

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(exprs = profile_mat)
    )

    n_values <- c(100, 200, 500, 1000)
    TISI <- pbapply::pblapply(n_values, function(n) {
      top_w <- sort(w, decreasing = TRUE)[1:n]

      param <- GSVA::ssgseaParam(
        exprData  = se,
        geneSets  = list(signature = names(top_w)),
        normalize = TRUE
      )

      res_se <- GSVA::gsva(
        param   = param,
        verbose = FALSE
      )

      es <- SummarizedExperiment::assay(res_se, "es")
      score <- es[1, ]
      return(score)
    }, cl = 1)

    names(TISI) <- as.character(n_values)

  } else if (m == "mean") {
    pb <- txtProgressBar(min = 0, max = 500, style = 3)
    TISI <- list()
    for (n in 1:500) {
      setTxtProgressBar(pb, n)
      top_w <- sort(w, decreasing = TRUE)[1:n]
      if (n == 1) {
        score <- profile_mat[names(top_w), ]
      } else {
        score <- colMeans(profile_mat[names(top_w), , drop = FALSE])
      }
      TISI[[as.character(n)]] <- score
    }
    close(pb)
  }

  return(TISI)
}

minmax <- function(x){
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

############################################################
## Step 1. loading exprData_filtered + meta
############################################################

objs <- load(in_rdata_drop)
exprData <- exprData_filtered
rm(exprData_filtered)
meta <- meta[colnames(exprData), , drop = FALSE]

############################################################
## Step 2. generate OOD splits by Cancer_type
############################################################

cancertypes <- unique(meta$Cancer_type)
set.seed(2024)
cancertypes_shuffled <- replicate(3, sample(cancertypes, size = length(cancertypes)))
chunks <- apply(cancertypes_shuffled, 2, function(x)
  split(x, cut(seq_along(x), 10, labels = FALSE)))

cancertype_train <- list()
for (t in 1:3) {
  for (k in 1:10) {
    s <- paste0("t_", t, ":k_", k)
    cancertype_train[[s]] <- cancertypes[!cancertypes %in% chunks[[t]][[k]]]
  }
}

train_test_data <- list()
for (i in seq_along(cancertype_train)) {
  nm <- names(cancertype_train)[i]
  idx <- which(meta$Cancer_type %in% cancertype_train[[i]])
  train_set <- exprData[, idx, drop = FALSE] %>% normalize_exp()
  test_set  <- exprData[, -idx, drop = FALSE] %>% normalize_exp()
  train_test_data[[nm]] <- list(train_set, test_set)
}


############################################################
## Step 3. train OCLR models
############################################################

lambda_grid <- c(0.001, 0.01, 0.1, 1, 5, 10)
SenOCLR_L2 <- list()

for (l2 in lambda_grid) {
  for (i in seq_along(train_test_data)) {
    nm <- names(train_test_data)[i]
    train_set <- train_test_data[[i]][[1]]
    sen_ids <- rownames(meta[meta$label == "senescent", ])
    train_set_sen <- train_set[, colnames(train_set) %in% sen_ids, drop = FALSE]
    if (ncol(train_set_sen) < 2) next
    SenOCLR <- gelnet(t(train_set_sen), NULL, 0, l2)
    n_model  <- paste0("l2_", as.character(l2), ":", nm)
    SenOCLR_L2[[n_model]] <- SenOCLR
  }
}


############################################################
## Step 4. score test set
############################################################

TISIs_list <- list()
for (n in names(SenOCLR_L2)) {

  parts <- strsplit(n, ":")[[1]]
  test_set <- train_test_data[[paste(parts[2], parts[3], sep=":")]][[2]]

  TISIs_list_m <- list()
  for (m in c("spearman","pearson","dot","logistic","gsva","mean")) {
    TISIs_list_m[[m]] <- scoreOCLR(test_set, SenOCLR_L2[[n]], m)
  }
  TISIs_list[[n]] <- TISIs_list_m
}


############################################################
## Step 5. AUC per Cancer_type
############################################################

auc_models <- list()
for (q in c("spearman","pearson","logistic","dot")) {

  TISIs_q <- lapply(TISIs_list, function(x) x[[q]][[1]])
  auc_l2 <- list()

  for (l2 in lambda_grid) {

    l2_name <- paste0("l2_", as.character(l2))
    nms <- names(TISIs_q)[grepl(paste0("^", l2_name, ":"), names(TISIs_q))]

    auc_celltypes <- list()
    for (nn in nms) {
      pred <- TISIs_q[[nn]]
      samples <- names(pred)
      test_label <- meta[samples, "label"]
      celltypes_sam <- meta[samples, "Cancer_type"]

      for (ct in unique(celltypes_sam)) {
        idx <- which(celltypes_sam == ct)
        pred_c <- pred[idx]
        test_label_c <- test_label[idx]
        if (length(unique(test_label_c)) < 2) next

        auc_val <- pROC::auc(
          pROC::roc(test_label_c, pred_c,
                    levels=c("non-senescent","senescent"), direction="<")
        )

        if (!ct %in% names(auc_celltypes))
          auc_celltypes[[ct]] <- as.numeric(auc_val)
        else
          auc_celltypes[[ct]] <- c(auc_celltypes[[ct]], as.numeric(auc_val))
      }
    }

    auc_l2[[l2_name]] <- auc_celltypes
  }

  auc_models[[q]] <- auc_l2
}


auc_df <- do.call(rbind, lapply(names(auc_models), function(q){
  do.call(rbind, lapply(names(auc_models[[q]]), function(l2){
    do.call(rbind, lapply(names(auc_models[[q]][[l2]]), function(ct){
      data.frame(method=q,
                 L2=as.numeric(gsub("l2_","",l2)),
                 Cancer_type=ct,
                 AUC=auc_models[[q]][[l2]][[ct]])
    }))
  }))
}))

write.csv(auc_df, file.path(outdir,"Step5_AUC_perCancerType_perL2_method.csv"), row.names=FALSE)

auc_summary <- auc_df %>%
  group_by(method,L2) %>%
  summarise(meanAUC=mean(AUC,na.rm=TRUE),
            sdAUC=sd(AUC,na.rm=TRUE), .groups="drop")

write.csv(auc_summary, file.path(outdir,"Step5_AUC_summary_by_L2_and_method.csv"), row.names=FALSE)

best_row <- auc_summary[which.max(auc_summary$meanAUC), ]
lambda_best <- best_row$L2
method_best <- as.character(best_row$method)

############################################################
## Step 6 – final model
############################################################

X_norm <- normalize_exp(exprData)
sen_ids <- rownames(meta[meta$label == "senescent", ])
idx_sen <- colnames(X_norm) %in% sen_ids
X_tr <- X_norm[, idx_sen, drop=FALSE]

SenOCLR_final <- gelnet(t(X_tr), NULL, 0, lambda_best)   
save(SenOCLR_final,
     file = file.path(outdir, paste0("Step6_SenOCLR_L2_",lambda_best,"_final.rdata")))

############################################################
## Step 7 – compare all quantification methods
############################################################

l2_tag_best <- paste0("l2_", lambda_best)

TISIs_list_best <- TISIs_list[grepl(paste0("^", l2_tag_best, ":"), names(TISIs_list))]

TISIs_list_best2 <- lapply(TISIs_list_best, function(x){
  y <- x$mean        # list of 500
  z <- x$gsva        # list of 4（100/200/500/1000）

  x$mean <- NULL
  x$gsva <- NULL

  ## 其余 spearman / pearson / dot / logistic 都是 list[[1]]
  x <- lapply(x, function(h) h[[1]])

  n <- names(x)
  x <- c(x, y, z)
  names(x) <- c(n, paste0("mean:", names(y)), paste0("gsva:", names(z)))
  return(x)
})


auc_models_final <- list()

for (q2 in names(TISIs_list_best2[[1]])) {

  TISIs_q2 <- lapply(TISIs_list_best2, function(x) x[[q2]])
  auc_celltypes <- list()

  for (nn in names(TISIs_q2)) {

    pred <- TISIs_q2[[nn]]
    samples <- names(pred)

    test_label <- meta[samples, "label"]
    cancer_type_sam <- meta[samples, "Cancer_type"]

    for (ct in unique(cancer_type_sam)) {

      idx <- which(cancer_type_sam == ct)
      pred_c <- pred[idx]
      test_label_c <- test_label[idx]

      if (length(unique(test_label_c)) < 2) next

      auc_val <- pROC::auc(
        pROC::roc(test_label_c, pred_c,
                  levels=c("non-senescent","senescent"), direction="<")
      )

      if (!ct %in% names(auc_celltypes))
        auc_celltypes[[ct]] <- as.numeric(auc_val)
      else
        auc_celltypes[[ct]] <- c(auc_celltypes[[ct]], as.numeric(auc_val))
    }
  }

  auc_models_final[[q2]] <- setNames(
    list(auc_celltypes),
    paste0("l2_", lambda_best)
  )
}

auc_final_df <- do.call(rbind,
  lapply(names(auc_models_final), function(q2){
    do.call(rbind,
      lapply(names(auc_models_final[[q2]]), function(l2){
        do.call(rbind,
          lapply(names(auc_models_final[[q2]][[l2]]), function(ct){
            data.frame(
              quant_method = q2,
              L2           = lambda_best,
              Cancer_type  = ct,
              AUC          = auc_models_final[[q2]][[l2]][[ct]]
            )
          })
        )
      })
    )
  })
)

write.csv(auc_final_df,
          file.path(outdir, "Step7_finalQuantMethods_AUC_detail.csv"),
          row.names=FALSE)

auc_final_summary <- auc_final_df %>%
  group_by(quant_method) %>%
  summarise(meanAUC=mean(AUC), sdAUC=sd(AUC), .groups="drop")

write.csv(
  auc_final_summary,
  file.path(outdir,"Step7_finalQuantMethods_AUC_summary.csv"),
  row.names=FALSE
)

############################################################
##  Step 8 Compute TISI scores
############################################################

TISI_scores_list <- scoreOCLR(profile = X_norm, model = SenOCLR_final, m = method_best)

TISI_raw <- TISI_scores_list[[1]]
TISI_scaled <- minmax(TISI_raw)

TISI_df <- data.frame(
  Sample      = colnames(X_norm),
  TISI_raw    = as.numeric(TISI_raw[colnames(X_norm)]),
  TISI_scaled = as.numeric(TISI_scaled[colnames(X_norm)]),
  Group       = meta$Group,
  label       = meta$label,
  Cancer_type = meta$Cancer_type,
  Cell_line   = meta$Cell_line,
  agent       = meta$agent,
  row.names   = colnames(X_norm)
)

write.csv(TISI_df,
          file.path(outdir,"Step8_TISI_scores_all_samples.csv"),
          row.names=FALSE)

