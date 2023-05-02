library(Seurat)
library(pdfCluster) # for adjusted rand index
library(fossil) # for adjusted rand index and rand index

library(dplyr)

library(knockoff)

source("estimate_zipoisson.R")



get_seurat_obj_with_knockoffs <- function(seurat_obj) {
  var.features <- VariableFeatures(seurat_obj)
  seurat_obj_data <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@counts)))
  
  seurat_obj_data <- seurat_obj_data[var.features]
  
  
  ml_estimates <- lapply(seurat_obj_data, estimate_zi_poisson)
  
  print("computing knockoffs")
  ko <- as.data.frame(lapply(ml_estimates, function(x) rzipoisson(dim(seurat_obj_data[1]), x$lambda.hat, x$pi.hat)))
  

  num_variable_features <- length(var.features)
  colnames(ko) <- paste0(rep('knockoff', num_variable_features), 1:num_variable_features)
  combined.data <- cbind(seurat_obj_data, ko)
  
  new_project_name <- paste0(seurat_obj@project.name, "_with_knockoffs")
  new_seurat_obj <- CreateSeuratObject(counts = t(combined.data), project = new_project_name)
  
  return(new_seurat_obj)
}


seurat_workflow <- function(seurat_obj, num_variable_features) {
  seurat_obj <- NormalizeData(seurat_obj)
 
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num_variable_features)
  
  all.genes <- rownames(seurat_obj)
  
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10) # todo check if i should use all dims for knockoffs
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)
  
  # todo differential expression
  
  
  
  return(seurat_obj)
}

compare_clusterings <- function(seurat_obj1, seurat_obj2) {
  # todo plot umap/tsne side by side
  
  # compute ARI
  cluster1 <- seurat_obj1@active.ident
  cluster2 <- seurat_obj2@active.ident
  
  ari <- adj.rand.index(cluster1, cluster2)
  return(ari)
}

compute_knockoff_filter_one_cluster <- function(seurat_obj, cluster, q) {
  markers <- FindMarkers(seurat_obj,
                         ident.1 = cluster,
                         logfc.threshold = 0,
                         min.pct = 0)
  
  knockoff.indices <- grepl("^knockoff", rownames(markers))
  original.indices <- !knockoff.indices
  
  knockoff.markers <- markers[knockoff.indices, ]
  original.markers <- markers[original.indices, ]
  
  all.genes <- rownames(seurat_obj)
  
  knockoff.indices.sorted <- grepl("^knockoff", all.genes)
  original.indices.sorted <- !knockoff.indices.sorted
  
  knockoff_names_sorted <- all.genes[knockoff.indices.sorted]
  original_names_sorted <- all.genes[original.indices.sorted]
  
  knockoff.markers.sorted <- knockoff.markers[knockoff_names_sorted, ]
  original.markers.sorted <- original.markers[original_names_sorted, ]
  
  original_p_values <- original.markers.sorted$p_val
  knockoff_p_values <- knockoff.markers.sorted$p_val
  
  log_original_p_values <- -log10(original_p_values)
  log_knockoff_p_values <- -log10(knockoff_p_values)
  
  W <- log_original_p_values - log_knockoff_p_values
  
  thres = knockoff.threshold(W, fdr=q, offset=1)
  
  print(paste("threshold", thres))
  
  selected_indices = which(W >= thres)
  
  selected_genes <- original_names_sorted[selected_indices]
  selected_Ws <- W[selected_indices]
  
  ret <- as.data.frame(list("selected_gene" = selected_genes, "W" = selected_Ws))
  
  ret <- ret[order(ret$W, decreasing = TRUE),]
  
  return(ret)
  
}


compute_knockoff_filter <- function(seurat_obj, cluster1, cluster2, q) {
  markers <- FindMarkers(seurat_obj,
                         ident.1 = cluster1,
                         ident.2 = cluster2,
                         logfc.threshold = 0,
                         min.pct = 0)
  

  knockoff.indices <- grepl("^knockoff", rownames(markers))
  original.indices <- !knockoff.indices
  
  knockoff.markers <- markers[knockoff.indices, ]
  original.markers <- markers[original.indices, ]
  
  all.genes <- rownames(seurat_obj)
  
  knockoff.indices.sorted <- grepl("^knockoff", all.genes)
  original.indices.sorted <- !knockoff.indices.sorted
  
  knockoff_names_sorted <- all.genes[knockoff.indices.sorted]
  original_names_sorted <- all.genes[original.indices.sorted]

  knockoff.markers.sorted <- knockoff.markers[knockoff_names_sorted, ]
  original.markers.sorted <- original.markers[original_names_sorted, ]
  
  original_p_values <- original.markers.sorted$p_val
  knockoff_p_values <- knockoff.markers.sorted$p_val
  
  log_original_p_values <- -log10(original_p_values)
  log_knockoff_p_values <- -log10(knockoff_p_values)
  
  W <- log_original_p_values - log_knockoff_p_values
  
  thres = knockoff.threshold(W, fdr=q, offset=1)
  
  hist(W, breaks = 100)
  abline(v=thres)
  abline(v=-thres)
  
  print(paste("threshold", thres))

  
  selected_indices = which(W >= thres)
  
  selected_genes <- original_names_sorted[selected_indices]
  selected_Ws <- W[selected_indices]
  
  ret <- as.data.frame(list("selected_gene" = selected_genes, "W" = selected_Ws))
  
  ret <- ret[order(ret$W, decreasing = TRUE),]
  
  return(ret)
}





# todo test this function
# todo create knockoffs from entire dataset before subsetting 
# todo subset dataset with knockoffs based on the chosen cluster identities
# todo redo clustering with K=2
# todo test markers between knew clusters
FindMarkersWithKnockoffs <- function(seurat_obj, ident.1, ident.2, q, num_var_features) {
  new_seurat_obj <- subset(seurat_obj, idents = c(ident.1, ident.2))
  
  # remove genes that are lowly expressed (i.e. all zeroes)
  new_seurat_obj <- CreateSeuratObject(new_seurat_obj@assays$RNA@counts, min.cells = 10)
  

  new_seurat_obj <- FindVariableFeatures(new_seurat_obj, selection.method = "vst", nfeatures = num_var_features)
  
  
  new_seurat_obj <- get_seurat_obj_with_knockoffs(new_seurat_obj)
  
  

  x <- t(as.matrix(new_seurat_obj@assays$RNA@counts))

  kmeans_fit <- kmeans(x, 2)
  kmeans_idents <- as.factor(kmeans_fit$cluster)
  
  Idents(new_seurat_obj) <- kmeans_idents
  
  cluster_labels <- levels(kmeans_idents)
  
  print(new_seurat_obj)
  
  markers <- compute_knockoff_filter(new_seurat_obj, cluster_labels[1], cluster_labels[2], q)
  
  return(markers)
}






# thres = knockoff.threshold(W, fdr=q, offset=1)

# todo implement this k-fwer threshold and the more conservative threshold with the min
# https://arxiv.org/abs/1505.06549
# Familywise error rate control via knockoffs
# https://github.com/amspector100/knockpy/blob/master/knockpy/knockoff_stats.py#L1324
# https://github.com/cran/knockoff/blob/master/R/knockoff_filter.R#L181
kfwer_knockoff_threshold <- function(W, k, fwer) {
  # todo report if there is no valid v parameter (i.e. if k is small)
  
}



optimize_louvain_resolution_parameter <- function(seurat_obj, original_num_clusters) {
}
