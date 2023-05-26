

# todo do this properly
source("/Users/alandenadel/Code/repos/PCKnockoffs/R/estimate_zipoisson.R")


#' @title Returns a Seurat object that contains additional (fake) RNA expression counts in the form of knockoffs.
#'
#' @description Given a Seurat object, returns a new Seurat object whose RNA expression counts includes the 
#' variable features from the original object and an equal number of knockoff features.
#'
#' @details 
#'
#' @param seurat_obj A Seurat object containing RNA expression counts
#' @returns A Seurat object that contains the original variable features and an equal number of knockoff features.
#' @examples
#' @name get_seurat_obj_with_knockoffs
#' @export
get_seurat_obj_with_knockoffs <- function(seurat_obj, subset_n=10000) {
  var.features <- Seurat::VariableFeatures(seurat_obj)
  #seurat_obj_data <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@counts)))
  
  #seurat_obj_data <- seurat_obj_data[var.features]

  print("Pulling data from Seurat object")
  seurat_obj_data <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@counts[VariableFeatures(seurat_obj),])))

  print("Subsetting cells")
  if (subset_n < nrow(seurat_obj_data)) {
    seurat_obj_data <- seurat_obj_data[sample(nrow(seurat_obj_data), subset_n), ]
  }

  print(dim(seurat_obj_data))
  
  
  print("Computing MLE for zero inflated poisson")
  ml_estimates <- lapply(seurat_obj_data, estimate_zi_poisson)
  
  print("computing knockoffs")
  ko <- as.data.frame(lapply(ml_estimates, function(x) rzipoisson(nrow(seurat_obj_data), x$lambda.hat, x$pi.hat)))
  

  num_variable_features <- length(var.features)
  colnames(ko) <- paste0(rep('knockoff', num_variable_features), 1:num_variable_features)
  combined.data <- cbind(seurat_obj_data, ko)
  
  new_project_name <- paste0(seurat_obj@project.name, "_with_knockoffs")
  new_seurat_obj <- Seurat::CreateSeuratObject(counts = t(combined.data), project = new_project_name)
  
  return(new_seurat_obj)
}




#' @name cluster_optimal_louvain_resolution_parameter
#' @export
cluster_optimal_louvain_resolution_parameter <- function(seurat_obj, original_num_clusters, num_variable_features) {

  # todo make this more efficient
  # todo move this function into seurat_workflow so that it happens automatically and doesn't repeat computation

  new_seurat_obj <- seurat_obj

  resolution_params <- seq(0.1, 2, 0.1)


  new_seurat_obj <- Seurat::NormalizeData(new_seurat_obj)
 
  new_seurat_obj <- Seurat::FindVariableFeatures(new_seurat_obj, selection.method = "vst", nfeatures = num_variable_features)
  
  all.genes <- rownames(new_seurat_obj)
  
  new_seurat_obj <- Seurat::ScaleData(new_seurat_obj, features = all.genes)
  
  new_seurat_obj <- Seurat::RunPCA(new_seurat_obj, features = VariableFeatures(object = new_seurat_obj))
  new_seurat_obj <- Seurat::FindNeighbors(new_seurat_obj, dims = 1:10) # todo check num dims

  for (resolution_param in resolution_params) {
    new_seurat_obj <- Seurat::FindClusters(new_seurat_obj, resolution = resolution_param)

    print("Original num clusters")
    print(original_num_clusters)

    print("Current clusters")
    print(length(levels(Idents(new_seurat_obj))))

    if (length(levels(Idents(new_seurat_obj))) > original_num_clusters) {
      stop("No resolution value possible")
    }

    if (length(levels(Idents(new_seurat_obj))) == original_num_clusters) {
      print(paste0("Resolution param:", resolution_param))
      return(resolution_param)
    }

  }
}


#' @title Runs a typical Seurat workflow on a Seurat object (up to dimensionality reduction and clustering).
#'
#' @description Given a Seurat object, returns a new Seurat that has been normalized, had variable features identified,
#' scaled, had principal components computed, had clusters identified, and had tSNE and UMAP embeddings determined.
#'
#' @details 
#'
#' @param seurat_obj A Seurat object that will be analyzed.
#' @param num_variable_features The number of variable features to use in the analysis.
#' @param resolution_param The resolution parameter to use when clustering.
#' @returns A Seurat object containing the relevant analysis results.
#' @examples
#' @name seurat_workflow
#' @export
seurat_workflow <- function(seurat_obj, num_variable_features, resolution_param=0.5) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
 
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num_variable_features)
  
  all.genes <- rownames(seurat_obj)
  
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all.genes)
  
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:10) # todo check if i should use all dims for knockoffs
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution_param)
  
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)
  seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = 1:10)
  
  # todo differential expression
  
  
  
  return(seurat_obj)
}


#' @title Returns the adjusted Rand index for the active identities in two Seurat objects
#'
#' @description Given two Seurat objects, returns the adjusted Rand index for the active identities of each of the objects.
#' @details 
#'
#' @param seurat_obj1 A Seurat object 
#' @param seurat_obj2 A Seurat object
#' @returns The adjusted Rand index for the active identities of each of the objects.
#' @examples
#' compare_Seurat_clusterings(seurat_obj1, seurat_obj2)
#' @name compare_Seurat_clusterings
#' @export
compare_Seurat_clusterings <- function(seurat_obj1, seurat_obj2) {
  # todo plot umap/tsne side by side
  
  # compute ARI
  cluster1 <- seurat_obj1@active.ident
  cluster2 <- seurat_obj2@active.ident
  
  ari <- pdfCluster::adj.rand.index(cluster1, cluster2)
  return(ari)
}

#' @name compute_knockoff_filter_one_cluster
#' @export
compute_knockoff_filter_one_cluster <- function(seurat_obj, cluster, q) {
  # todo test this function

  markers <- Seurat::FindMarkers(seurat_obj,
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
  
  thres = knockoff::knockoff.threshold(W, fdr=q, offset=1)
  
  print(paste("threshold", thres))
  
  selected_indices = which(W >= thres)
  
  selected_genes <- original_names_sorted[selected_indices]
  selected_Ws <- W[selected_indices]
  
  ret <- as.data.frame(list("selected_gene" = selected_genes, "W" = selected_Ws))
  
  ret <- ret[order(ret$W, decreasing = TRUE),]
  
  return(ret)
  
}

#' @title Returns the genes selected by the knockoff filter
#'
#' @description Given two Seurat objects, returns the  the genes selected by the knockoff filter and their W statistics.
#' @details 
#'
#' @param seurat_obj1 A Seurat object 
#' @param seurat_obj2 A Seurat object
#' @param cluster1 The Idents of the cluster of interest in seurat_obj1
#' @param cluster2 The Idents of the cluster of interest in seurat_obj2
#' @param q The desired rate to control the FDR at
#' @param return_all Determines if the returned object will contain all genes or just the selected genes.

#' @returns The adjusted Rand index for the active identities of each of the objects.
#' @examples
#' compute_knockoff_filter(seurat_obj1, seurat_obj2, 1, 2, 0.05)
#' @name compute_knockoff_filter
#' @export
compute_knockoff_filter <- function(seurat_obj, cluster1, cluster2, q, return_all=FALSE) {
  markers <- Seurat::FindMarkers(seurat_obj,
                         ident.1 = cluster1,
                         ident.2 = cluster2,
                         logfc.threshold = 0,
                         min.pct = 0)
  

  # FindMarkers orders by p-value, so we can't rely on position to know which genes are which
  knockoff.indices <- grepl("^knockoff", rownames(markers))
  original.indices <- !knockoff.indices
  
  # subset the markers data.frame into originals and knockoffs
  knockoff.markers <- markers[knockoff.indices, ]
  original.markers <- markers[original.indices, ]
  
  all.genes <- rownames(seurat_obj)
  
  # get indices of knockoffs and originals from seurat_obj, should be [FALSE, ..., FALSE, TRUE, ..., TRUE]
  knockoff.indices.sorted <- grepl("^knockoff", all.genes)
  original.indices.sorted <- !knockoff.indices.sorted
  
  knockoff_names_sorted <- all.genes[knockoff.indices.sorted]
  original_names_sorted <- all.genes[original.indices.sorted]

  # sort markers data.frames by their original orderings
  knockoff.markers.sorted <- knockoff.markers[knockoff_names_sorted, ]
  original.markers.sorted <- original.markers[original_names_sorted, ]
  
  original_p_values <- original.markers.sorted$p_val
  knockoff_p_values <- knockoff.markers.sorted$p_val
  
  log_original_p_values <- -log10(original_p_values)
  log_knockoff_p_values <- -log10(knockoff_p_values)
  
  W <- log_original_p_values - log_knockoff_p_values
  
  thres = knockoff::knockoff.threshold(W, fdr=q, offset=1)
  
  hist(W, breaks = 30)
  abline(v=thres)
  abline(v=-thres)

  W_hist <- ggplot2::ggplot() + aes(W) + 
    geom_histogram(binwidth=1, colour="black", fill="white") + 
    geom_vline(xintercept=thres) +
    geom_vline(xintercept=-thres) +
    theme(axis.line = element_line(colour = "black", size=2),
        panel.background = element_blank(),
        axis.title.x = element_text("W", size = 8), 
        axis.title.y = element_text("Count", size = 8), 
        #axis.text.x = element_blank(), 
        #axis.text.y = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_blank(),
        legend.position = "none")


  ggsave("B_cell_W_hist.png")

  
  print(paste("threshold", thres))

  if (return_all) {
    ret <- as.data.frame(list("gene" = original_names_sorted, "W" = W))

    return(ret)
  }
  selected_indices = which(W >= thres)
  
  selected_genes <- original_names_sorted[selected_indices]
  selected_Ws <- W[selected_indices]
  
  ret <- as.data.frame(list("selected_gene" = selected_genes, "W" = selected_Ws))
  
  ret <- ret[order(ret$W, decreasing = TRUE),]
  
  return(ret)
}






#' @name FindMarkersWithKnockoffs
#' @export
FindMarkersWithKnockoffs <- function(seurat_obj, ident.1, ident.2, q, num_var_features) {
  # todo test this function
  # todo create knockoffs from entire dataset before subsetting 
  # todo subset dataset with knockoffs based on the chosen cluster identities
  # todo redo clustering with K=2
  # todo test markers between knew clusters
  new_seurat_obj <- subset(seurat_obj, idents = c(ident.1, ident.2))
  
  # remove genes that are lowly expressed (i.e. all zeroes)
  new_seurat_obj <- Seurat::CreateSeuratObject(new_seurat_obj@assays$RNA@counts, min.cells = 10)
  

  new_seurat_obj <- Seurat::FindVariableFeatures(new_seurat_obj, selection.method = "vst", nfeatures = num_var_features)
  
  
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






#' @title Returns the Jaccard coefficient of two vectors.
#'
#' @description Given two vectors, returns their Jaccard coefficient.
#' @details 
#'
#' @param a A vector 
#' @param b A vector
#' @returns The Jaccard coefficient of two vectors.
#' @examples
#' jaccard(c(1,2,3), c(3,4,5))
#' @name jaccard
#' @export
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#' @name compare_markers_jaccard
#' @export
compare_markers_jaccard <- function(orig_seurat_obj, knock_seurat_obj, q,
                                    orig_ident_1, orig_ident_2, 
                                    knock_ident_1, knock_ident_2) {
  
  original.markers <- Seurat::FindMarkers(orig_seurat_obj, ident.1 = orig_ident_1, ident.2 = orig_ident_2, features = VariableFeatures(orig_seurat_obj))
  #knockoff.markers <- FindMarkers(knock_seurat_obj, ident.1 = knock_ident_1, ident.2 = knock_ident_2)
  
  markers.selected <- compute_knockoff_filter(knock_seurat_obj, knock_ident_1, knock_ident_2, q)
  
  top.original <- rownames(head(original.markers, 100))
  top.knockoffs <- head(markers.selected, 100)$selected_gene
  
  
  ret_list <- list("jaccard" = jaccard(top.original, top.knockoffs),
                   "original_selected" = dim(original.markers)[1],
                   "knockoff_selected" = dim(markers.selected)[1])
  return(ret_list)
}