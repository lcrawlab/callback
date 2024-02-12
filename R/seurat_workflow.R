

#' @title Returns a Seurat object that contains additional (fake) RNA expression counts in the form of knockoffs.
#'
#' @description Given a Seurat object, returns a new Seurat object whose RNA expression counts includes the 
#' variable features from the original object and an equal number of knockoff features.
#'
#' @details 
#'
#' @param seurat_obj A Seurat object containing RNA expression counts.
#' @param assay The assay to generate knockoffs from.
#' @returns A Seurat object that contains the original variable features and an equal number of knockoff features.
#' @name get_seurat_obj_with_knockoffs
get_seurat_obj_with_knockoffs <- function(seurat_obj, assay = "RNA") {
  var.features <- Seurat::VariableFeatures(seurat_obj)
  #seurat_obj_data <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@counts)))
  
  #seurat_obj_data <- seurat_obj_data[var.features]

  message("Pulling data from Seurat object")
  #seurat_obj_data <- as.data.frame(t(as.matrix(seurat_obj@assays$RNA@counts[Seurat::VariableFeatures(seurat_obj),])))
  seurat_obj_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seurat_obj, assay=assay, slot='counts')[Seurat::VariableFeatures(seurat_obj),])))

  
  
  message("Computing MLE for zero-inflated poisson")
  ml_estimates <- lapply(seurat_obj_data, estimate_zi_poisson)
  
  message("Computing knockoffs")
  ko <- as.data.frame(lapply(ml_estimates, function(x) rzipoisson(nrow(seurat_obj_data), x$lambda.hat, x$pi.hat)))
  

  num_variable_features <- length(var.features)
  colnames(ko) <- paste0(rep('knockoff', num_variable_features), 1:num_variable_features)
  combined.data <- cbind(seurat_obj_data, ko)
  
  new_project_name <- paste0(seurat_obj@project.name, "_with_knockoffs")
  new_seurat_obj <- Seurat::CreateSeuratObject(counts = t(combined.data), project = new_project_name)
  
  return(new_seurat_obj)
}




#' @name cluster_optimal_louvain_resolution_parameter
cluster_optimal_louvain_resolution_parameter <- function(seurat_obj,
                                                         original_num_clusters,
                                                         num_variable_features,
                                                         res_start = 0.1,
                                                         res_end = 3,
                                                         res_increment = 0.05) {

  # todo make this more efficient
  # todo move this function into  flow so that it happens automatically and doesn't repeat computation

  new_seurat_obj <- seurat_obj

  resolution_params <- seq(res_start, res_end, res_increment)


  new_seurat_obj <- Seurat::NormalizeData(new_seurat_obj)
 
  new_seurat_obj <- Seurat::FindVariableFeatures(new_seurat_obj, selection.method = "vst", nfeatures = num_variable_features)
  
  all.genes <- rownames(new_seurat_obj)
  
  new_seurat_obj <- Seurat::ScaleData(new_seurat_obj, features = all.genes)
  
  new_seurat_obj <- Seurat::RunPCA(new_seurat_obj, features = Seurat::VariableFeatures(object = new_seurat_obj))
  new_seurat_obj <- Seurat::FindNeighbors(new_seurat_obj, dims = 1:10) # todo check num dims

  for (resolution_param in resolution_params) {
    new_seurat_obj <- Seurat::FindClusters(new_seurat_obj, resolution = resolution_param)

    print("Original num clusters")
    print(original_num_clusters)

    print("Current clusters")
    print(length(levels(Seurat::Idents(new_seurat_obj))))

    if (length(levels(Seurat::Idents(new_seurat_obj))) > original_num_clusters) {
      stop("No resolution value possible")
    }

    if (length(levels(Seurat::Idents(new_seurat_obj))) == original_num_clusters) {
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
#' @param visualization_method Either "umap" or "tsne".
#' @param num_dims The number of principal components to use.
#' @param algorithm The clustering algorithm to use, either "louvain" or "leiden".
#' @returns A Seurat object containing the relevant analysis results.
#' @export
#' @name seurat_workflow
seurat_workflow <- function(seurat_obj,
                            num_variable_features,
                            resolution_param = 0.5,
                            visualization_method = "umap",
                            num_dims = 10,
                            algorithm = "louvain") {
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
   
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num_variable_features)
    
  all.genes <- rownames(seurat_obj)
    
  #seurat_obj <- Seurat::ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
    
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(object = seurat_obj))
    
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:num_dims) # todo check if i should use all dims for knockoffs

  if (algorithm == "louvain") {
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution_param)
  }

  if (algorithm == "leiden") {
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution_param, algorithm=4, method = "igraph")
  }

  if (visualization_method == "umap") {
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:num_dims)
  }
  if (visualization_method == "tsne") {
    seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = 1:num_dims)
  }

if (visualization_method == "both") {
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:num_dims)
    seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = 1:num_dims)
  } 
  
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
#' compare_Seurat_clusterings(seurat_obj1, seurat_obj2)
#' @name compare_Seurat_clusterings
compare_Seurat_clusterings <- function(seurat_obj1, seurat_obj2) {
  # todo plot umap/tsne side by side
  
  # compute ARI
  cluster1 <- seurat_obj1@active.ident
  cluster2 <- seurat_obj2@active.ident
  
  ari <- pdfCluster::adj.rand.index(cluster1, cluster2)
  return(ari)
}

#' @name compute_knockoff_filter_one_cluster
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

  #selected_indices = which(W >= thres)
  selected_indices = which(W > thres)

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
#' @param seurat_obj A Seurat object 
#' @param cluster1 The Idents of the cluster of interest in seurat_obj1
#' @param cluster2 The Idents of the cluster of interest in seurat_obj2
#' @param q The desired rate to control the FDR at
#' @param return_all Determines if the returned object will contain all genes or just the selected genes.
#' @param threshold One of "fdr", "kfwer", or "heuristic".
#' @param num_cores The number of cores for computing marker genes in parallel.
#' @returns The 
#' @name compute_knockoff_filter
compute_knockoff_filter <- function(seurat_obj, cluster1, cluster2, q, return_all=FALSE, threshold="fdr", num_cores=1) {
  #library(future)
  options(future.globals.maxSize = 8000 * 1024^2) # todo note what this is for, figure this out as a parameter or programmatically
  future::plan("multicore", workers = as.numeric(num_cores)) # todo log number of cores being used
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

  if (threshold == "fdr") {
    thres = knockoff::knockoff.threshold(W, fdr=q, offset=1)
  }

  if (threshold == "kfwer") {
    k <- 10
    alpha <- q
    thres <- knockoff.kfwer.threshold(W, k, alpha)
  }

  if (threshold == "heuristic") {
    k <- 10
    alpha <- q
    thres <- knockoff.heuristic.threshold(W, q, k, alpha)
  }

  #print(paste("threshold:", thres))

  if (return_all) {
    all_features <- as.data.frame(list("gene" = original_names_sorted, "W" = W))

    ret <-  list("all_features"=all_features, "threshold"=thres)

    return(ret)
  }
  selected_indices = which(W >= thres) # todo check if this should be > (case where threshold is Inf, but there are still some Inf -log p)
  #selected_indices = which(W > thres) # todo check if this should be > (case where threshold is Inf, but there are still some Inf -log p)
  

  #print("Num selected indices:")
  #print(length(selected_indices))
  selected_genes <- original_names_sorted[selected_indices]
  selected_Ws <- W[selected_indices]
  
  selected_features <- as.data.frame(list("selected_gene" = selected_genes, "W" = selected_Ws))
  
  selected_features <- selected_features[order(selected_features$W, decreasing = TRUE),]

  ret <-  list("selected_features"=selected_features, "threshold"=thres)
  
  return(ret)
}






#' @name FindMarkersWithKnockoffs
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

  kmeans_fit <- stats::kmeans(x, 2)
  kmeans_idents <- as.factor(kmeans_fit$cluster)

  Seurat::Idents(new_seurat_obj) <- kmeans_idents

  cluster_labels <- levels(kmeans_idents)

  print(new_seurat_obj)

  markers <- compute_knockoff_filter(new_seurat_obj, cluster_labels[1], cluster_labels[2], q)

  return(markers)
}






#' @title Returns the Jaccard coefficient of two vectors.
#'
#' @description Given two vectors, returns their Jaccard coefficient.
#' @details blank
#'
#' @param a A vector 
#' @param b A vector
#' @returns The Jaccard coefficient of two vectors.
#' @name jaccard
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#' @name compare_markers_jaccard
compare_markers_jaccard <- function(orig_seurat_obj, knock_seurat_obj, q,
                                    orig_ident_1, orig_ident_2, 
                                    knock_ident_1, knock_ident_2) {
  
  original.markers <- Seurat::FindMarkers(orig_seurat_obj,
                                          ident.1 = orig_ident_1,
                                          ident.2 = orig_ident_2, 
                                          features = Seurat::VariableFeatures(orig_seurat_obj),
                                          logfc.threshold = 0,
                                          min.pct = 0)
  #knockoff.markers <- FindMarkers(knock_seurat_obj, ident.1 = knock_ident_1, ident.2 = knock_ident_2)
  
  markers.selected <- compute_knockoff_filter(knock_seurat_obj, knock_ident_1, knock_ident_2, q)

  thres <- markers.selected$threshold
  markers.selected <- markers.selected$selected_features

  top.original <- rownames(utils::head(original.markers, 50))
  top.knockoffs <- utils::head(markers.selected, 50)$selected_gene

  ret_list <- list("jaccard" = jaccard(top.original, top.knockoffs),
                   "original_selected" = dim(original.markers)[1],
                   "knockoff_selected" = dim(markers.selected)[1],
                   "intersection" = length(intersect(top.original, top.knockoffs)),
                   "top_original" = top.original,
                   "top_knockoffs" = top.knockoffs
                   )
  return(ret_list)
}



cluster_seurat_heuristically_find_num_clusters_hierarchical_clustering <- function(seurat_obj,
                                                                                   starting_num_clusters,
                                                                                   num_dims=10,
                                                                                   cores=1) {
  knockoff_seurat_obj <- get_seurat_obj_with_knockoffs(seurat_obj)

  num_variable_features <- 2 * length(Seurat::VariableFeatures(seurat_obj))


  # Pre-process data
  knockoff_seurat_obj <- Seurat::NormalizeData(knockoff_seurat_obj)
   
  knockoff_seurat_obj <- Seurat::FindVariableFeatures(knockoff_seurat_obj, selection.method = "vst", nfeatures = num_variable_features)
    
  all.genes <- rownames(knockoff_seurat_obj)
    
  knockoff_seurat_obj <- Seurat::ScaleData(knockoff_seurat_obj)
    
  knockoff_seurat_obj <- Seurat::RunPCA(knockoff_seurat_obj, features = Seurat::VariableFeatures(object = knockoff_seurat_obj))
    
  #knockoff_seurat_obj <- Seurat::FindNeighbors(knockoff_seurat_obj, dims = 1:num_dims) # todo check if i should use all dims for knockoffs

  hierarchical_clust_results <- stats::hclust(stats::dist(knockoff_seurat_obj@reductions$pca@cell.embeddings[,1:num_dims]), method = "ward.D")

  num_clusters <- starting_num_clusters + 1
  while(TRUE) {
    num_clusters <- num_clusters - 1
    print("####################################################################")
    print("Num Clusters:")
    print(num_clusters)
    print("####################################################################")

    k <- num_clusters

    knock_idents <- 1:k

    num_selected_matrix <- matrix(nrow=k, ncol=k)

    found_no_sign_diff <- FALSE

    Seurat::Idents(knockoff_seurat_obj) <- stats::cutree(hierarchical_clust_results, k)
    
    m <- 0
    for (i in 1:length(knock_idents)) {
      for (j in 1:length(knock_idents)) {
        if (j >= i) {
          next
        }

        m <- m + 1

        print("Pair:")
        print(paste(i,j))

        print("Knockoff Pair:")
        print(paste(knock_idents[i], knock_idents[j]))

        markers.selected <- compute_knockoff_filter(knockoff_seurat_obj, knock_idents[i], knock_idents[j], 0.05, num_cores=cores)

        num.selected <- nrow(markers.selected$selected_features)

        if (num.selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }

        num_selected_matrix[i, j] <- num.selected
        num_selected_matrix[j, i] <- num.selected

      }
      if (found_no_sign_diff) {
        print("Found clusters with no significant differences. Progressing to next clustering iteration.")
        break
      }
    }

    #num.zero.findings <- sum(num_selected_matrix == 0, na.rm = TRUE)

    #print(num.zero.findings)

    #print(num_selected_matrix)

    #if (num.zero.findings > 0) {
    if (found_no_sign_diff) {
      next
    }

    break


  }

    ret <-  list("knockoff_seurat_obj" = knockoff_seurat_obj,
                 "DEG_results" = "",
                 "num_selected_matrix" = num_selected_matrix)

  return(ret)

}


cluster_seurat_heuristically_find_num_clusters_k_means <- function(seurat_obj,
                                                                   starting_num_clusters,
                                                                   num_dims = 10,
                                                                   cores = 1) {
  knockoff_seurat_obj <- get_seurat_obj_with_knockoffs(seurat_obj)

  num_variable_features <- 2 * length(Seurat::VariableFeatures(seurat_obj))


  # Pre-process data
  knockoff_seurat_obj <- Seurat::NormalizeData(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::FindVariableFeatures(knockoff_seurat_obj,
                                                      selection.method = "vst",
                                                      nfeatures = num_variable_features)

  all.genes <- rownames(knockoff_seurat_obj)
  knockoff_seurat_obj <- Seurat::ScaleData(knockoff_seurat_obj)
  knockoff_seurat_obj <- Seurat::RunPCA(knockoff_seurat_obj, 
                                        features = Seurat::VariableFeatures(object = knockoff_seurat_obj))
  # todo check if i should use all dims for knockoffs
  knockoff_seurat_obj <- Seurat::FindNeighbors(knockoff_seurat_obj,
                                               dims = 1:num_dims)

  num_clusters <- starting_num_clusters + 1
  while(TRUE) {
    num_clusters <- num_clusters - 1
    print("####################################################################")
    print("Num Clusters:")
    print(num_clusters)
    print("####################################################################")

    k <- num_clusters

    knock_idents <- 1:k

    num_selected_matrix <- matrix(nrow = k, ncol = k)

    found_no_sign_diff <- FALSE

    Seurat::Idents(knockoff_seurat_obj) <- stats::kmeans(knockoff_seurat_obj@reductions$pca@cell.embeddings[`1:num_dims], k, iter.max = 20, nstart = 10)

    m <- 0
    for (i in 1:length(knock_idents)) {
      for (j in 1:length(knock_idents)) {
        if (j >= i) {
          next
        }

        m <- m + 1

        print("Pair:")
        print(paste(i, j))

        print("Knockoff Pair:")
        print(paste(knock_idents[i], knock_idents[j]))

        markers.selected <- compute_knockoff_filter(knockoff_seurat_obj, knock_idents[i], knock_idents[j], 0.05, num_cores=cores)

        num.selected <- nrow(markers.selected$selected_features)

        if (num.selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }

        num_selected_matrix[i, j] <- num.selected
        num_selected_matrix[j, i] <- num.selected

      }
      if (found_no_sign_diff) {
        print("Found clusters with no significant differences. Progressing to next clustering iteration.")
        break
      }
    }
    
    #num.zero.findings <- sum(num_selected_matrix == 0, na.rm = TRUE)
    
    #print(num.zero.findings)

    #print(num_selected_matrix)

    #if (num.zero.findings > 0) {
    if (found_no_sign_diff) {
      next
    }

    break


  }

    ret <-  list("knockoff_seurat_obj" = knockoff_seurat_obj,
                 "DEG_results" = "",
                 "num_selected_matrix" = num_selected_matrix)

  return(ret)

}





cluster_seurat_heuristically_find_num_clusters <- function(seurat_obj,
                                                           resolution_param_start,
                                                           resolution_param_decrement,
                                                           visualization_method,
                                                           num_dims = 10,
                                                           algorithm = "louvain",
                                                           cores = 1) {
  
  # bump up starting point by the decrement because we subtract the decrement in the loop
  resolution_param <- resolution_param_start + resolution_param_decrement


  knockoff_seurat_obj <- get_seurat_obj_with_knockoffs(seurat_obj)

  num_variable_features <- 2 * length(Seurat::VariableFeatures(seurat_obj))


  # Pre-process data

  #library(future)
  options(future.globals.maxSize = 8000 * 1024^2)
  future::plan("multicore", workers = as.numeric(cores)) # todo log number of cores being used

  knockoff_seurat_obj <- Seurat::NormalizeData(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::FindVariableFeatures(knockoff_seurat_obj,
                                                      selection.method = "vst",
                                                      nfeatures = num_variable_features)

  all.genes <- rownames(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::ScaleData(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::RunPCA(knockoff_seurat_obj,
                                        features = Seurat::VariableFeatures(object = knockoff_seurat_obj))

  # todo check if i should use all dims for knockoffs
  knockoff_seurat_obj <- Seurat::FindNeighbors(knockoff_seurat_obj,
                                               dims = 1:num_dims)

  while(TRUE) {
    #resolution_param <- resolution_param - resolution_param_decrement

    resolution_param <- 0.8 * resolution_param

    print("####################################################################")
    print("Resolution param:")
    print(resolution_param)
    print("####################################################################")


    print("Finding clusters")

    if (algorithm == "louvain") {
      print("Louvain")
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj, resolution = resolution_param)
    }

    if (algorithm == "leiden") {
      print("Leiden")
      #plan("sequential") # todo log number of cores being used # this is a weird one because leiden has a forked job hanging
      knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj,
                                                  resolution = resolution_param,
                                                  algorithm = 4,
                                                  method = "igraph")
    }
        print("Found clusters")

    k <- length(levels(Seurat::Idents(knockoff_seurat_obj)))
    #knock_idents <- 0:(k-1)

    print("Num clusters:")
    print(k)

    knock_idents <- levels(Seurat::Idents(knockoff_seurat_obj))


    num_selected_matrix <- matrix(nrow=k, ncol=k)

    found_no_sign_diff <- FALSE

    m <- 0
    for (i in 1:length(knock_idents)) {
      for (j in 1:length(knock_idents)) {
        if (j >= i) {
          next
        }

        m <- m + 1

        print("Pair:")
        print(paste(i, j))

        print("Knockoff Pair:")
        print(paste(knock_idents[i], knock_idents[j]))

        markers.selected <- compute_knockoff_filter(knockoff_seurat_obj,
                                                    knock_idents[i],
                                                    knock_idents[j],
                                                    0.05,
                                                    num_cores = cores)

        num.selected <- nrow(markers.selected$selected_features)

        if (num.selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }

        num_selected_matrix[i, j] <- num.selected
        num_selected_matrix[j, i] <- num.selected

      }
      if (found_no_sign_diff) {
        print("Found clusters with no significant differences. Progressing to next clustering iteration.")
        break
      }
    }

    if (found_no_sign_diff) {
      next
    }

    break

  }

  ret <-  list("knockoff_seurat_obj" = knockoff_seurat_obj,
               "DEG_results" = "",
               "num_selected_matrix" = num_selected_matrix)

  return(ret)
}



get_cluster_centroid <- function(seurat_obj, num_PCs, cluster_ident) {
  smaller_obj <- subset(seurat_obj, idents = c(cluster_ident))

  data <- smaller_obj@reductions$pca@cell.embeddings[,1:num_PCs]

  centroid <- colSums(data)

  return(centroid)
}

get_cluster_centroid_distance_matrix <- function(seurat_obj, num_PCs) {
  cluster_names <- levels(Seurat::Idents(seurat_obj))

  centroids <- lapply(cluster_names, FUN = function(cluster_name) { get_cluster_centroid(seurat_obj, num_PCs, cluster_name) })

  centroid.df <- t(as.data.frame(centroids))
  rownames(centroid.df) <- cluster_names


  distance_mat <- stats::dist(centroid.df)

  return(as.matrix(distance_mat))
}


get_cluster_comparison_order_sorted_by_centroid_distance <- function(seurat_obj, num_PCs) {
  idents <- levels(Seurat::Idents(seurat_obj))
  comparisons <- utils::combn(idents, 2)
  num_comparisons <- ncol(comparisons)

  dist.mat <- get_cluster_centroid_distance_matrix(seurat_obj, num_PCs)

  unsorted_comparisons <- unlist(lapply(1:num_comparisons, function(comparison) { dist.mat[comparisons[1,comparison], comparisons[2,comparison]] }))
  names(unsorted_comparisons) <- 1:num_comparisons
  sorted_comparisons <- sort(unsorted_comparisons)

  new_order <- names(sorted_comparisons)

  return(new_order)
}

cluster_seurat_heuristically_find_num_clusters_sorted <- function(seurat_obj,
                                                           resolution_param_start,
                                                           resolution_param_decrement,
                                                           visualization_method,
                                                           num_dims=10) {
  
  # bump up starting point by the decrement because we subtract the decrement in the loop
  resolution_param <- resolution_param_start + resolution_param_decrement


  knockoff_seurat_obj <- get_seurat_obj_with_knockoffs(seurat_obj)

  num_variable_features <- 2 * length(Seurat::VariableFeatures(seurat_obj))


  # Pre-process data
  knockoff_seurat_obj <- Seurat::NormalizeData(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::FindVariableFeatures(knockoff_seurat_obj, selection.method = "vst", nfeatures = num_variable_features)

  all.genes <- rownames(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::ScaleData(knockoff_seurat_obj)

  knockoff_seurat_obj <- Seurat::RunPCA(knockoff_seurat_obj, features = Seurat::VariableFeatures(object = knockoff_seurat_obj))

  knockoff_seurat_obj <- Seurat::FindNeighbors(knockoff_seurat_obj, dims = 1:num_dims) # todo check if i should use all dims for knockoffs

  while(TRUE) {


    #resolution_param <- resolution_param - resolution_param_decrement

    resolution_param <- 0.8 * resolution_param

    print("####################################################################")
    print("Resolution param:")
    print(resolution_param)
    print("####################################################################")


    knockoff_seurat_obj <- Seurat::FindClusters(knockoff_seurat_obj, resolution = resolution_param)


    k <- length(levels(Seurat::Idents(knockoff_seurat_obj)))
    knock_idents <- 0:(k-1)


    num_selected_matrix <- matrix(nrow=k, ncol=k)

    found_no_sign_diff <- FALSE


    idents <- levels(Seurat::Idents(knockoff_seurat_obj))

    if (length(idents) == 1) {
      break
    }
    sorted_comparisons <- get_cluster_comparison_order_sorted_by_centroid_distance(knockoff_seurat_obj, num_dims)




    comparisons <- utils::combn(idents, 2)

    print("Idents")
    print(idents)

    print("comparisons")
    print(comparisons)
    m <- 0

    for (comparison in sorted_comparisons) {
      i <- as.numeric(comparisons[1, as.numeric(comparison)])
      j <- as.numeric(comparisons[2, as.numeric(comparison)])
      m <- m + 1

      print("Comparison:")
      print(m)

      print("Sorted Comparison:")
      print(comparison)

      print("Pair:")
      print(paste(i,j))



      print("Knockoff Pair:")
      print(paste(knock_idents[i], knock_idents[j]))

      #markers.selected <- compute_knockoff_filter(knockoff_seurat_obj, knock_idents[i], knock_idents[j], 0.05)
      markers.selected <- compute_knockoff_filter(knockoff_seurat_obj, i, j, 0.05)

      num.selected <- nrow(markers.selected$selected_features)

      if (num.selected == 0) {
        found_no_sign_diff <- TRUE
        print("Found clusters with no significant differences. Progressing to next clustering iteration.")

        break
      }

      num_selected_matrix[i, j] <- num.selected
      num_selected_matrix[j, i] <- num.selected
    }

    if (found_no_sign_diff) {
      next
    }
    break

  }

    ret <-  list("knockoff_seurat_obj"=knockoff_seurat_obj, "DEG_results"="", "num_selected_matrix"=num_selected_matrix)

  return(ret)
}